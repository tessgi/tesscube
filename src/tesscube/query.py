"""Mixin class to enable queries. This is in a separate class for readability."""

import asyncio
import struct
import warnings
from functools import lru_cache
from io import BytesIO
from typing import Optional, Tuple, List, Union

import numpy as np
from aiobotocore.session import get_session
from aiobotocore.client import AioBaseClient
from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning
from botocore import UNSIGNED
from botocore.config import Config

from . import (
    BUCKET_NAME,
    BYTES_PER_PIX,
    DATA_OFFSET,
    HDR_SIZE,
    MAX_CONCURRENT_DOWNLOADS,
)
from .fits import _fix_primary_hdu
from .utils import _sync_call, convert_coordinates_to_runs


class QueryMixin:
    def find_byte_offset(self, row: int, column: int, frame: int = 0) -> int:
        """
        Returns the byte offset of a specific pixel position.

        Parameters
        ----------
        row : int
            The row index of the pixel (1-indexed).
        column : int
            The column index of the pixel (1-indexed).
        frame : int, optional
            The frame index, by default 0.

        Returns
        -------
        offset: int
            The byte offset of the specified pixel.
        """
        # if (row < 1) | (row > 2049) | (column < 1) | (column > 2049):
        #     raise ValueError(f"`(row, column)` position `({row}, {column})` is outside of range.")
        # Subtract 1 from column and row because the byte location assumes zero-indexing,
        # whereas the TESS convention is to address column and row number with one-indexing.
        pixel_offset = (column - 1) * self.nframes * self.nsets + (
            row - 1
        ) * self.ncolumns * self.nframes * self.nsets
        pixel_offset += frame * self.nsets
        return DATA_OFFSET + BYTES_PER_PIX * pixel_offset

    def find_byte_range_subframes(
        self, row: int, column: int, frame_range: Optional[Tuple[int, int]] = None
    ) -> tuple:
        """
        Returns the byte range for a given row, column, and frame range.

        Parameters
        ----------
        row : int
            The row index of the pixel (1-indexed).
        column : int
            The column index of the pixel (1-indexed).
        frame_range : tuple of int, optional
            The range of frames to consider, by default None.

        Returns
        -------
        tuple of int
            The byte range (start, end) for the specified parameters.
        """
        if frame_range is None:
            frame_range = (0, self.nframes)
        begin = self.find_byte_offset(row=row, column=column, frame=frame_range[0])
        end = self.find_byte_offset(row=row, column=column, frame=frame_range[1])
        return (begin, end)

    def find_byte_range(self, row: int, column_range: tuple) -> tuple:
        """
        Returns the byte range for a given row and column range.

        Parameters
        ----------
        row : int
            The row index of the pixel (1-indexed).
        column_range : tuple of int
            The range of columns to consider.

        Returns
        -------
        tuple of int
            The byte range (start, end) for the specified parameters.
        """
        begin = self.find_byte_offset(row=row, column=column_range[0])
        end = self.find_byte_offset(row=row, column=column_range[1])
        return (begin, end)

    async def async_get_data_per_pixel(
        self,
        s3: AioBaseClient,
        semaphore: asyncio.Semaphore,
        column: int,
        row: int,
        frame_range: Optional[Tuple[int, int]] = None,
    ) -> np.ndarray:
        """
        Asynchronously retrieves data for a specific pixel over a range of frames.

        Parameters
        ----------
        s3 : aiobotocore.client.AioBaseClient
            The S3 client used to fetch the data.
        semaphore : asyncio.Semaphore
            The semaphore to limit the number of concurrent downloads.
        column : int
            The column index of the pixel (1-indexed).
        row : int
            The row index of the pixel (1-indexed).
        frame_range : tuple of int, optional
            The range of frames to consider, by default None.

        Returns
        -------
        data: np.ndarray
            The pixel values for the specified range of frames.
        """

        if frame_range is None:
            frame_range = (0, self.nframes)
        nframes = frame_range[1] - frame_range[0]
        if (nframes <= 0) | (nframes > self.nframes):
            raise ValueError(f"Can not parse `frame_range` of `{frame_range}`")

        async with semaphore:
            byte_range = self.find_byte_range_subframes(
                row=row, column=column, frame_range=frame_range
            )
            range_string = f"bytes={byte_range[0]}-{byte_range[1] - 1}"
            response = await s3.get_object(
                Bucket=BUCKET_NAME, Key=self.object_key, Range=range_string
            )
            data_bytes = await response["Body"].read()
            n_pixels = len(data_bytes) // BYTES_PER_PIX
            values = np.asarray(
                struct.unpack(">" + "f" * n_pixels, data_bytes)
            ).reshape((nframes, self.nsets))
            return values

    async def async_get_data_per_row(
        self,
        s3: AioBaseClient,
        semaphore: asyncio.Semaphore,
        start_column: int,
        ncolumns: int,
        row: int = 1,
    ):
        """
        Asynchronously retrieves data for a specific row over a range of columns.

        Parameters
        ----------
        s3 : aiobotocore.client.AioBaseClient
            The S3 client used to fetch the data.
        semaphore : asyncio.Semaphore
            The semaphore to limit the number of concurrent downloads.
        start_column : int
            The starting column index (1-indexed).
        ncolumns : int
            The number of columns to retrieve.
        row : int, optional
            The row index of the pixel (1-indexed), by default 1.

        Returns
        -------
        data: np.ndarray
            The pixel values for the specified row and columns.
        """
        async with semaphore:
            byte_range = self.find_byte_range(
                row=row, column_range=(start_column, start_column + ncolumns)
            )
            range_string = f"bytes={byte_range[0]}-{byte_range[1] - 1}"
            response = await s3.get_object(
                Bucket=BUCKET_NAME, Key=self.object_key, Range=range_string
            )
            data_bytes = await response["Body"].read()
            n_pixels = len(data_bytes) // BYTES_PER_PIX
            values = np.asarray(
                struct.unpack(">" + "f" * n_pixels, data_bytes)
            ).reshape((ncolumns, self.nframes, self.nsets))
            return values

    async def async_get_data_per_rows(
        self, runs: List[Tuple[int, int, int]]
    ) -> np.ndarray:
        """
        Asynchronously retrieves data for multiple rows.

        Parameters
        ----------
        runs : list of tuples
            A list of tuples containing parameters for each run.

        Returns
        -------
        data: np.ndarray
            The pixel values for the specified rows.
        """
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_DOWNLOADS)
        async with get_session().create_client(
            "s3", config=Config(signature_version=UNSIGNED)
        ) as s3:
            tasks = []
            for run in runs:
                task = asyncio.create_task(
                    self.async_get_data_per_row(s3, semaphore, **run)
                )
                tasks.append(task)
            values = await asyncio.gather(*tasks)
            return values

    async def async_get_data_per_block(
        self,
        start_column: int,
        ncolumns: int,
        start_row: int,
        nrows: int,
        frame_range: Optional[Tuple[int, int]] = None,
    ) -> np.ndarray:
        """
        Asynchronously retrieves data for a block of pixels.

        Parameters
        ----------
        start_column : int
            The starting column index (1-indexed).
        ncolumns : int
            The number of columns to retrieve.
        start_row : int
            The starting row index (1-indexed).
        nrows : int
            The number of rows to retrieve.
        frame_range : tuple of int, optional
            The range of frames to consider, by default None.

        Returns
        -------
        data: numpy.ndarray
            The pixel values for the specified block.
        """
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_DOWNLOADS)
        if frame_range is not None:
            async with get_session().create_client(
                "s3", config=Config(signature_version=UNSIGNED)
            ) as s3:
                tasks = []
                for row in np.arange(start_row, start_row + nrows):
                    for column in np.arange(start_column, start_column + ncolumns):
                        task = asyncio.create_task(
                            self.async_get_data_per_pixel(
                                s3,
                                semaphore,
                                column=column,
                                row=row,
                                frame_range=frame_range,
                            )
                        )
                        tasks.append(task)
                values = await asyncio.gather(*tasks)
                return values
        else:
            async with get_session().create_client(
                "s3", config=Config(signature_version=UNSIGNED)
            ) as s3:
                tasks = []
                for row in np.arange(start_row, start_row + nrows):
                    task = asyncio.create_task(
                        self.async_get_data_per_row(
                            s3,
                            semaphore,
                            start_column=start_column,
                            ncolumns=ncolumns,
                            row=row,
                        )
                    )
                    tasks.append(task)
                values = await asyncio.gather(*tasks)
                return values

    def get_pixel_timeseries(
        self,
        coordinates: Union[Tuple[int, int], List[Tuple[int, int]]],
        return_errors: bool = False,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Retrieves the pixel time series for the specified coordinates.

        Parameters
        ----------
        coordinates : tuple or list of tuple
            The coordinates (row, column) to retrieve time series for.
        return_errors : bool, optional
            Whether to return the errors alongside the flux, by default False.

        Returns
        -------
        data: numpy.ndarray or tuple of numpy.ndarray
            The flux (and optionally errors) time series for the specified coordinates.
        """
        if isinstance(coordinates, tuple):
            coordinates = [coordinates]
        runs = convert_coordinates_to_runs(coordinates)
        flux, flux_err = np.vstack(
            _sync_call(self.async_get_data_per_rows, runs=runs)
        ).transpose([2, 1, 0])
        if return_errors:
            return flux, flux_err
        return flux

    async def async_get_flux(
        self,
        start_column: int,
        ncolumns: int,
        start_row: int,
        nrows: int,
        frame_range: Optional[Tuple[int, int]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Asynchronously retrieves flux data for a block of pixels.

        Parameters
        ----------
        start_column : int
            The starting column index (1-indexed).
        ncolumns : int
            The number of columns to retrieve.
        start_row : int
            The starting row index (1-indexed).
        nrows : int
            The number of rows to retrieve.
        frame_range : tuple of int, optional
            The range of frames to consider, by default None.

        Returns
        -------
        tuple of numpy.ndarray
            The flux and flux error for the specified block.
        """
        if frame_range is not None:
            nframes = frame_range[1] - frame_range[0]
            if (nframes <= 0) | (nframes > self.nframes):
                raise ValueError(f"Can not parse `frame_range` of `{frame_range}`")
        else:
            nframes = self.nframes

        values = await self.async_get_data_per_block(
            start_column=start_column,
            ncolumns=ncolumns,
            start_row=start_row,
            nrows=nrows,
            frame_range=frame_range,
        )
        values = np.asarray(values).reshape((nrows, ncolumns, nframes, self.nsets))
        flux, flux_err = np.asarray(values).transpose([3, 2, 0, 1])
        return flux, flux_err

    @lru_cache(maxsize=128)
    def get_flux(
        self,
        corner: Tuple[int, int] = (1014, 1058),
        shape: Tuple[int, int] = (20, 21),
        frame_range: Optional[Tuple[int, int]] = None,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Retrieves the flux and flux error for a specific block of pixels.

        This method caches the result to avoid redundant calculations. The flux
        and flux error are retrieved for a rectangular block of pixels specified
        by the `corner` and `shape` parameters. The data can also be restricted
        to a specific range of frames using the `frame_range` parameter.

        Parameters
        ----------
        corner : tuple of int, optional
            The (row, column) coordinates of the top-left corner of the pixel block,
            by default (1014, 1058).
        shape : tuple of int, optional
            The (number of rows, number of columns) shape of the pixel block,
            by default (20, 21).
        frame_range : tuple of int, optional
            The (start, end) frame range to retrieve the data, by default None.

        Returns
        -------
        tuple of numpy.ndarray
            A tuple containing:
            - flux: A 3D numpy array of shape (nframes, nrows, ncolumns)
            containing the flux data.
            - flux_err: A 3D numpy array of shape (nframes, nrows, ncolumns)
            containing the flux error data.
        """
        flux, flux_err = _sync_call(
            self.async_get_flux,
            start_column=corner[1],
            ncolumns=shape[1],
            start_row=corner[0],
            nrows=shape[0],
            frame_range=frame_range,
        )
        return flux, flux_err


@lru_cache()
def get_primary_hdu(object_key):
    return _fix_primary_hdu(_sync_call(async_get_primary_hdu, object_key=object_key))


@lru_cache()
def get_last_hdu(object_key, end):
    return _sync_call(async_get_last_hdu, object_key=object_key, end=end)


async def async_get_primary_hdu(object_key: str) -> fits.PrimaryHDU:
    """
    Asynchronously retrieves the primary HDU of a cube FITS file.

    Parameters
    ----------
    object_key : str
        The S3 object key of the FITS file (i.e. the file name)

    Returns
    -------
    astropy.io.fits.PrimaryHDU
        The primary HDU of the FITS file, with corrected headers.
    """
    async with get_session().create_client(
        "s3", config=Config(signature_version=UNSIGNED)
    ) as s3:
        # Retrieve the cube header
        response = await s3.get_object(
            Bucket=BUCKET_NAME,
            Key=object_key,
            Range=f"bytes=0-{HDR_SIZE * 2-1}",
        )
        first_bytes = await response["Body"].read()
        with warnings.catch_warnings():
            # Ignore "File may have been truncated" warning
            warnings.simplefilter("ignore", AstropyUserWarning)
            with fits.open(BytesIO(first_bytes)) as hdulist:
                (
                    hdulist[0].header["NAXIS1"],
                    hdulist[0].header["NAXIS2"],
                    hdulist[0].header["NAXIS3"],
                    hdulist[0].header["NAXIS4"],
                ) = (
                    hdulist[1].header["NAXIS1"],
                    hdulist[1].header["NAXIS2"],
                    hdulist[1].header["NAXIS3"],
                    hdulist[1].header["NAXIS4"],
                )
                return hdulist[0]


async def async_get_last_hdu(object_key: str, end: int) -> fits.PrimaryHDU:
    """
    Asynchronously retrieves the last HDU of a cube FITS file.

    Parameters
    ----------
    object_key : str
        The S3 object key of the FITS file  (i.e. the file name)
    end : int
        The byte position to start retrieving data from.
    Returns
    -------
    astropy.io.fits.PrimaryHDU
        The last HDU of the FITS file.
    """
    # This seems wrong for a lot of sets
    async with get_session().create_client(
        "s3", config=Config(signature_version=UNSIGNED)
    ) as s3:
        # Retrieve the cube header
        response = await s3.get_object(
            Bucket=BUCKET_NAME, Key=object_key, Range=f"bytes={end}-"
        )
        first_bytes = await response["Body"].read()
        # n_pixels = len(first_bytes) // BYTES_PER_PIX
        # values = np.asarray(struct.unpack(">" + "f" * n_pixels, first_bytes))
        # return values
        return fits.open(
            BytesIO(first_bytes.lstrip(b"\x00")),
            ignore_missing_simple=True,
            lazy_load_hdus=False,
        )[0]


async def async_get_ffi(ffi_name: str) -> fits.HDUList:
    """
    Asynchronously retrieves a full-frame image (FFI) FITS file.

    Parameters
    ----------
    ffi_name : str
        The name of the full-frame image (FFI) as stored in s3.

    Returns
    -------
    astropy.io.fits.HDUList
        The HDU list of the retrieved FFI FITS file.
    """
    date, sector_str, camera, ccd, _, _ = ffi_name[4:].split("-")
    object_key = (
        f"tess/public/ffi/{sector_str}/{date[:4]}/{date[4:7]}/{camera}-{ccd}/{ffi_name}"
    )
    async with get_session().create_client(
        "s3", config=Config(signature_version=UNSIGNED)
    ) as s3:
        # Retrieve the cube header
        response = await s3.get_object(Bucket=BUCKET_NAME, Key=object_key)
        data_bytes = await response["Body"].read()
        return fits.open(
            BytesIO(data_bytes),
            ignore_missing_simple=True,
            lazy_load_hdus=False,
        )
