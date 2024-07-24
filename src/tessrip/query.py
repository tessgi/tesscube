import asyncio
import bz2
import json
import os
import struct
import warnings
from functools import lru_cache
from io import BytesIO
from typing import Optional, Union, Tuple, List

import numpy as np
from aiobotocore.session import get_session
from astropy.io import fits
from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from botocore import UNSIGNED
from botocore.config import Config

from . import (
    BUCKET_NAME,
    BYTES_PER_PIX,
    DATA_OFFSET,
    HDR_SIZE,
    MAX_CONCURRENT_DOWNLOADS,
    PACKAGEDIR,
    get_logger,
)
from .utils import (
    WCS_ATTRS,
    _extract_average_WCS,
    _fix_primary_hdu,
    _sync_call,
    convert_coordinates_to_runs,
    convert_to_native_types,
)

log = get_logger()


class Rip(object):
    """Ripper object to obtain portions of TESS data cube from MASTs AWS bucket

    This object will grab particular bytes out of the TESS data cubes to try to make obtaining data from the cubes efficient.
    """

    def __init__(self, sector: int, camera: int, ccd: int):
        self.sector, self.camera, self.ccd = sector, camera, ccd
        self.object_key = f"tess/public/mast/tess-s{sector:04}-{camera}-{ccd}-cube.fits"
        (
            self.nsets,
            self.nframes,
            self.ncolumns,
            self.nrows,
        ) = (
            self.primary_hdu.header["NAXIS1"],
            self.primary_hdu.header["NAXIS2"],
            self.primary_hdu.header["NAXIS3"],
            self.primary_hdu.header["NAXIS4"],
        )

    @property
    def center(self):
        if hasattr(self, "corner"):
            return (
                self.corner[0] + self.shape[0] / 2,
                self.corner[1] + self.shape[1] / 2,
            )
        else:
            raise ValueError("Run `get_flux` to specify a region.")

    @property
    def primary_hdu(self):
        return _primary_hdu(object_key=self.object_key)

    @property
    def last_hdu(self):
        end = (
            DATA_OFFSET
            + (self.ncolumns * self.nframes * self.nsets * self.nrows) * BYTES_PER_PIX
        )
        return _last_hdu(object_key=self.object_key, end=end)

    @property
    def ffi_names(self):
        return list(self.last_hdu.data["FFI_FILE"])

    def find_byte_offset(self, row: int, column: int, frame: int = 0) -> int:
        """Returns the byte offset of a specific pixel position."""
        # if (row < 1) | (row > 2049) | (column < 1) | (column > 2049):
        #     raise ValueError(f"`(row, column)` position `({row}, {column})` is outside of range.")
        # Subtract 1 from column and row because the byte location assumes zero-indexing,
        # whereas the TESS convention is to address column and row number with one-indexing.
        pixel_offset = ((column - 1) * self.nframes * self.nsets
                        + (row - 1) * self.ncolumns * self.nframes * self.nsets)
        pixel_offset += frame * self.nsets
        return DATA_OFFSET + BYTES_PER_PIX * pixel_offset

    def find_byte_range_subframes(
        self, row: int, column: int, frame_range: Optional[tuple] = None
    ) -> tuple:
        if frame_range is None:
            frame_range = (0, self.nframes)
        begin = self.find_byte_offset(row=row, column=column, frame=frame_range[0])
        end = self.find_byte_offset(row=row, column=column, frame=frame_range[1])
        return (begin, end)

    def find_byte_range(self, row: int, column_range: tuple) -> tuple:
        begin = self.find_byte_offset(row=row, column=column_range[0])
        end = self.find_byte_offset(row=row, column=column_range[1])
        return (begin, end)

    async def _async_get_data_per_pixel(
        self, s3, semaphore, column: int, row: int, frame_range: Optional[tuple] = None
    ):
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

    async def _async_get_data_per_row(
        self, s3, semaphore, start_column: int, ncolumns: int, row: int = 1
    ):
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

    async def _async_get_data_per_rows(self, runs):
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_DOWNLOADS)
        async with get_session().create_client(
            "s3", config=Config(signature_version=UNSIGNED)
        ) as s3:
            tasks = []
            for run in runs:
                task = asyncio.create_task(
                    self._async_get_data_per_row(s3, semaphore, **run)
                )
                tasks.append(task)
            values = await asyncio.gather(*tasks)
            return values

    async def _async_get_data_per_block(
        self,
        start_column: int,
        ncolumns: int,
        start_row: int,
        nrows: int,
        frame_range: Optional[tuple] = None,
    ):
        semaphore = asyncio.Semaphore(MAX_CONCURRENT_DOWNLOADS)
        if frame_range is not None:
            async with get_session().create_client(
                "s3", config=Config(signature_version=UNSIGNED)
            ) as s3:
                tasks = []
                for row in np.arange(start_row, start_row + nrows):
                    for column in np.arange(start_column, start_column + ncolumns):
                        task = asyncio.create_task(
                            self._async_get_data_per_pixel(
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
                        self._async_get_data_per_row(
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

    async def _async_get_data_per_coordinates(
        self,
        coordinates: list[tuple],
        frame_range: tuple
    ):

        semaphore = asyncio.Semaphore(MAX_CONCURRENT_DOWNLOADS)
        async with get_session().create_client(
            "s3", config=Config(signature_version=UNSIGNED)
        ) as s3:
            tasks = []
            for coordinate in coordinates:
                task = asyncio.create_task(
                    self._async_get_data_per_pixel(
                        s3,
                        semaphore,
                        column=coordinate[1],
                        row=coordinate[0],
                        frame_range=frame_range,
                        )
                )
                tasks.append(task)
            values = await asyncio.gather(*tasks)
            return values

    @staticmethod
    def to_column(flux, flux_err):
        tform = str(flux[0].size) + "E"
        dims = str(flux[0].shape[::-1])

        flux = fits.Column(
            name="FLUX",
            format=tform,
            dim=dims,
            unit="e-/s",
            disp="E14.7",
            array=flux
        )
        flux_err = fits.Column(
            name="FLUX_ERR",
            format=tform,
            dim=dims,
            unit="e-/s",
            disp="E14.7",
            array=flux_err,
        )
        return flux, flux_err

    def get_pixel_timeseries(
        self,
        coordinates: Union[List[Tuple[int, int]], SkyCoord],
        frame_range: Optional[Tuple[int, int]]=None
        ) -> Tuple[fits.Column, fits.Column]:
        if isinstance(coordinates, SkyCoord):
            coordinates = np.round(
                self.wcs.world_to_pixel(coordinates)
            ).astype(int)
        if frame_range is not None:
            npix = len(coordinates)
            nframes = frame_range[1] - frame_range[0]
            flux, flux_err = np.vstack(
                _sync_call(
                    self._async_get_data_per_coordinates,
                    coordinates=coordinates,
                    frame_range=frame_range
                    )
            ).reshape(npix, nframes, 2).transpose([2,1,0])

        else:
            runs = convert_coordinates_to_runs(coordinates)
            flux, flux_err = np.vstack(
                _sync_call(self._async_get_data_per_rows, runs=runs)
            ).transpose([2, 1, 0])

        return self.to_column(flux, flux_err)

    async def _async_get_flux(
        self,
        start_column: int,
        ncolumns: int,
        start_row: int,
        nrows: int,
        frame_range: Optional[tuple] = None,
    ):
        if frame_range is not None:
            nframes = frame_range[1] - frame_range[0]
            if (nframes <= 0) | (nframes > self.nframes):
                raise ValueError(f"Can not parse `frame_range` of `{frame_range}`")
        else:
            nframes = self.nframes

        values = await self._async_get_data_per_block(
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
        corner: tuple = (1014, 1058),
        shape: tuple = (20, 21),
        frame_range: Optional[tuple] = None,
    ):
        self.corner = corner
        self.shape = shape
        flux, flux_err = _sync_call(
            self._async_get_flux,
            start_column=corner[1],
            ncolumns=shape[1],
            start_row=corner[0],
            nrows=shape[0],
            frame_range=frame_range,
        )

        return self.to_column(flux, flux_err)

    def get_tpf(
        self,
        corner: tuple = (1014, 1014),
        shape: tuple = (20, 21),
        frame_range: Optional[tuple] = None,
    ):
        flux, flux_err = self.get_flux(
            corner=corner, shape=shape, frame_range=frame_range
        )

        cols = [
            self.time,
            self.timecorr,
            self.cadence_number,
            self.quality,
            flux,
            flux_err,
        ]

        table_hdu = fits.BinTableHDU.from_columns(cols)
        table_hdu.header["EXTNAME"] = "PIXELS"

        aperture_hdu = fits.ImageHDU(data=np.ones(shape))
        aperture_hdu.header["EXTNAME"] = "APERTURE"
        for kwd, val, cmt in self.wcs.to_header(relax=True).cards:
            aperture_hdu.header.set(kwd, val, cmt)

        # Adding extra aperture keywords (TESS specific)
        aperture_hdu.header.set(
            "NPIXMISS", None, "Number of op. aperture pixels not collected"
        )
        aperture_hdu.header.set("NPIXSAP", None, "Number of pixels in optimal aperture")
        # Need to fix NAXIS in primary hdu
        hdulist = fits.HDUList([self.primary_hdu, table_hdu, aperture_hdu])
        return hdulist

    @property
    def time(self):
        return fits.Column(
            name="TIME",
            format="D",
            unit="BJD - 2457000, days",
            disp="D14.7",
            array=(self.last_hdu.data["TSTART"] + self.last_hdu.data["TSTOP"]) / 2,
        )

    @property
    def timecorr(self):
        return fits.Column(
            name="TIMECORR",
            format="E",
            unit="d",
            disp="E14.7",
            array=self.last_hdu.data["BARYCORR"],
        )

    @property
    def quality(self):
        return fits.Column(
            name="QUALITY",
            format="J",
            disp="B16.16",
            array=self.last_hdu.data["DQUALITY"],
        )

    @property
    def n_frames(self):
        return fits.Column(
            name="NUM_FRM",
            format="I",
            array=self.last_hdu.data["NUM_FRM"],
        )

    @property
    def cadence_number(self):
        cadence_number = np.cumsum(
            np.round(
                np.diff(self.last_hdu.data["TSTART"])
                / np.median(self.last_hdu.data["TELAPSE"])
            ).astype(int)
        )
        return fits.Column(name="CADENCENO", format="I", array=cadence_number)

    @property
    def wcs(self):
        return _extract_average_WCS(self.last_hdu)

    def _save_wcss(self, dir=None):
        if dir is None:
            dir = f"{PACKAGEDIR}/data/s{self.sector:04}/"
        os.makedirs(dir, exist_ok=True)
        hdu = self.last_hdu
        wcs_dict = {
            attr: hdu.data[attr].tolist()
            if isinstance(hdu.data[attr][0], (float, np.integer, int))
            else hdu.data[attr][0]
            for attr in WCS_ATTRS(hdu)
        }
        wcs_dict = convert_to_native_types(wcs_dict)
        filename = f"{dir}TESS_wcs_sector{self.sector:04}_cam{self.camera}_ccd{self.ccd}.json.bz2"
        wcs_dict["CTYPE1"] = "RA---TAN-SIP"
        wcs_dict["CTYPE2"] = "DEC--TAN-SIP"
        json_data = json.dumps(wcs_dict)
        with bz2.open(filename, "wt", encoding="utf-8") as f:
            f.write(json_data)

    def _load_wcss(self, dir=None):
        if dir is None:
            dir = f"{PACKAGEDIR}/data/s{self.sector:04}/"
        filename = f"{dir}TESS_wcs_sector{self.sector:04}_cam{self.camera}_ccd{self.ccd}.json.bz2"
        if not os.path.isfile(filename):
            self._save_wcss()
        with bz2.open(filename, "rt", encoding="utf-8") as f:
            loaded_dict = json.load(f)
        wcs_attrs = WCS_ATTRS(self.last_hdu)
        hdr = fits.PrimaryHDU().header

        def _get_wcs(idx):
            wcs_dict = {
                attr: loaded_dict[attr]
                if (not isinstance(loaded_dict[attr], list))
                else loaded_dict[attr][idx]
                for attr in wcs_attrs
            }
            for attr in wcs_attrs:
                if not isinstance(loaded_dict[attr], list):
                    hdr[attr] = loaded_dict[attr]
                else:
                    if not np.isfinite(loaded_dict[attr][idx]):
                        return None
                    hdr[attr] = loaded_dict[attr][idx]
            return WCS(wcs_dict, relax=True)

        return {idx: _get_wcs(idx) for idx in range(self.nframes)}

    @lru_cache(maxsize=128)
    def _wcss(self):
        return self._load_wcss()

    @property
    def wcss(self):
        return self._wcss()


@lru_cache()
def _primary_hdu(object_key):
    return _fix_primary_hdu(_sync_call(async_get_primary_hdu, object_key=object_key))


@lru_cache()
def _last_hdu(object_key, end):
    return _sync_call(async_get_last_hdu, object_key=object_key, end=end)


async def async_get_primary_hdu(object_key):
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


async def async_get_last_hdu(object_key, end):
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
