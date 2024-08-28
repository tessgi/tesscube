from functools import lru_cache, cached_property
from typing import Optional, Union, Dict

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs.utils import fit_wcs_from_points

from . import BYTES_PER_PIX, DATA_OFFSET, log
from .fits import (
    get_header_dict,
    get_output_first_extention_header,
    get_output_primary_hdu,
    get_output_second_extension_header,
    get_wcs_header_by_extension,
)
from .query import QueryMixin, async_get_ffi, get_last_hdu, get_primary_hdu
from .utils import _sync_call, validate_tuple
from .wcs import WCSMixin, WCS_ATTRS


class TESSCube(QueryMixin, WCSMixin):
    """
    A class to obtain portions of TESS data cubes from MAST's AWS bucket.

    This class provides methods to efficiently retrieve specific portions
    of TESS (Transiting Exoplanet Survey Satellite) FFI data cubes, enabling
    users to extract and analyze the desired segments of data.

    Parameters
    ----------
    sector : int
        The TESS observation sector number.
    camera : int
        The camera number (1-4).
    ccd : int
        The CCD number (1-4).
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

        self.output_primary_ext = get_output_primary_hdu(self)
        self.output_first_header = get_output_first_extention_header(self)
        self.output_secondary_header = get_output_second_extension_header(self)
        self.wcs_attrs_no_sip = WCS_ATTRS(self.last_hdu, sip=False)
        self.wcs_attrs_no_sip = list(
            set(self.wcs_attrs_no_sip) - set(["WCSAXES", "WCSAXESP"])
        )

    def __repr__(self):
        return f"TESSCube [Sector {self.sector}, Camera {self.camera}, CCD {self.ccd}]"

    @cached_property
    def header_dict(self) -> Dict:
        """
        Get the header keywords from the MAST cube.

        Returns
        -------
        dict
            A dictionary of header keywords from the MAST cube.
        """
        return get_header_dict(self)

    @staticmethod
    def from_skycoord(coord: SkyCoord, sector: int):
        """
        Create a TESSCube object from a SkyCoord object.

        This method identifies the appropriate TESS cube that contains the
        given SkyCoord and returns a TESSCube object.

        Parameters
        ----------
        coord : SkyCoord
            The celestial coordinates of the target.
        sector : int
            The TESS observation sector number.

        Returns
        -------
        TESSCube
            The TESSCube object containing the specified coordinates.
        """
        # Grabs a cube given a SkyCoord by querying the WCS.
        for camera in np.arange(1, 5):
            for ccd in np.arange(1, 5):
                cube = TESSCube(sector=sector, camera=camera, ccd=ccd)
                if cube.wcs.footprint_contains(coord):
                    return cube

    @cached_property
    def primary_hdu(self):
        """The primary HDU of the cube file."""
        return get_primary_hdu(object_key=self.object_key)

    @cached_property
    def last_hdu(self):
        """The last HDU of the cube file."""
        end = (
            DATA_OFFSET
            + (self.ncolumns * self.nframes * self.nsets * self.nrows) * BYTES_PER_PIX
        )
        return get_last_hdu(object_key=self.object_key, end=end)

    @cached_property
    def ffi_names(self):
        """The FFI names used to make the cube."""
        return list(self.last_hdu.data["FFI_FILE"])

    @cached_property
    def tstart(self):
        return self.last_hdu.data["TSTART"]

    @cached_property
    def tstop(self):
        return self.last_hdu.data["TSTOP"]

    @cached_property
    def telapse(self):
        return self.last_hdu.data["TELAPSE"]

    @lru_cache(maxsize=4)
    def get_ffi(
        self,
        index: Optional[int] = None,
        time: Optional[Time] = None,
        ffi_name: Optional[str] = None,
        raw: bool = False,
    ) -> fits.HDUList:
        """
        Retrieve a Full-Frame Image (FFI) from the data cube.

        One of `time`, `index`, or `ffi_name` must be provided to identify the FFI.
        This method caches up to 4 FFIs for efficient reuse.

        Parameters
        ----------
        index : int, optional
            The index of the FFI to retrieve.
        time : astropy.time.Time, optional
            The time at which the FFI was captured.
        ffi_name : str, optional
            The name of the FFI file.
        raw : bool, optional
            Whether to retrieve the raw FFI (default is False).

        Returns
        -------
        astropy.io.fits.HDUList
            The HDUList containing the FFI data.
        """
        provided_args = [arg is not None for arg in (time, index, ffi_name)]
        if sum(provided_args) != 1:
            raise ValueError(
                "You must provide exactly one of `time`, `index`, or `ffi_name`."
            )

        if time is not None:
            start = Time(self.tstart + 2457000, format="jd")
            end = Time(self.tstop + 2457000, format="jd")
            t = Time.now()
            t = Time(2459343.87182313, format="jd")
            if not ((t > start).any() & (t < end).any()):
                raise ValueError(
                    f"Input time is not during Sector {self.sector} observations."
                )
            ffi_name = self.ffi_names[np.where((t > start))[0][-1]]

        if index is not None:
            ffi_name = self.ffi_names[index]
        if raw:
            return _sync_call(async_get_ffi, ffi_name.replace("ffic", "ffir"))
        return _sync_call(async_get_ffi, ffi_name)

    def get_tpf(
        self,
        target: Union[tuple, SkyCoord] = (1014, 1014),
        shape: tuple = (20, 21),
        frame_range: Optional[tuple] = None,
        frame_bin: int = 1,
        calculate_poscorr=True,
    ) -> fits.HDUList:
        """
        Retrieve a Target Pixel File (TPF) from the data cube.

        The TPF contains the flux data for a specific target region, defined by
        either pixel coordinates or a SkyCoord object. The data can be binned
        over multiple frames (cadences) using `frame_bin`.

        Parameters
        ----------
        target : tuple or SkyCoord, optional
            The (row, column) coordinates or a SkyCoord object for the target region.
        shape : tuple of int, optional
            The (number of rows, number of columns) shape of the target region.
        frame_range : tuple of int, optional
            The (start, end) frame range to retrieve the data.
        frame_bin : int, optional
            The number of frames to bin together (default is 1, i.e. no binning).
        calculate_poscorr : bool
            Whether to calculate the `POS_CORR` columns. These are present in SPOC
            data but not present in some HLSPs. This takes extra time to calculate,
            but provides extra position information. Setting to False will result in
            POS_CORR1 and POS_CORR2 being set to zero.

        Returns
        -------
        astropy.io.fits.HDUList
            The HDUList containing the TPF data.
        """
        if isinstance(target, tuple):
            corner = validate_tuple(target)
            if calculate_poscorr:
                target = SkyCoord(*self.wcs.all_pix2world([corner], 0)[0], unit="deg")
        elif isinstance(target, SkyCoord):
            if not self.wcs.footprint_contains(target):
                raise ValueError(
                    f"Target {target} not in Sector {self.sector}, Camera {self.camera}, CCD {self.ccd}."
                )
            corner = (
                np.asarray(self.wcs.world_to_pixel(target))[::-1]
                - np.asarray(shape) // 2
            )
            corner = np.floor(corner).astype(int)
            corner = (corner[0] + 1, corner[1] + 1)
        else:
            raise ValueError("Pass an origin coordinate or a SkyCoord object")
        shape = validate_tuple(shape)
        R, C = np.meshgrid(
            np.arange(corner[0], corner[0] + shape[0]),
            np.arange(corner[1], corner[1] + shape[1]),
        )
        c = self.wcs.pixel_to_world(C.ravel(), R.ravel())
        wcs = fit_wcs_from_points(
            (C.ravel() - corner[1] + 1, R.ravel() - corner[0] + 1), c
        )
        wcs_header = wcs.to_header(relax=True)

        # Adding the physical wcs keywords
        wcs_header.set("CRVAL1P", corner[1], "value at reference CCD column")
        wcs_header.set("CRVAL2P", corner[0], "value at reference CCD row")

        wcs_header.set(
            "WCSNAMEP", "PHYSICAL", "name of world coordinate system alternate P"
        )
        wcs_header.set("WCSAXESP", 2, "number of WCS physical axes")

        wcs_header.set("CTYPE1P", "RAWX", "physical WCS axis 1 type CCD col")
        wcs_header.set("CUNIT1P", "PIXEL", "physical WCS axis 1 unit")
        wcs_header.set("CRPIX1P", 1, "reference CCD column")
        wcs_header.set("CDELT1P", 1.0, "physical WCS axis 1 step")

        wcs_header.set("CTYPE2P", "RAWY", "physical WCS axis 2 type CCD col")
        wcs_header.set("CUNIT2P", "PIXEL", "physical WCS axis 2 unit")
        wcs_header.set("CRPIX2P", 1, "reference CCD row")
        wcs_header.set("CDELT2P", 1.0, "physical WCS axis 2 step")

        flux, flux_err = self.get_flux(
            corner=corner, shape=shape, frame_range=frame_range
        )

        tform = str(flux[0].size) + "E"
        dims = str(flux[0].shape[::-1])

        if calculate_poscorr:
            pos_corr1, pos_corr2 = self.get_poscorr(target)
        else:
            pos_corr1, pos_corr2 = np.zeros((2, flux.shape[0]))

        if frame_range is not None:
            k = np.zeros(len(self), bool)
            k[frame_range[0] : frame_range[1]] = True
        else:
            k = np.ones(len(self), bool)

        time, timecorr, cadenceno, quality, pos_corr1, pos_corr2 = (
            self.time[k],
            self.timecorr[k],
            self.cadence_number[k],
            self.quality[k],
            pos_corr1[k],
            pos_corr2[k],
        )

        if frame_bin > 1:
            idxs = np.hstack(
                [
                    r[: frame_bin * np.floor(len(r) / frame_bin).astype(int)]
                    for r in np.array_split(
                        cadenceno, np.where(np.diff(cadenceno) != 1)[0]
                    )
                ]
            )
            mask = np.in1d(cadenceno, idxs)
            time = (
                self.tstop[k][mask][::frame_bin]
                + self.tstart[k][mask][(frame_bin - 1) :: frame_bin]
            ) / 2
            timecorr = np.nanmean(
                np.asarray(
                    [timecorr[mask][idx::frame_bin] for idx in range(frame_bin)]
                ),
                axis=0,
            )
            cadenceno = cadenceno[mask][::frame_bin]
            quality = np.bitwise_or.reduce(
                np.asarray([quality[mask][idx::frame_bin] for idx in range(frame_bin)]),
                axis=0,
            )
            pos_corr1 = np.nanmean(
                np.asarray(
                    [pos_corr1[mask][idx::frame_bin] for idx in range(frame_bin)]
                ),
                axis=0,
            )
            pos_corr2 = np.nanmean(
                np.asarray(
                    [pos_corr2[mask][idx::frame_bin] for idx in range(frame_bin)]
                ),
                axis=0,
            )
            flux = np.nansum(
                np.asarray([flux[mask][idx::frame_bin] for idx in range(frame_bin)]),
                axis=0,
            )
            flux_err = (
                np.nansum(
                    np.asarray(
                        [
                            flux_err[mask][idx::frame_bin] ** 2
                            for idx in range(frame_bin)
                        ]
                    ),
                    axis=0,
                )
                ** 0.5
            )

        cols = [
            fits.Column(
                name="TIME",
                format="D",
                unit="BJD - 2457000, days",
                disp="D14.7",
                array=time,
            ),
            fits.Column(
                name="TIMECORR",
                format="E",
                unit="d",
                disp="E14.7",
                array=timecorr,
            ),
            fits.Column(name="CADENCENO", format="I", array=cadenceno),
            # This is included to give the files the same structure as the SPOC files
            fits.Column(
                name="RAW_CNTS",
                format=tform,
                dim=dims,
                unit="e-/s",
                disp="E14.7",
                array=np.zeros_like(flux),
            ),
            fits.Column(
                name="FLUX",
                format=tform,
                dim=dims,
                unit="e-/s",
                disp="E14.7",
                array=flux,
            ),
            fits.Column(
                name="FLUX_ERR",
                format=tform,
                dim=dims,
                unit="e-/s",
                disp="E14.7",
                array=flux_err,
            ),
            fits.Column(
                name="FLUX_BKG",
                format=tform,
                dim=dims,
                unit="e-/s",
                disp="E14.7",
                array=np.zeros_like(flux),
            ),
            fits.Column(
                name="FLUX_BKG_ERR",
                format=tform,
                dim=dims,
                unit="e-/s",
                disp="E14.7",
                array=np.zeros_like(flux_err),
            ),
            fits.Column(
                name="QUALITY",
                format="J",
                disp="B16.16",
                array=quality,
            ),
            # This can be updated using the WCS arrays
            fits.Column(
                name="POS_CORR1",
                format="E",
                unit="pixel",
                disp="E14.7",
                array=pos_corr1,
            ),
            fits.Column(
                name="POS_CORR2",
                format="E",
                unit="pixel",
                disp="E14.7",
                array=pos_corr2,
            ),
        ]

        table_hdu = fits.BinTableHDU.from_columns(
            cols,
            header=fits.Header(
                [
                    *self.output_first_header.cards,
                    *get_wcs_header_by_extension(wcs_header, ext=4).cards,
                    *get_wcs_header_by_extension(wcs_header, ext=5).cards,
                    *get_wcs_header_by_extension(wcs_header, ext=6).cards,
                    *get_wcs_header_by_extension(wcs_header, ext=7).cards,
                    *get_wcs_header_by_extension(wcs_header, ext=8).cards,
                ]
            ),
        )
        table_hdu.header["EXTNAME"] = "PIXELS"
        if frame_bin > 1:
            table_hdu.header.set("EXPOSURE", table_hdu.header["EXPOSURE"] * frame_bin)
            table_hdu.header.set("NUM_FRM", table_hdu.header["NUM_FRM"] * frame_bin)
            table_hdu.header.set("TIMEDEL", table_hdu.header["TIMEDEL"] * frame_bin)

        aperture_hdu = fits.ImageHDU(
            data=np.ones(shape),
            header=fits.Header(
                [*self.output_secondary_header.cards, *wcs_header.cards]
            ),
        )
        aperture_hdu.header["EXTNAME"] = "APERTURE"
        aperture_hdu.header.set(
            "NPIXMISS", None, "Number of op. aperture pixels not collected"
        )
        aperture_hdu.header.set("NPIXSAP", None, "Number of pixels in optimal aperture")

        # Need to fix NAXIS in primary hdu
        hdulist = fits.HDUList([self.output_primary_ext, table_hdu, aperture_hdu])
        return hdulist

    @cached_property
    def time(self):
        """Time of the frames in BTJD. Note this is the time at the center of the FFI."""
        return (self.tstart + self.tstop) / 2

    @cached_property
    def timecorr(self):
        """Barycentric time correction for the center of the FFI."""
        return self.last_hdu.data["BARYCORR"]

    @cached_property
    def quality(self):
        """SPOC provided quality flags for each cadence."""
        return self.last_hdu.data["DQUALITY"]

    @cached_property
    def cadence_number(self):
        """Cadence number for each frame. Note this is not the same as the SPOC provided cadence numbers."""
        cadence_number = np.cumsum(
            np.round(np.diff(self.tstop) / np.median(self.telapse)).astype(int)
        )
        return np.hstack([cadence_number, cadence_number[-1] + 1])

    @cached_property
    def exposure_time(self):
        """Exposure time in days"""
        return self.last_hdu.data["EXPOSURE"][0]

    @property
    def shape(self):
        return (self.nframes, self.nrows, self.ncolumns)

    def __getitem__(self, other):
        if isinstance(other, int):
            return self.get_ffi(index=other)
        elif isinstance(other, (slice, np.ndarray, list, range)):
            raise ValueError("Can not return time-slices of entire cube.")
        elif isinstance(other, tuple):
            if len(other) != 3:
                other = (*other, slice(0, self.ncolumns))
            if not isinstance(other[0], slice):
                if isinstance(other[0], int):
                    other = (slice(other[0], other[0] + 1), *other[1:])
                else:
                    raise ValueError("Pass a time slice.")
            if (other[0].start is None) & (other[0].stop is None):
                frame_range = None
            else:
                time_index = (
                    other[0].start if other[0].start is not None else 0,
                    other[0].stop if other[0].stop is not None else self.nframes,
                )
                if time_index[1] - time_index[0] > 20:
                    log.warn(
                        "Time slices are inefficient, consider downloading all data and slicing locally."
                    )
                frame_range = (time_index[0], time_index[1])
            frame_bin = other[0].step if other[0].step is not None else 1
            if isinstance(other[1], slice):
                row_index = (
                    other[1].start if other[1].start is not None else 0,
                    other[1].stop if other[1].stop is not None else self.nrows,
                )
            else:
                raise ValueError("Pass columns as a slice.")
            if isinstance(other[2], slice):
                col_index = (
                    other[2].start if other[2].start is not None else 0,
                    other[2].stop if other[2].stop is not None else self.ncolumns,
                )
            else:
                raise ValueError("Pass rows as a slice.")
            for sl in other[1:]:
                if (sl.step is not None) & (sl.step != 1):
                    raise ValueError("Can not skip rows and columns in slicing.")
            return self.get_tpf(
                target=(row_index[0] + 1, col_index[0] + 1),
                shape=(row_index[1] - row_index[0], col_index[1] - col_index[0]),
                frame_range=frame_range,
                frame_bin=frame_bin,
            )
        else:
            raise ValueError(f"Can not index into `TESSCube` using {other}.")

    def __len__(self):
        return self.nframes
