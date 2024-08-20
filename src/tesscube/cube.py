from functools import lru_cache
from typing import Optional, Union

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.wcs.utils import fit_wcs_from_points

from . import BYTES_PER_PIX, DATA_OFFSET, log
from .fits import (get_header_dict, get_output_first_extention_header,
                   get_output_primary_hdu, get_output_second_extension_header,
                   get_wcs_header_by_extension)
from .query import QueryMixin, async_get_ffi, get_last_hdu, get_primary_hdu
from .utils import _sync_call
from .wcs import WCSMixin


class TESSCube(QueryMixin, WCSMixin):
    """Cube object to obtain portions of TESS data cube from MASTs AWS bucket

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

    def __repr__(self):
        return f"TESSCube [Sector {self.sector}, Camera {self.camera}, CCD {self.ccd}]"

    @property
    def header_dict(self):
        """Get the header keywords from the MAST cube."""
        return get_header_dict(self)

    @staticmethod
    def from_skycoord(coord: SkyCoord, sector: int):
        # Grabs a cube given a SkyCoord by querying the WCS.
        for camera in np.arange(1, 5):
            for ccd in np.arange(1, 5):
                cube = TESSCube(sector=sector, camera=camera, ccd=ccd)
                if cube.wcs.footprint_contains(coord):
                    return cube

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
        return get_primary_hdu(object_key=self.object_key)

    @property
    def last_hdu(self):
        end = (
            DATA_OFFSET
            + (self.ncolumns * self.nframes * self.nsets * self.nrows) * BYTES_PER_PIX
        )
        return get_last_hdu(object_key=self.object_key, end=end)

    @property
    def ffi_names(self):
        return list(self.last_hdu.data["FFI_FILE"])

    @lru_cache(maxsize=4)
    def get_ffi(
        self, index: int = None, time: Time = None, ffi_name: str = None, raw=False
    ):
        provided_args = [arg is not None for arg in (time, index, ffi_name)]
        if sum(provided_args) != 1:
            raise ValueError(
                "You must provide exactly one of `time`, `index`, or `ffi_name`."
            )

        if time is not None:
            start = Time(self.last_hdu.data["TSTART"] + 2457000, format="jd")
            end = Time(self.last_hdu.data["TSTART"] + 2457000, format="jd")
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
    ):
        if isinstance(target, tuple):
            corner = target
            target = SkyCoord(*self.wcs.all_pix2world([corner], 0)[0], unit="deg")
        elif isinstance(target, SkyCoord):
            if not self.wcs.footprint_contains(target):
                raise ValueError(
                    f"Target {target} not in Sector {self.sector}, Camera {self.camera}, CCD {self.ccd}."
                )
            corner = (
                np.asarray(self.wcs.world_to_pixel(target))[::-1] - np.asarray(shape) // 2
            )
            corner = np.floor(corner).astype(int)
            corner = (corner[0] + 1, corner[1] + 1)
        else:
            raise ValueError("Pass an origin coordinate or a SkyCoord object")

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

        pos_corr1, pos_corr2 = self.get_poscorr(target)

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
            time = (self.last_hdu.data['tstop'][k][mask][::frame_bin] + self.last_hdu.data['tstart'][k][mask][(frame_bin-1)::frame_bin])/2
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
                    *get_output_first_extention_header(self).cards,
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
                [*get_output_second_extension_header(self).cards, *wcs_header.cards]
            ),
        )
        aperture_hdu.header["EXTNAME"] = "APERTURE"
        aperture_hdu.header.set(
            "NPIXMISS", None, "Number of op. aperture pixels not collected"
        )
        aperture_hdu.header.set("NPIXSAP", None, "Number of pixels in optimal aperture")

        # Need to fix NAXIS in primary hdu
        hdulist = fits.HDUList([get_output_primary_hdu(self), table_hdu, aperture_hdu])
        return hdulist

    @property
    def time(self):
        return (self.last_hdu.data["TSTART"] + self.last_hdu.data["TSTOP"]) / 2

    @property
    def timecorr(self):
        return self.last_hdu.data["BARYCORR"]

    @property
    def quality(self):
        return self.last_hdu.data["DQUALITY"]

    @property
    def cadence_number(self):
        cadence_number = np.cumsum(
            np.round(
                np.diff(self.last_hdu.data["TSTART"])
                / np.median(self.last_hdu.data["TELAPSE"])
            ).astype(int)
        )
        return np.hstack([cadence_number, cadence_number[-1] + 1])

    @property
    def exposure_time(self):
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