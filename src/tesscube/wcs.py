"""Tools to work with WCS"""

import bz2
import json
import os
from functools import lru_cache

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from . import PACKAGEDIR
from .utils import convert_to_native_types

WCS_ATTRS_STARTS = [
    "CTYPE",
    "CRVAL",
    "CRPIX",
    "CUNIT",
    "NAXIS",
    "CD1",
    "CD2",
    "CDELT",
    "WCS",
    "1P",
    "2P",
    "A_",
    "AP_",
    "B_",
    "BP_",
]


def WCS_ATTRS(hdu, sip=True):
    return np.hstack(
        [
            *[
                [key for key in hdu.header.keys() if key.startswith(keystart)]
                for keystart in [WCS_ATTRS_STARTS if sip else WCS_ATTRS_STARTS[:-4]][0]
            ],
        ]
    ).tolist()


def _extract_average_WCS(hdu):
    wcs_hdu = fits.PrimaryHDU()
    for attr in WCS_ATTRS(hdu):
        if not hdu.columns[attr].format.endswith("A"):
            wcs_hdu.header[attr] = np.nanmedian(hdu.data[attr])
        else:
            data = hdu.data[attr]
            value = ""
            idx = 0
            while value == "":
                value = data[idx]
                idx += 1
            wcs_hdu.header[attr] = value
    wcs_hdu.header["WCSAXES"] = int(wcs_hdu.header["WCSAXES"])
    wcs_hdu.header["WCSAXESP"] = int(wcs_hdu.header["WCSAXESP"])
    return WCS(wcs_hdu.header)


# def _extract_all_WCS(hdu):
#     """Extract all the WCSs from a TableHDU from the TESS cube"""
#     wcss = []
#     for idx in np.arange(len(hdu.data["CRPIX1"])):
#         try:
#             wcs_hdu = fits.PrimaryHDU()
#             for attr in WCS_ATTRS(hdu):
#                 wcs_hdu.header[attr] = hdu[0].data[attr][idx]
#             wcs_hdu.header["WCSAXES"] = int(wcs_hdu.header["WCSAXES"])
#             wcs_hdu.header["WCSAXESP"] = int(wcs_hdu.header["WCSAXESP"])
#             wcss.append(WCS(wcs_hdu.header))
#         except:
#             wcss.append(None)
#     return wcss


class WCSMixin:
    """Mixins to use the WCS"""

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
            if isinstance(np.nanmedian(hdu.data[attr]), (float, np.integer, int))
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

    def get_poscorr(self, coord):
        hdu = self.last_hdu
        wcs0 = WCS(
            {
                attr: hdu.data[attr][0]
                if isinstance(hdu.data[attr][0], str)
                else np.nan_to_num(hdu.data[attr][0])
                for attr in WCS_ATTRS(hdu, sip=False)
            }
        )
        pos_corr1, pos_corr2 = np.zeros((2, len(self))) * np.nan
        for idx in range(self.shape[0]):
            crval = np.asarray([hdu.data["CRVAL1"][idx], hdu.data["CRVAL2"][idx]])
            cd = np.asarray(
                [
                    [hdu.data["CD1_1"][idx], hdu.data["CD2_1"][idx]],
                    [hdu.data["CD1_2"][idx], hdu.data["CD2_2"][idx]],
                ]
            ).T
            if np.isfinite(cd).all() & np.isfinite(crval).all():
                wcs0.wcs.crval = crval
                wcs0.wcs.cd = cd
                pos_corr1[idx], pos_corr2[idx] = wcs0.wcs_world2pix(
                    [(coord.ra.deg, coord.dec.deg)], 0
                )[0]
        pos_corr1_0, pos_corr2_0 = self.wcs.wcs_world2pix(
            [(coord.ra.deg, coord.dec.deg)], 0
        )[0]
        return pos_corr1 - pos_corr1_0, pos_corr2 - pos_corr2_0
