"""Utilities to help rip data"""
import asyncio
import functools
import sys
from io import BytesIO

import numpy as np
from aiobotocore.session import get_session
from astropy.io import fits
from astropy.wcs import WCS
from botocore import UNSIGNED
from botocore.config import Config

from . import BUCKET_NAME, __version__
from .config import get_logger

log = get_logger()


def convert_to_native_types(obj):
    """
    Recursively convert objects in a data structure to native Python types.
    Handles dictionaries, lists, and NumPy data types.
    """
    if isinstance(obj, np.integer):
        return int(obj)  # Convert NumPy scalar to a native Python type
    elif isinstance(obj, np.ndarray):
        return obj.tolist()  # Convert NumPy arrays to Python list
    elif isinstance(obj, dict):
        return {key: convert_to_native_types(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_to_native_types(item) for item in obj]
    else:
        return obj  # Return the object as is for native Python types


WCS_ATTRS_STARTS = [
    "CTYPE",
    "CRVAL",
    "CRPIX",
    "CUNIT",
    "NAXIS",
    "CD1",
    "CD2",
    "CDELT",
    "A_",
    "AP_",
    "B_",
    "BP_",
    "WCS",
    "1P",
    "2P",
]


def WCS_ATTRS(hdu):
    return np.hstack(
        [
            *[
                [key for key in hdu.header.keys() if key.startswith(keystart)]
                for keystart in WCS_ATTRS_STARTS
            ],
        ]
    ).tolist()


# Flag to indicate if nest_asyncio has been applied
_nest_asyncio_applied = False


def _sync_call(func, *args, **kwargs):
    global _nest_asyncio_applied
    # Check if we're in a Jupyter notebook environment
    if "ipykernel" in sys.modules and not _nest_asyncio_applied:
        # We are in Jupyter, check for nest_asyncio
        try:
            import nest_asyncio

            nest_asyncio.apply()
            _nest_asyncio_applied = True  # Set the flag so we don't apply it again
        except ImportError:
            log.warn(
                "nest_asyncio is required in a Jupyter environment. Please install with `!pip install nest_asyncio`."
            )
            return None
        # Run the async function with the current event loop
        return asyncio.get_event_loop().run_until_complete(func(*args, **kwargs))
    else:
        # We are not in Jupyter or nest_asyncio has already been applied, use asyncio.run()
        return asyncio.run(func(*args, **kwargs))


def convert_coordinates_to_runs(coordinates):
    """Converts a list of (row, column) coordinates to a list of (start_row, end_row, column) coordinates."""
    # Step 1: Sort the coordinates by row and then by column.
    sorted_coords = sorted(coordinates, key=lambda x: (x[0], x[1]))

    result = []
    start_column, count, row = None, 0, None

    for r, col in sorted_coords:
        if row is None:
            # Initialize the first run
            row, start_column, count = r, col, 1
        elif r == row and col == start_column + count:
            # Continuation of a run
            count += 1
        else:
            # End of a run, add to result and start a new run
            result.append({"start_column": start_column, "ncolumns": count, "row": row})
            row, start_column, count = r, col, 1

    # Add the last run
    result.append({"start_column": start_column, "ncolumns": count, "row": row})

    return result


def _fix_primary_hdu(hdu):
    """Fixes the primary HDU header to mimic TESS official products.

    This is largely taken from MASTs astrocut package under the BSD license.
    """
    log.debug("Fixing primary HDU.")
    hdu0 = hdu.copy()
    hdr = hdu0.header
    hdr["CREATOR"] = ("tessrip", "software used to produce this file")
    hdr["PROCVER"] = (__version__, "software version")
    hdr["FFI_TYPE"] = ("SPOC", "the FFI type used to make the cutouts")
    # We can't update this without RA and Dec
    # hdr["RA_OBJ"] = (ra, "[deg] right ascension")
    # hdr["DEC_OBJ"] = (dec, "[deg] declination")
    hdr["TIMEREF"] = (
        "SOLARSYSTEM",
        "barycentric correction applied to times",
    )
    hdr["TASSIGN"] = ("SPACECRAFT", "where time is assigned")
    hdr["TIMESYS"] = (
        "TDB",
        "time system is Barycentric Dynamical Time (TDB)",
    )
    hdr["BJDREFI"] = (2457000, "integer part of BTJD reference date")
    hdr["BJDREFF"] = (
        0.00000000,
        "fraction of the day in BTJD reference date",
    )
    hdr["TIMEUNIT"] = ("d", "time unit for TIME, TSTART and TSTOP")
    hdr["TELAPSE "] = (
        hdr.get("TSTOP", 0) - hdr.get("TSTART", 0),
        "[d] TSTOP - TSTART",
    )

    # Updating card comment to be more explicit
    if hasattr(hdr, "DATE"):
        hdr["DATE"] = (hdr["DATE"], "FFI cube creation date")
    if hasattr(hdr, "TSTART"):
        hdr["TSTART"] = (
            hdr["TSTART"],
            "observation start time in TJD of first FFI",
        )
    if hasattr(hdr, "TSTOP"):
        hdr["TSTOP"] = (
            hdr["TSTOP"],
            "observation stop time in TJD of last FFI",
        )
    if hasattr(hdr, "DATE-OBS"):
        hdr["DATE-OBS"] = (
            hdr["DATE-OBS"],
            "TSTART as UTC calendar date of first FFI",
        )
    if hasattr(hdr, "DATE-END"):
        hdr["DATE-END"] = (
            hdr["DATE-END"],
            "TSTOP as UTC calendar date of last FFI",
        )
    if hasattr(hdr, "FFIINDEX"):
        hdr["FFIINDEX"] = (
            hdr["FFIINDEX"],
            "number of FFI cadence interval of first FFI",
        )

    # These are all the things in the TESS pipeline tpfs about the object that we can't fill
    hdr["OBJECT"] = ("", "string version of target id ")
    hdr["TCID"] = (0, "unique tess target identifier")
    hdr["PXTABLE"] = (0, "pixel table id")
    hdr["PMRA"] = (0.0, "[mas/yr] RA proper motion")
    hdr["PMDEC"] = (0.0, "[mas/yr] Dec proper motion")
    hdr["PMTOTAL"] = (0.0, "[mas/yr] total proper motion")
    hdr["TESSMAG"] = (0.0, "[mag] TESS magnitude")
    hdr["TEFF"] = (0.0, "[K] Effective temperature")
    hdr["LOGG"] = (0.0, "[cm/s2] log10 surface gravity")
    hdr["MH"] = (0.0, "[log10([M/H])] metallicity")
    hdr["RADIUS"] = (0.0, "[solar radii] stellar radius")
    hdr["TICVER"] = (0, "TICVER")
    hdr["TICID"] = (None, "unique tess target identifier")

    delete_kwds_wildcards = [
        "SC_*",
        "RMS*",
        "A_*",
        "AP_*",
        "B_*",
        "BP*",
        "GAIN*",
        "TESS_*",
        "CD*",
        "CT*",
        "CRPIX*",
        "CRVAL*",
        "MJD*",
    ]
    delete_kwds = [
        "COMMENT",
        "FILTER",
        "TIME",
        "EXPTIME",
        "ACS_MODE",
        "DEC_TARG",
        "FLXWIN",
        "RA_TARG",
        "CCDNUM",
        "CAMNUM",
        "WCSGDF",
        "UNITS",
        "CADENCE",
        "SCIPIXS",
        "INT_TIME",
        "PIX_CAT",
        "REQUANT",
        "DIFF_HUF",
        "PRIM_HUF",
        "QUAL_BIT",
        "SPM",
        "STARTTJD",
        "ENDTJD",
        "CRM",
        "TJD_ZERO",
        "CRM_N",
        "ORBIT_ID",
        "MIDTJD",
    ]

    # Bulk removal with wildcards. Most of these should only live in EXT 1 header.
    for kwd in delete_kwds_wildcards:
        if kwd in hdr:
            del hdr[kwd]

    # Removal of specific kwds not relevant for cutouts.
    # Most likely these describe a single FFI, and not
    # the whole cube, which is misleading because we are
    # working with entire stacks of FFIs. Other keywords
    # are analogs to ones that have already been added
    # to the primary header in the lines above.
    for kwd in delete_kwds:
        if kwd in hdr:
            del hdr[kwd]
    return hdu0


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


def _extract_all_WCS(hdu):
    """Extract all the WCSs from a TableHDU from the TESS cube"""
    wcss = []
    for idx in np.arange(len(hdu.data["CRPIX1"])):
        try:
            wcs_hdu = fits.PrimaryHDU()
            for attr in WCS_ATTRS(hdu):
                wcs_hdu.header[attr] = hdu[0].data[attr][idx]
            wcs_hdu.header["WCSAXES"] = int(wcs_hdu.header["WCSAXES"])
            wcs_hdu.header["WCSAXESP"] = int(wcs_hdu.header["WCSAXESP"])
            wcss.append(WCS(wcs_hdu.header))
        except:
            wcss.append(None)
    return wcss


@functools.lru_cache(maxsize=4)
def get_FFI(filename):
    """Downloads a Full Frame Image from AWS"""
    return _sync_call(_async_get_FFI, filename)


async def _async_get_FFI(filename):
    """Downloads a Full Frame Image from AWS"""
    sector, camera, ccd = filename.split("-")[1:4]
    year, day = filename[4:8], filename[8:11]
    object_key = f"tess/public/ffi/{sector}/{year}/{day}/{camera}-{ccd}/{filename}"
    async with get_session().create_client(
        "s3", config=Config(signature_version=UNSIGNED)
    ) as s3:
        response = await s3.get_object(Bucket=BUCKET_NAME, Key=object_key)
        data_bytes = await response["Body"].read()
    return fits.open(
        BytesIO(data_bytes),
        ignore_missing_simple=True,
        lazy_load_hdus=False,
    )
