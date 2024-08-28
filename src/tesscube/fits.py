"""Functions to fix fits formats"""

import numpy as np
from astropy.io import fits
from astropy.time import Time

from . import __version__, log


def get_header_dict(cube):
    delete_kwds = [
        "SIMPLE",
        "BITPIX",
        "NAXIS",
        "EXTEND",
        "NEXTEND",
        "EXTNAME",
        "EXTVER",
        "COMMENT",
        "FILTER",
        "NAXIS1",
        "NAXIS2",
        "NAXIS3",
        "NAXIS4",
        "XTENSION",
        "PCOUNT",
        "GCOUNT",
        "TFIELDS",
    ]

    hdu = cube.last_hdu
    hdr = cube.last_hdu.header
    header_dict = {
        k: fits.Card(
            k,
            hdu.data[k][0]
            if isinstance(hdu.data[k][0], str)
            else np.nanmedian(hdu.data[k]),
            c,
        )
        if (k in hdu.data.columns.names)
        else fits.Card(k, hdr[k], c)
        for k, c in zip(hdr.keys(), hdr.comments)
        if (
            (k not in delete_kwds)
            and (not k.startswith("TTYPE") and (not k.startswith("TFORM")))
        )
    }
    hdr = cube.primary_hdu.header
    header_dict.update(
        {
            k: fits.Card(k, hdr[k], c)
            for k, c in zip(hdr.keys(), hdr.comments)
            if (k not in delete_kwds)
        }
    )
    header_dict["CREATOR"] = fits.Card(
        "CREATOR", "tesscube", "software used to produce this file"
    )
    header_dict["PROCVER"] = fits.Card("PROCVER", __version__, "software version")
    header_dict["TSTART"] = fits.Card(
        "TSTART",
        cube.last_hdu.data["TSTART"][0],
        "observation start time in TJD of first FFI",
    )
    header_dict["TSTOP"] = fits.Card(
        "TSTOP",
        cube.last_hdu.data["TSTOP"][-1],
        "observation stop time in TJD of last FFI",
    )
    header_dict["DATE-OBS"] = fits.Card(
        "DATE-OBS",
        Time(cube.last_hdu.data["TSTART"][0] + 2457000, format="jd").isot,
        "TSTART as UTC calendar date",
    )
    header_dict["DATE-END"] = fits.Card(
        "DATE-END",
        Time(cube.last_hdu.data["TSTOP"][-1] + 2457000, format="jd").isot,
        "TSTOP as UTC calendar date",
    )
    return header_dict


def get_output_primary_hdu(cube):
    keys = [
        "SIMDATA",
        "ORIGIN",
        "DATE",
        "TSTART",
        "TSTOP",
        "DATE-OBS",
        "DATE-END",
        "CREATOR",
        "PROCVER",
        "FILEVER",
        "TIMVERSN",
        "TELESCOP",
        "INSTRUME",
        "DATA_REL",
        "OBJECT",
        "TICID",
        "SECTOR",
        "CAMERA",
        "CCD",
        "PXTABLE",
        "RADESYS",
        "RA_OBJ",
        "DEC_OBJ",
        "EQUINOX",
        "PMRA",
        "PMDEC",
        "PMTOTAL",
        "TESSMAG",
        "TEFF",
        "LOGG",
        "MH",
        "RADIUS",
        "TICVER",
        "CRMITEN",
        "CRBLKSZ",
        "CRSPOC",
    ]
    header_dict = cube.header_dict
    return fits.PrimaryHDU(header=fits.Header(cards=[header_dict[key] for key in keys]))


def get_output_first_extention_header(cube):
    keys = [
        "SIMDATA",
        "TELESCOP",
        "INSTRUME",
        "OBJECT",
        "TICID",
        "RADESYS",
        "RA_OBJ",
        "DEC_OBJ",
        "EQUINOX",
        "EXPOSURE",
        "TIMEREF",
        "TASSIGN",
        "TIMESYS",
        "BJDREFI",
        "BJDREFF",
        "TIMEUNIT",
        "TELAPSE",
        "LIVETIME",
        "TSTART",
        "TSTOP",
        "DATE-OBS",
        "DATE-END",
        "DEADC",
        "TIMEPIXR",
        "TIERRELA",
        "INT_TIME",
        "READTIME",
        "FRAMETIM",
        "NUM_FRM",
        "TIMEDEL",
        "BACKAPP",
        "DEADAPP",
        "VIGNAPP",
        "GAINA",
        "GAINB",
        "GAINC",
        "GAIND",
        "READNOIA",
        "READNOIB",
        "READNOIC",
        "READNOID",
        "NREADOUT",
        "FXDOFF",
        "MEANBLCA",
        "MEANBLCB",
        "MEANBLCC",
        "MEANBLCD",
    ]
    header_dict = cube.header_dict
    return fits.Header(cards=[header_dict[key] for key in keys])


def get_output_second_extension_header(cube):
    keys = [
        "SIMDATA",
        "TELESCOP",
        "INSTRUME",
        "OBJECT",
        "TICID",
        "RADESYS",
        "RA_OBJ",
        "DEC_OBJ",
        "EQUINOX",
    ]
    header_dict = cube.header_dict
    return fits.Header(cards=[header_dict[key] for key in keys])


def get_wcs_header_by_extension(wcs_header, ext=5):
    """
    Converts WCS header into the equivalent for each extension
    """
    return fits.Header(
        [
            fits.Card(f"WCAX{ext}", wcs_header["WCSAXES"], "number of WCS axes"),
            fits.Card(
                f"1CTYP{ext}", wcs_header["CTYPE1"], "right ascension coordinate type"
            ),
            fits.Card(
                f"2CTYP{ext}", wcs_header["CTYPE2"], "declination coordinate type"
            ),
            fits.Card(
                f"1CRPX{ext}",
                wcs_header["CRPIX1"],
                "[pixel] reference pixel along image axis 1",
            ),
            fits.Card(
                f"2CRPX{ext}",
                wcs_header["CRPIX2"],
                "[pixel] reference pixel along image axis 2",
            ),
            fits.Card(
                f"1CRVL{ext}",
                wcs_header["CRVAL1"],
                "[deg] right ascension at reference pixel",
            ),
            fits.Card(
                f"2CRVL{ext}",
                wcs_header["CRVAL2"],
                "[deg] declination at reference pixel",
            ),
            fits.Card(
                f"1CUNI{ext}", wcs_header["CUNIT1"], "physical unit in column dimension"
            ),
            fits.Card(
                f"2CUNI{ext}", wcs_header["CUNIT2"], "physical unit in row dimension"
            ),
            fits.Card(
                f"1CDLT{ext}", wcs_header["CDELT1"], "[deg] pixel scale in RA dimension"
            ),
            fits.Card(
                f"2CDLT{ext}",
                wcs_header["CDELT1"],
                "[deg] pixel scale in DEC dimension",
            ),
            fits.Card(
                f"11PC{ext}",
                wcs_header["PC1_1"],
                "Coordinate transformation matrix element",
            ),
            fits.Card(
                f"12PC{ext}",
                wcs_header["PC1_2"],
                "Coordinate transformation matrix element",
            ),
            fits.Card(
                f"21PC{ext}",
                wcs_header["PC2_1"],
                "Coordinate transformation matrix element",
            ),
            fits.Card(
                f"22PC{ext}",
                wcs_header["PC2_2"],
                "Coordinate transformation matrix element",
            ),
            fits.Card(f"WCSN{ext}P", "PHYSICAL", "table column WCS name"),
            fits.Card(f"WCAX{ext}P", 2, "table column physical WCS dimensions"),
            fits.Card(
                f"1CTY{ext}P", "RAWX", "table column physical WCS axis 1 type, CCD col"
            ),
            fits.Card(
                f"2CTY{ext}P", "RAWY", "table column physical WCS axis 2 type, CCD row"
            ),
            fits.Card(f"1CUN{ext}P", "PIXEL", "table column physical WCS axis 1 unit"),
            fits.Card(f"2CUN{ext}P", "PIXEL", "table column physical WCS axis 2 unit"),
            fits.Card(
                f"1CRV{ext}P", wcs_header["CRVAL1P"], "value at reference CCD column"
            ),
            fits.Card(
                f"2CRV{ext}P", wcs_header["CRVAL2P"], "value at reference CCD row"
            ),
            fits.Card(f"1CDL{ext}P", 1.0, "table column physical WCS a1 step"),
            fits.Card(f"2CDL{ext}P", 1.0, "table column physical WCS a2 step"),
            fits.Card(f"1CRP{ext}P", 1, "table column physical WCS a1 reference"),
            fits.Card(f"2CRP{ext}P", 1, "table column physical WCS a2 reference"),
        ]
    )


def _fix_primary_hdu(hdu):
    """Fixes the primary HDU header to mimic TESS official products.
    This is largely taken from MASTs astrocut package under the BSD license.
    """
    log.debug("Fixing primary HDU.")
    hdu0 = hdu.copy()
    hdr = hdu0.header
    hdr["CREATOR"] = ("tesscube", "software used to produce this file")
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

    hdr["TSTART"] = (hdr["TSTART"], "observation start time in TJD of first FFI")
    hdr["TSTOP"] = (hdr["TSTOP"], "observation stop time in TJD of last FFI")
    hdr["DATE-OBS"] = (hdr["DATE-OBS"], "TSTART as UTC calendar date of first FFI")
    hdr["DATE-END"] = (hdr["DATE-END"], "TSTOP as UTC calendar date of last FFI")
    hdr["FFIINDEX"] = (hdr["FFIINDEX"], "number of FFI cadence interval of first FFI")

    # These are all the things in the TESS pipeline tpfs about the object that we can't fill
    hdr["OBJECT"] = ("", "string version of target id ")
    hdr["RA_OBJ"] = (0, "[deg] right ascension")
    hdr["DEC_OBJ"] = (0, "[deg] declination")
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
        "NAXIS*" "SC_*",
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
        "FLXWIN",
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
