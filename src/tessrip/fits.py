"""Functions to fix fits formats"""

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
