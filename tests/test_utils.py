from astropy.io import fits

from tessrip.utils import _fix_primary_hdu, convert_coordinates_to_runs


def test_fix_primary():
    """Test that primary headers can be fixed"""
    hdu = fits.PrimaryHDU()
    hdu = _fix_primary_hdu(hdu)

    assert "BJDREFI" in hdu.header
    assert "FFI_TYPE" in hdu.header
    assert "TICID" in hdu.header


def test_runs():
    """Test if we can convert individual coordinates to row wise runs."""
    coords = [(2, 2), (2, 3), (2, 4)]
    runs = convert_coordinates_to_runs(coords)
    assert runs == [{"start_column": 2, "ncolumns": 3, "row": 2}]
