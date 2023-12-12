import numpy as np
import pytest
from astropy.io import fits

from tessrip import Rip
from tessrip.query import _last_hdulist


# This test needs to be run first!
@pytest.mark.remote_data
def test_caching():
    # Check that we are doing caching correctly
    r = Rip(sector=1, camera=1, ccd=1)
    r.last_hdulist  # First call
    r.last_hdulist  # Second call

    # Check if the cache has one entry for this function
    assert _last_hdulist.cache_info().hits == 1


@pytest.mark.remote_data
def test_rip():
    """Test that we can initialize a Rip object"""
    for ntimes, sector in zip([1282, 3464], [1, 28]):
        r = Rip(sector=sector, camera=1, ccd=1)
        flux, flux_err = r.get_flux(corner=(200, 201), shape=(2, 3))
        for ar in [flux, flux_err]:
            assert isinstance(ar, fits.Column)
            assert ar.array.ndim == 3
            assert ar.array.shape == (ntimes, 2, 3)
        flux, flux_err = r.get_flux(corner=(200, 201), shape=(2, 3), frame_range=(0, 1))
        for ar in [flux, flux_err]:
            assert isinstance(ar, fits.Column)
            assert ar.array.ndim == 3
            assert ar.array.shape == (1, 2, 3)
        assert hasattr(r, "center")
        assert hasattr(r, "primary_hdu")
        assert hasattr(r, "last_hdulist")
        assert isinstance(r.primary_hdu, fits.PrimaryHDU)
        assert isinstance(r.last_hdulist, fits.HDUList)
        for attr in ["time", "timecorr", "wcs", "cadence_number"]:
            assert hasattr(r, attr)
