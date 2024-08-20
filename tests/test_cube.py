import pytest

from tesscube import TESSCube
from tesscube.query import get_last_hdu

# This test needs to be run first!
@pytest.mark.remote_data
def test_caching():
    # Check that we are doing caching correctly
    r = TESSCube(sector=1, camera=1, ccd=1)
    r.last_hdu  # First call
    r.last_hdu  # Second call

    # Check if the cache has one entry for this function
    assert get_last_hdu.cache_info().hits == 1
