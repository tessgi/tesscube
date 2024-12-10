import pytest

from tesscube import TESSCube


# This test needs to be run first!
@pytest.mark.remote_data
def test_caching():
    # Check that we are doing caching correctly
    TESSCube(sector=1, camera=1, ccd=1)
