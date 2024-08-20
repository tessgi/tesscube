from tesscube.utils import convert_coordinates_to_runs


def test_runs():
    """Test if we can convert individual coordinates to row wise runs."""
    coords = [(2, 2), (2, 3), (2, 4)]
    runs = convert_coordinates_to_runs(coords)
    assert runs == [{"start_column": 2, "ncolumns": 3, "row": 2}]
