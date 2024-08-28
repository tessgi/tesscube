"""Utilities to help work with cube data"""

from typing import List, Tuple
import asyncio
import sys

import numpy as np

from . import log


def validate_tuple(t: tuple):
    if not len(t) == 2:
        raise ValueError("Pass a tuple with length 2.")
    return (int(t[0]), int(t[1]))


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


def convert_coordinates_to_runs(
    coordinates: List[Tuple[int, int]],
) -> List[Tuple[int, int, int]]:
    """
    Converts a list of (row, column) coordinates to a list of (start_row, end_row, column) coordinates.
    """
    # Step 1: Sort the coordinates by row and then by column.

    # create unique coordinates
    coordinates = [
        (int(i.split(", ")[0]), int(i.split(", ")[1]))
        for i in np.unique([f"{r}, {c}" for r, c in coordinates])
    ]
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
