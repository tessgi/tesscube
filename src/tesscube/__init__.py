__version__ = "1.0.2dev"
# Standard library
import os  # noqa
import tempfile

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

HDR_SIZE = 2880  # bytes
BYTES_PER_PIX = 4  # float32
MAX_CONCURRENT_DOWNLOADS = 10
DATA_OFFSET = HDR_SIZE * 2
BUCKET_NAME = "stpubdata"


import logging  # noqa: E402

# This library lets us have log messages with syntax highlighting
from rich.logging import RichHandler  # noqa: E402


def get_logger():
    """Configure and return a logger with RichHandler."""
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.WARN)

    # Add RichHandler
    rich_handler = RichHandler()
    rich_handler.setFormatter(
        logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    )

    logger.addHandler(rich_handler)
    return logger


# Make sure that we have a directory we can put package config in
def check_package_cache_dir(package_cache_dir):
    if not os.path.exists(package_cache_dir):
        try:
            os.makedirs(package_cache_dir)  # Try creating the package directory
        except OSError as e:
            raise PermissionError(
                f"Unable to create package cache directory at {package_cache_dir}."
            ) from e

    # Check if the directory is writable
    try:
        with tempfile.NamedTemporaryFile(dir=package_cache_dir, delete=False) as fp:
            fp.write(b"Hello world!")
            fp.close()
            with open(fp.name, mode="rb") as f:
                f.read()
        os.unlink(fp.name)
    except OSError as e:
        raise PermissionError(
            f"The directory {package_cache_dir} is not writable."
        ) from e
    return package_cache_dir


_package_cache_dir = check_package_cache_dir(
    os.path.join(os.path.expanduser("~"), ".tesscube")
)

log = get_logger()

# Make sure we can load the config
from .config import load_config  # noqa

load_config()
from .cube import TESSCube  # noqa
