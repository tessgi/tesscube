__version__ = "0.1.0dev"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

# Standard library
import logging  # noqa: E402

# This library lets us have log messages with syntax highlighting
from rich.logging import RichHandler  # noqa: E402

log = logging.getLogger("tessrip")
log.addHandler(RichHandler(markup=True))
log.setLevel("INFO")

HDR_SIZE = 2880  # bytes
BYTES_PER_PIX = 4  # float32
MAX_CONCURRENT_DOWNLOADS = 10
DATA_OFFSET = HDR_SIZE * 2
BUCKET_NAME = "stpubdata"

from .query import Rip
from .utils import get_FFI
