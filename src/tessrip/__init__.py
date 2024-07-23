__version__ = "0.1.0dev"
# Standard library
import os  # noqa

PACKAGEDIR = os.path.abspath(os.path.dirname(__file__))

HDR_SIZE = 2880  # bytes
BYTES_PER_PIX = 4  # float32
MAX_CONCURRENT_DOWNLOADS = 10
DATA_OFFSET = HDR_SIZE * 2
BUCKET_NAME = "stpubdata"

from .query import Rip  # noqa
from .utils import get_FFI  # noqa
