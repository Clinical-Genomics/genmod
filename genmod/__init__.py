from logging import getLogger

try:
    from importlib.metadata import version
except ImportError:  # Backport support for importlib metadata on Python 3.7
    from importlib_metadata import version

logger = getLogger(__name__)

__version__ = version("genmod")
