from importlib.metadata import version
from logging import getLogger

logger = getLogger(__name__)

__version__ = version("genmod")
