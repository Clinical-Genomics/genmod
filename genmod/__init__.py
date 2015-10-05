from __future__ import (print_function, absolute_import)
from logging import getLogger
from pkg_resources import get_distribution

logger = getLogger(__name__)

__version__ = get_distribution("genmod").version
