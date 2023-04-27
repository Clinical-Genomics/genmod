from __future__ import (print_function, absolute_import)
from logging import getLogger
from pkg_resources import get_distribution, DistributionNotFound
from os import environ

logger = getLogger(__name__)

try:
    __version__ = get_distribution("genmod").version
except DistributionNotFound as error:
    # Genmod is not installed (as the case for docker test suite)
    if 'GENMOD_DOCKER_TEST' not in environ.keys():
        raise error
    __version__ = '-1.-1.-1'
