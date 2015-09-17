from __future__ import (print_function, unicode_literals, absolute_import)
from logging import getLogger
from pkg_resources import require

logger = getLogger(__name__)

__version__ = require("genmod")[0].version

# from genmod.annotation_parser import AnnotationParser
# from genmod.variant_annotation import (check_genetic_models, VariantAnnotator)
# from .variant_scorer import VariantScorer
# from genmod.get_batches import get_batches
# from genmod.utils import (sort_variants, load_annotations, is_number,
#     collectKeys, print_variants, print_headers, add_metadata, PairGenerator)
# from genmod.errors import warning
#
# from genmod.models import (check_dominant, check_recessive, check_compounds,
#                         check_X_recessive, check_X_dominant, score_variants)
#
# from genmod.log import init_log
