# Here is the code for annotating variants with what genetic region they
# belong to.
# In general this includes building structures for region annotation,
# annotating variants based on these structures etc

from __future__ import absolute_import

from .check_overlap import check_overlap
from .load_annotation import load_annotations
from .parse_annotations import parse_annotations
from .get_features import get_genes, check_exonic