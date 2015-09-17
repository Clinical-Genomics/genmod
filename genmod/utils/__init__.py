# -*- coding: utf-8 -*-
<<<<<<< HEAD
from __future__ import print_function
from genmod.utils.load_annotation import load_annotations
from genmod.utils.is_number import is_number
from genmod.utils.plugin_reader import collectKeys
from genmod.utils.sort_variants import sort_variants
from genmod.utils.print_variants import print_variants, print_variant_for_sorting
from genmod.utils.print_headers import print_headers
from genmod.utils.add_metadata import add_metadata
from genmod.utils.pair_generator import PairGenerator
=======
from __future__ import (print_function, absolute_import)

INTERESTING_SO_TERMS = set(
            [
                u'transcript_ablation',
                u'splice_donor_variant',
                u'splice_acceptor_variant',
                u'stop_gained',
                u'frameshift_variant',
                u'stop_lost',
                u'initiator_codon_variant',
                u'transcript_amplification',
                u'inframe_insertion',
                u'inframe_deletion',
                u'missense_variant',
                u'splice_region_variant',
                u'incomplete_terminal_codon_variant',
                u'stop_retained_variant',
                u'synonymous_variant',
                u'coding_sequence_variant',
                u'mature_miRNA_variant',
                u'5_prime_UTR_variant',
                u'3_prime_UTR_variant',
                u'non_coding_exon_variant',
                u'non_coding_transcript_exon_variant',
                u'non_coding_transcript_variant',
                u'nc_transcript_variant',
                u'intron_variant',
                u'NMD_transcript_variant',
                u'non_coding_transcript_variant',
            ]
)

EXONIC_SO_TERMS = set(
            [
                u'transcript_ablation',
                u'splice_donor_variant',
                u'splice_acceptor_variant',
                u'stop_gained',
                u'frameshift_variant',
                u'stop_lost',
                u'initiator_codon_variant',
                u'transcript_amplification',
                u'inframe_insertion',
                u'inframe_deletion',
                u'missense_variant',
                u'splice_region_variant',
                u'incomplete_terminal_codon_variant',
                u'stop_retained_variant',
                u'synonymous_variant',
                u'coding_sequence_variant',
            ]
)

from .is_number import is_number
from .check_individuals import check_individuals
from .pair_generator import generate_pairs
from .get_features import (get_annotation, check_vep_annotation, check_exonic)
from .get_batches import get_batches
from .get_priority import (get_chromosome_priority, get_rank_score)
from .variant_printer import VariantPrinter
>>>>>>> feature/fix_compounds_single
