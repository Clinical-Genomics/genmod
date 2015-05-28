# -*- coding: utf-8 -*-
from __future__ import (print_function, unicode_literals, absolute_import)

INTERESTING_SO_TERMS = set(
            [
                'transcript_ablation',
                'splice_donor_variant',
                'splice_acceptor_variant',
                'stop_gained',
                'frameshift_variant',
                'stop_lost',
                'initiator_codon_variant',
                'transcript_amplification',
                'inframe_insertion',
                'inframe_deletion',
                'missense_variant',
                'splice_region_variant',
                'incomplete_terminal_codon_variant',
                'stop_retained_variant',
                'synonymous_variant',
                'coding_sequence_variant',
                'mature_miRNA_variant',
                '5_prime_UTR_variant',
                '3_prime_UTR_variant',
                'non_coding_exon_variant',
                'non_coding_transcript_exon_variant',
                'non_coding_transcript_variant',
                'nc_transcript_variant',
                'intron_variant',
                'NMD_transcript_variant',
                'non_coding_transcript_variant',
            ]
)

EXONIC_SO_TERMS = set(
            [
                'transcript_ablation',
                'splice_donor_variant',
                'splice_acceptor_variant',
                'stop_gained',
                'frameshift_variant',
                'stop_lost',
                'initiator_codon_variant',
                'transcript_amplification',
                'inframe_insertion',
                'inframe_deletion',
                'missense_variant',
                'splice_region_variant',
                'incomplete_terminal_codon_variant',
                'stop_retained_variant',
                'synonymous_variant',
                'coding_sequence_variant',
            ]
)

from .load_annotation import load_annotations
from .is_number import is_number
from genmod.utils.plugin_reader import collectKeys
from .sort_variants import sort_variants
from .print_variants import print_variants
from .print_headers import print_headers
from .add_metadata import add_metadata
from .pair_generator import generate_pairs

from .get_features import (get_annotation, check_vep_annotation)
from .get_batches import get_batches

