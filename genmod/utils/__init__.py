# -*- coding: utf-8 -*-
from __future__ import (print_function, absolute_import)

INTERESTING_SO_TERMS = set(
            [
                u'transcript_ablation',
                u'splice_donor_variant',
                u'splice_acceptor_variant',
                u'stop_gained',
                u'start_lost',
                u'frameshift_variant',
                u'stop_lost',
                u'initiator_codon_variant',
                u'transcript_amplification',
                u'inframe_insertion',
                u'inframe_deletion',
                u'protein_altering_variant',
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
                u'upstream_gene_variant',
                u'downstream_gene_variant',
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
                u'start_lost',
                u'initiator_codon_variant',
                u'transcript_amplification',
                u'inframe_insertion',
                u'inframe_deletion',
                u'missense_variant',
                u'protein_altering_variant',
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
