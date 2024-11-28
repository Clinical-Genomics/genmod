# -*- coding: utf-8 -*-

from .check_individuals import check_individuals
from .get_batches import get_batches
from .get_features import check_vep_annotation, get_annotation
from .get_priority import get_chromosome_priority, get_rank_score
from .is_number import is_number
from .pair_generator import generate_pairs
from .variant_printer import VariantPrinter

INTERESTING_SO_TERMS = {"transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained",
                        "start_lost", "frameshift_variant", "stop_lost", "initiator_codon_variant",
                        "transcript_amplification", "inframe_insertion", "inframe_deletion", "protein_altering_variant",
                        "missense_variant", "splice_region_variant", "incomplete_terminal_codon_variant",
                        "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
                        "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "non_coding_exon_variant",
                        "non_coding_transcript_exon_variant", "non_coding_transcript_variant", "nc_transcript_variant",
                        "intron_variant", "NMD_transcript_variant", "non_coding_transcript_variant",
                        "upstream_gene_variant", "downstream_gene_variant"}

EXONIC_SO_TERMS = {"transcript_ablation", "splice_donor_variant", "splice_acceptor_variant", "stop_gained",
                   "frameshift_variant", "stop_lost", "start_lost", "initiator_codon_variant",
                   "transcript_amplification", "inframe_insertion", "inframe_deletion", "missense_variant",
                   "protein_altering_variant", "splice_region_variant", "incomplete_terminal_codon_variant",
                   "stop_retained_variant", "synonymous_variant", "coding_sequence_variant"}
