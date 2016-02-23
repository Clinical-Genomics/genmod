from genmod.annotate_regions.parse_annotations import gene_pred_parser

def test_build_ref_seq_genes(ref_seq_lines):
    """
    Test if the parser builds the ccds gene dict correct
    """
    gene_trees, exon_trees = gene_pred_parser(ref_seq_lines)
    assert set(gene_trees.keys()) == set(['1', '6', '11', '13'])
