from genmod.annotate_regions import parse_annotations

def test_parse_annotations(ccds_file):
    """
    Test parse_annotations
    """
    gene_trees, exon_trees = parse_annotations(infile=ccds_file, annotation_type="ccds")
    assert type(gene_trees) == type({})
    assert set(gene_trees.keys()) == set(['1', '2'])
