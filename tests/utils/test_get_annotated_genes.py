from genmod.utils import SNPEFF, VEP, get_annotated_genes


def test_get_none():
    """
    Test to get wrong annotation terms
    """
    annotations = [{"Consequence": "Two&Different", "Gene": "ADK"}]

    # The result should be empty since the terms do not exist
    assert get_annotated_genes(annotations, VEP) == set()


def test_get_annotation():
    """
    Test to get true annotation terms
    """
    annotations = [{"Consequence": "transcript_ablation&splice_donor_variant", "Gene": "ADK"}]

    assert get_annotated_genes(annotations, VEP) == set(["ADK"])


def test_empty_gene_is_skipped():
    # GIVEN a VEP annotation with an interesting consequence but no gene
    annotations = [{"Consequence": "missense_variant", "Gene": ""}]

    # WHEN getting the genes
    genes = get_annotated_genes(annotations, VEP)

    # THEN no empty gene is added
    assert genes == set()


def test_snpeff_gene_ids():
    # GIVEN SnpEff annotations over two genes, one with a SnpEff-only SO term
    annotations = [
        {"Annotation": "missense_variant&splice_region_variant", "Gene_ID": "ENSG01"},
        {"Annotation": "conservative_inframe_deletion", "Gene_ID": "ENSG02"},
    ]

    # WHEN getting the genes with the SnpEff source
    genes = get_annotated_genes(annotations, SNPEFF)

    # THEN the gene ids of both annotations are returned
    assert genes == {"ENSG01", "ENSG02"}


def test_snpeff_intergenic_region_has_no_gene():
    # GIVEN a SnpEff intergenic annotation, where Gene_ID holds the two closest genes
    annotations = [{"Annotation": "intergenic_region", "Gene_ID": "ENSG01-ENSG02"}]

    # WHEN getting the genes with the SnpEff source
    genes = get_annotated_genes(annotations, SNPEFF)

    # THEN the variant does not belong to a gene
    assert genes == set()
