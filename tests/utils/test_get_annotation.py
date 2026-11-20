from genmod.utils import SNPEFF, VEP, get_annotation

VEP_COLUMNS = ["Allele", "Consequence", "Gene"]
ANN_COLUMNS = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID"]


def get_variant(chrom="1", pos="1", ref="A", alt="G", annotation=["ADK"]):
    """
    Return a variant dictionary
    """
    variant_id = "_".join([chrom, pos, ref, alt])
    variant = {
        "CHROM": chrom,
        "POS": pos,
        "INFO": "Annotation={0}".format(annotation),
        "info_dict": {
            "Annotation": ",".join(annotation),
        },
        "variant_id": variant_id,
    }

    return variant


def get_vep_variant(chrom="1", pos="1", ref="A", alt="G", annotation="ADK"):
    """
    Return a variant dictionary
    """
    variant_id = "_".join([chrom, pos, ref, alt])
    variant = {
        "CHROM": chrom,
        "POS": pos,
        "info_dict": {"CSQ": "G|transcript_ablation|{0}".format(annotation)},
        "variant_id": variant_id,
    }
    return variant


def test_empty():
    """Test if get_features behave as suspected"""
    variant = {"CHROM": "1", "POS": "12", "ALT": "A"}

    assert get_annotation(variant) == set()


def test_simple():
    """Test if get_annotation behave as suspected"""
    variant = get_variant()
    assert get_annotation(variant) == set(["ADK"])


def test_double_region():
    """Test if get_annotation behave as suspected"""
    variant = get_variant(annotation=["ADK", "DDD"])
    assert get_annotation(variant) == set(["ADK", "DDD"])


def test_get_vep_region():
    """docstring for test_get_vep_region"""
    variant = get_vep_variant()
    assert get_annotation(variant, source=VEP, columns=VEP_COLUMNS) == set(["ADK"])


def test_get_snpeff_region():
    # GIVEN a variant with a SnpEff ANN annotation
    variant = {
        "variant_id": "1_1_A_G",
        "info_dict": {"ANN": "G|missense_variant|MODERATE|ADK|ENSG00000156110"},
    }

    # WHEN getting the annotation with the SnpEff source
    annotation = get_annotation(variant, source=SNPEFF, columns=ANN_COLUMNS)

    # THEN the gene id is returned
    assert annotation == {"ENSG00000156110"}


def test_source_annotation_missing_from_variant():
    # GIVEN a variant without the ANN field
    variant = {"variant_id": "1_1_A_G", "info_dict": {"MQ": "1"}}

    # WHEN getting the annotation with the SnpEff source
    annotation = get_annotation(variant, source=SNPEFF, columns=ANN_COLUMNS)

    # THEN the variant has no annotation, instead of a KeyError
    assert annotation == set()
