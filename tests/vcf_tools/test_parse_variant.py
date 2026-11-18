from genmod.vcf_tools import get_variant_id, parse_annotations

ANN_COLUMNS = ["Allele", "Annotation", "Annotation_Impact", "Gene_Name", "Gene_ID", "Distance"]


class TestParseAnnotations:
    def test_one_annotation_per_allele(self):
        # GIVEN a SnpEff ANN value of a multiallelic site, one annotation per ALT allele
        ann = "G|missense_variant|MODERATE|ADK|ENSG01|,T|stop_gained|HIGH|ADK|ENSG01|"

        # WHEN parsing it
        annotations = parse_annotations(ann, ANN_COLUMNS)

        # THEN there is one dictionary per annotation, in order, with empty fields kept
        assert [a["Allele"] for a in annotations] == ["G", "T"]
        assert annotations[1]["Annotation_Impact"] == "HIGH"
        assert annotations[0]["Distance"] == ""

    def test_non_reference_allele(self):
        # GIVEN an annotation against a non-standard reference (cancer), Allele is ALT-REF
        annotations = parse_annotations("G-C|missense_variant|MODERATE|ADK|ENSG01|", ANN_COLUMNS)

        # THEN the allele is kept as written
        assert annotations[0]["Allele"] == "G-C"

    def test_empty(self):
        # A variant without the annotation (None) or a flag without value ([]) has none
        assert parse_annotations(None, ANN_COLUMNS) == []
        assert parse_annotations([], ANN_COLUMNS) == []
        assert parse_annotations("", ANN_COLUMNS) == []


class TestGetVariantId:
    def test_get_variant_id(self):
        variant = {"CHROM": "1", "POS": "10", "REF": "A", "ALT": "G"}
        assert get_variant_id(variant) == "1_10_A_G"

    def test_get_variant_id_sv_ins(self):
        variant = {"CHROM": "1", "POS": "10", "REF": "N", "ALT": "<INS>"}
        assert get_variant_id(variant) == "1_10_N_INS"

    def test_get_variant_id_sv_end(self):
        variant = {
            "CHROM": "1",
            "POS": "10",
            "REF": "N",
            "ALT": "<DEL>",
            "info_dict": {"END": "20"},
        }
        assert get_variant_id(variant) == "1_10_N_DEL_END20"

    def test_get_variant_id_sv_svlen(self):
        variant = {
            "CHROM": "1",
            "POS": "10",
            "REF": "N",
            "ALT": "<DEL>",
            "info_dict": {"SVLEN": "-10"},
        }
        assert get_variant_id(variant) == "1_10_N_DEL_SVLEN-10"

    def test_get_variant_id_sv_end_over_svlen(self):
        variant = {
            "CHROM": "1",
            "POS": "10",
            "REF": "N",
            "ALT": "<DEL>",
            "info_dict": {"END": "20", "SVLEN": "-10"},
        }
        assert get_variant_id(variant) == "1_10_N_DEL_END20"

    def test_get_variant_id_sv_dup_tandem(self):
        variant = {"CHROM": "1", "POS": "10", "REF": "N", "ALT": "<DUP:TANDEM>"}
        assert get_variant_id(variant) == "1_10_N_DUPTANDEM"

    def test_get_variant_id_sv_bdn(self):
        variant = {"CHROM": "1", "POS": "10", "REF": "A", "ALT": "T[6:134717462["}
        assert get_variant_id(variant) == "1_10_A_T6134717462"
