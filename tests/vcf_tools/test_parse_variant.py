from genmod.vcf_tools import get_variant_id

class TestGetVariantId:
    
    def test_get_variant_id(self):
        variant = {
            'CHROM': '1',
            'POS': '10',
            'REF': 'A',
            'ALT': 'G'
        }
        assert get_variant_id(variant) == "1_10_A_G"
    
    def test_get_variant_id_sv_ins(self):
        variant = {
            'CHROM': '1',
            'POS': '10',
            'REF': 'N',
            'ALT': '<INS>'
        }
        assert get_variant_id(variant) == "1_10_N_INS"
    
    def test_get_variant_id_sv_dup_tandem(self):
        variant = {
            'CHROM': '1',
            'POS': '10',
            'REF': 'N',
            'ALT': '<DUP:TANDEM>'
        }
        assert get_variant_id(variant) == "1_10_N_DUPTANDEM"
    
    def test_get_variant_id_sv_bdn(self):
        variant = {
            'CHROM': '1',
            'POS': '10',
            'REF': 'A',
            'ALT': 'T[6:134717462['
        }
        assert get_variant_id(variant) == "1_10_A_T6134717462"