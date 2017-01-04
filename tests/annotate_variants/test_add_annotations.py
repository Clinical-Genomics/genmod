from genmod.annotate_variants.add_annotations import (add_regions, add_exac, 
    add_exac_max, add_thousandg, add_thousandg_max, add_spidex, add_cadd, 
    add_cadd_raw, add_cosmic)

def test_add_regions(header):
    key = 'Annotation'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_regions(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_exac(header):
    key = 'EXACAF'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_exac(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_exac_max(header):
    key = 'EXAC_MAX_AF'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_exac_max(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_thousand_g(header):
    key = '1000GAF'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_thousandg(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_thousand_g_max(header):
    key = '1000G_MAX_AF'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_thousandg_max(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_spidex(header):
    key = 'SPIDEX'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_spidex(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_cadd(header):
    key = 'CADD'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_cadd(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_cadd_raw(header):
    key = 'CADD_raw'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_cadd_raw(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict

def test_add_cosmic(header):
    key = 'COSMIC'
    # GIVEN a header without key
    assert key not in header.info_dict
    # WHEN adding annotation
    add_cosmic(header)
    # THEN the key should have been added to the header
    assert key in header.info_dict
