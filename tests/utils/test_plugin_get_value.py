from genmod.plugin import Plugin

csq_format = [
    "Allele",
    "Gene",
    "Feature",
    "Feature_type",
    "Consequence",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "SYMBOL",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "BIOTYPE",
    "TSL",
    "CCDS",
    "ENSP",
    "SWISSPROT",
    "TREMBL",
    "UNIPARC",
    "SIFT",
    "PolyPhen",
    "EXON",
    "INTRON",
    "DOMAINS",
    "HGVSc",
    "HGVSp",
    "MOTIF_NAME",
    "MOTIF_POS",
    "HIGH_INF_POS",
    "MOTIF_SCORE_CHANGE",
]

vcf_header = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "ADM1003A3"]


def get_variant_line(
    chrom="2",
    pos="191964633",
    db_id="rs7574865",
    ref="T",
    alt="G",
    qual="1258.77",
    filt="PASS",
    info="",
    form="GT:AD:GQ:PL",
    gt_calls={"ADM1003A3": "1/1:0,35:99:1287,105,0"},
):
    """
    Return a vcf formated variant line
    """
    if not info:
        info = (
            "1000GAF=0.744609;AC=2;AF=1.00;AN=2;CLNSIG=255,10;"
            "CSQ=G|ENSG00000138378|ENST00000358470|Transcript|intron_variant"
            "||||||||-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000351255"
            "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&C9JFG0_HUMAN"
            "&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000358470.4:c.274-23582A>C|||||"
            ",G|ENSG00000138378|ENST00000392320|Transcript|intron_variant|||||||"
            "|-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000376134"
            "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&"
            "C9JFG0_HUMAN&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000392320.2:c.274-23582A>C|||||"
            ";Clinical_db_gene_annotation=OMIM-141201,FullList;DB;DP=35;Dbsnp129LCAF=0.2493"
        )

    variant_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}".format(
        chrom, pos, db_id, ref, alt, qual, filt, info, form
    )
    for ind_id in gt_calls:
        variant_line = "{0}\t{1}".format(variant_line, gt_calls[ind_id])

    return variant_line


def get_variant_dict(
    chrom="2",
    pos="191964633",
    db_id="rs7574865",
    ref="T",
    alt="G",
    qual="1258.77",
    filt="PASS",
    info="",
    form="GT:AD:GQ:PL",
    gt_calls={"ADM1003A3": "1/1:0,35:99:1287,105,0"},
):
    """
    Return a vcf formated variant line
    """
    if not info:
        info = (
            "1000GAF=0.744609;AC=2;AF=1.00;AN=2;CLNSIG=255,10;"
            "CSQ=G|ENSG00000138378|ENST00000358470|Transcript|intron_variant"
            "||||||||-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000351255"
            "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&C9JFG0_HUMAN"
            "&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000358470.4:c.274-23582A>C|||||"
            ",G|ENSG00000138378|ENST00000392320|Transcript|intron_variant|||||||"
            "|-1|STAT4|HGNC|11365|protein_coding||CCDS2310.1|ENSP00000376134"
            "|STAT4_HUMAN|Q53S87_HUMAN&Q53RU2_HUMAN&E9PBE2_HUMAN&C9JM11_HUMAN&"
            "C9JFG0_HUMAN&B7SIX5_HUMAN|UPI00000015F2||||3/23||ENST00000392320.2:c.274-23582A>C|||||"
            ";Clinical_db_gene_annotation=OMIM-141201,FullList;DB;DP=35;Dbsnp129LCAF=0.2493"
        )

    info_dict = {}

    for info_entry in info.split(";"):
        splitted_info = info_entry.split("=")
        if len(splitted_info) > 1:
            info_dict[splitted_info[0]] = splitted_info[1]
        else:
            info_dict[splitted_info[0]] = []

    variant_dict = {
        "CHROM": chrom,
        "POS": pos,
        "ID": db_id,
        "REF": ref,
        "ALT": alt,
        "QUAL": qual,
        "FILTER": filt,
        "INFO": info,
        "FORMAT": form,
        "info_dict": info_dict,
    }

    for ind_id in gt_calls:
        variant_dict[ind_id] = gt_calls[ind_id]

    return variant_dict


def test_chrom():
    """Test to get the raw chromosome"""
    plugin = Plugin(name="Chrom", field="CHROM")
    test_chrom = "10"
    variant_line = get_variant_line(chrom=test_chrom)
    variant_dict = get_variant_dict(chrom=test_chrom)

    assert plugin.get_value(variant_line=variant_line) == test_chrom
    assert plugin.get_value(variant_dict=variant_dict) == test_chrom


def test_pos():
    """Test to get the raw chromosome"""
    plugin = Plugin(name="Pos", field="POS")
    test_pos = "1000"
    variant_line = get_variant_line(pos=test_pos)
    variant_dict = get_variant_dict(pos=test_pos)

    assert plugin.get_value(variant_line=variant_line) == test_pos
    assert plugin.get_value(variant_dict=variant_dict) == test_pos


def test_id():
    """Test to get the raw chromosome"""
    plugin = Plugin(name="ID", field="ID")
    db_id = "rs001"
    variant_line = get_variant_line(db_id=db_id)
    variant_dict = get_variant_dict(db_id=db_id)

    assert plugin.get_value(variant_line=variant_line) == db_id
    assert plugin.get_value(variant_dict=variant_dict) == db_id


def test_id_flag():
    """Test to get the raw chromosome"""
    plugin = Plugin(name="ID", field="ID")
    db_id = "rs001"
    variant_line = get_variant_line(db_id=db_id)
    variant_dict = get_variant_dict(db_id=db_id)

    assert plugin.get_value(variant_line=variant_line) == db_id
    assert plugin.get_value(variant_dict=variant_dict) == db_id


def test_multiple_id_no_rule():
    """Test to get the raw chromosome"""
    plugin = Plugin(name="ID", field="ID", data_type="flag")
    db_id_1 = "rs001"
    db_id_2 = "rs002"
    db_id = "{0};{1}".format(db_id_1, db_id_2)
    variant_line = get_variant_line(db_id=db_id)
    variant_dict = get_variant_dict(db_id=db_id)

    assert plugin.get_value(variant_line=variant_line) is True
    assert plugin.get_value(variant_dict=variant_dict) is True


# def test_non_existing():
#     """Test to get the raw chromosome"""
#     plugin = Plugin(name='non_existing', field='INFO', info_key="non", separators=[','])
#
#     variant_line = get_variant_line()
#     variant_dict = get_variant_dict()
#
#     assert plugin.get_entry(variant_line=variant_line) == []
#     assert plugin.get_entry(variant_dict=variant_dict) == []
#


def test_filter_min_rule():
    """Test to get the raw chromosome"""
    plugin = Plugin(
        name="Filter",
        field="FILTER",
        data_type="string",
        string_rules={"PASS": 2, "NOT_PASS": 1},
        record_rule="min",
    )
    filt = "PASS;NOT_PASS"

    variant_line = get_variant_line(filt=filt)
    variant_dict = get_variant_dict(filt=filt)

    assert plugin.get_value(variant_line=variant_line) == "NOT_PASS"
    assert plugin.get_value(variant_dict=variant_dict) == "NOT_PASS"


def test_filter_max_rule():
    """Test to get the raw chromosome"""
    plugin = Plugin(
        name="Filter",
        field="FILTER",
        data_type="string",
        string_rules={"PASS": 2, "NOT_PASS": 1},
        record_rule="max",
    )

    filt = "PASS;NOT_PASS"

    variant_line = get_variant_line(filt=filt)
    variant_dict = get_variant_dict(filt=filt)

    assert plugin.get_value(variant_line=variant_line) == "PASS"
    assert plugin.get_value(variant_dict=variant_dict) == "PASS"


# def test_1000G():
#     """Test to get the raw chromosome"""
#     plugin = Plugin(name='thousand_g', field='INFO', info_key="1000GAF", separators=[','])
#     assert plugin.get_value(variant) == "0.744609"
#
def test_1000G_float():
    """Test to get the raw chromosome"""
    plugin = Plugin(
        name="thousand_g",
        field="INFO",
        info_key="1000GAF",
        separators=[","],
        data_type="float",
    )

    test_value = 0.744609
    variant_line = get_variant_line()
    variant_dict = get_variant_dict()

    assert plugin.get_value(variant_line=variant_line) == test_value
    assert plugin.get_value(variant_dict=variant_dict) == test_value


def test_1000G_record_rule():
    """Test to get the raw chromosome"""
    plugin = Plugin(
        name="thousand_g",
        field="INFO",
        info_key="1000GAF",
        separators=[","],
        data_type="float",
        record_rule="min",
    )
    info = "1000GAF=0.744609,0.02;AC=2;AF=1.00;AN=2"

    variant_line = get_variant_line(info=info)
    variant_dict = get_variant_dict(info=info)

    assert plugin.get_value(variant_line=variant_line) == 0.02
    assert plugin.get_value(variant_dict=variant_dict) == 0.02


def test_flag():
    """Test to get the raw chromosome"""
    plugin = Plugin(
        name="DB",
        field="INFO",
        info_key="DB",
        data_type="flag",
    )

    variant_line = get_variant_line()
    variant_dict = get_variant_dict()

    assert plugin.get_value(variant_line=variant_line) is True
    assert plugin.get_value(variant_dict=variant_dict) is True


def test_rank_score_no_key():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="RankScore=1:12,2:11")
    variant_line = get_variant_line(info="RankScore=1:12,2:11")

    plugin = Plugin(
        name="rank_score",
        data_type="integer",
        field="INFO",
        info_key="RankScore",
        separators=[",", ":"],
        dict_entry=True,
    )
    dict_entry = plugin.get_value(variant_dict=variant_dict)
    line_entry = plugin.get_value(variant_line=variant_line)

    assert dict_entry == 11
    assert line_entry == 11


def test_rank_score_with_key():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="RankScore=1:12,2:11")
    variant_line = get_variant_line(info="RankScore=1:12,2:11")

    plugin = Plugin(
        name="rank_score",
        data_type="integer",
        field="INFO",
        info_key="RankScore",
        separators=[",", ":"],
        dict_entry=True,
    )
    dict_entry = plugin.get_value(variant_dict=variant_dict, dict_key="1")
    line_entry = plugin.get_value(variant_line=variant_line, dict_key="1")

    assert dict_entry == 12
    assert line_entry == 12


def test_rank_models():
    """Test to get the raw chromosome"""
    variant_dict = get_variant_dict(info="GeneticModels=1:AD|AD_dn,2:AR_hom")
    variant_line = get_variant_line(info="GeneticModels=1:AD|AD_dn,2:AR_hom")

    plugin = Plugin(
        name="genetic_models",
        field="INFO",
        data_type="string",
        info_key="GeneticModels",
        separators=[",", ":", "|"],
        dict_entry=True,
        string_rules={"AD_dn": 1, "AD": 2},
    )

    dict_entry = plugin.get_value(variant_dict=variant_dict, dict_key="1")
    line_entry = plugin.get_value(variant_line=variant_line, dict_key="1")

    assert dict_entry == "AD"
    assert line_entry == "AD"
