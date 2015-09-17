from genmod.annotate_regions.parse_annotations import ccds_parser

def get_ccds_lines():
    """
    Create a small ccds file
    """
    ccds_lines = [
        "#chromosome\tnc_accession\tgene\tgene_id\tccds_id\tccds_status\t cds_strand\tcds_from\tcds_to\tcds_locations\tmatch_type",
        "1\tNC_000001.8\tLINC00115\t79854\tCCDS1.1\tWithdrawn\t-\t801942\t802433\t[801942-802433]\tIdentical",
        "2\tNC_000001.10\tSAMD11\t148398\tCCDS2.2\tPublic\t+\t861321\t879532\t[861321-861392, 865534-865715, 866418-866468, 871151-871275, 874419-874508, 874654-874839, 876523-876685, 877515-877630, 877789-877867, 877938-878437, 878632-878756, 879077-879187, 879287-879532]\tIdentical"
    ]
    return ccds_lines

def test_build_ccds_genes():
    """
    Test if the parser builds the ccds gene dict correct
    """
    ccds_lines = get_ccds_lines()
    gene_trees, exon_trees = ccds_parser(ccds_lines)
    assert set(gene_trees['1'].keys()) == set(["LINC00115"])
    assert gene_trees['1']["LINC00115"] == [801942, 802433, "LINC00115"]
    assert exon_trees['1'] == {'801942802433':[801940, 802435, '801942802433']}
