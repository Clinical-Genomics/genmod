from tempfile import NamedTemporaryFile

from genmod.annotate_regions import parse_annotations



def create_ccds_file():
    """
    Create a small ccds file
    """
    ccds_lines = [
        "#chromosome\tnc_accession\tgene\tgene_id\tccds_id\tccds_status\t cds_strand\tcds_from\tcds_to\tcds_locations\tmatch_type",
        "1\tNC_000001.8\tLINC00115\t79854\tCCDS1.1\tWithdrawn\t-\t801942\t802433\t[801942-802433]\tIdentical",
        "2\tNC_000001.10\tSAMD11\t148398\tCCDS2.2\tPublic\t+\t861321\t879532\t[861321-861392, 865534-865715, 866418-866468, 871151-871275, 874419-874508, 874654-874839, 876523-876685, 877515-877630, 877789-877867, 877938-878437, 878632-878756, 879077-879187, 879287-879532]\tIdentical"
    ]
    
    annotations_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.txt')
    annotations_file.writelines(ccds_lines)
    annotations_file.seek(0)
    annotations_file.close()
    
    return annotations_file.name

def test_parse_annotations():
    """
    Test parse_annotations
    """
    ccds_file = create_ccds_file()
    gene_trees, exon_trees = parse_annotations(infile=ccds_file, annotation_type="ccds")
    # assert gene_trees == {'hej':0}
    # assert set(gene_trees['1'].keys()) == set(["LINC00115"])
    # assert gene_trees['1']["LINC00115"] == [801942, 802433, "LINC00115"]
    # assert exon_trees['1'] == {'801942802433':[801940, 802435, '801942802433']}
