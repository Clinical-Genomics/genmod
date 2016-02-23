import os
import pytest
import tabix
from tempfile import NamedTemporaryFile

@pytest.fixture(scope='function')
def thousand_g_path(request):
    """Return the path to a bgzipped 1000G file"""
    thousand_g = "tests/fixtures/annotate_variant/small_1000G_maxAF.vcf.gz"
    return thousand_g

@pytest.fixture(scope='function')
def thousand_g_handle(request, thousand_g_path):
    """Return a tabix handle with a 1000G file"""
    thousand_g = tabix.open(thousand_g_path)
    return thousand_g

@pytest.fixture(scope='function')
def ccds_lines(request):
    """
    Create a small ccds file
    """
    ccds = [
        "#chromosome\tnc_accession\tgene\tgene_id\tccds_id\tccds_status\t"\
        "cds_strand\tcds_from\tcds_to\tcds_locations\tmatch_type",
        "1\tNC_000001.8\tLINC00115\t79854\tCCDS1.1\tWithdrawn\t-\t801942\t"\
        "802433\t[801942-802433]\tIdentical",
        "2\tNC_000001.10\tSAMD11\t148398\tCCDS2.2\tPublic\t+\t861321\t879532"\
        "\t[861321-861392, 865534-865715, 866418-866468, 871151-871275, 874419"\
        "-874508, 874654-874839, 876523-876685, 877515-877630, 877789-877867,"\
        " 877938-878437, 878632-878756, 879077-879187, 879287-879532]\tIdentical"
    ]
    return ccds


@pytest.fixture(scope='function')
def ref_seq_lines(request):
    ref_seq = [
        "970\tNM_213590\tchr13\t+\t50571142\t50592603\t50586076\t50587300"\
        "\t2\t50571142,50586070,\t50571899,50592603,\t0\tTRIM13\tcmpl\tcmpl\t-1,0,",
        "2166\tNM_001017364\tchr1\t+\t207262583\t207273337\t207262876\t207273274"\
        "\t6\t207262583,207263655,207264988,207269866,207271494,207273133,"\
        "\t207262934,207263826,207265165,207269960,207271609,207273337,\t0\tC4BPB"\
        "\tcmpl\tcmpl\t0,1,1,1,2,0,",
        "1044\tNM_152866\tchr11\t+\t60223281\t60238225\t60229847\t60235941"\
        "\t8\t60223281,60228544,60229657,60230474,60231760,60233393,60234431"\
        ",60235722,\t60223418,60228633,60230006,60230594,60231817,60233630"\
        ",60234533,60238225,\t0\tMS4A1\tcmpl\tcmpl\t-1,-1,0,0,0,0,0,0,",
        "1274\tNM_020466\tchr6\t-\t90341942\t90348474\t90346991\t90348435"\
        "\t3\t90341942,90347460,90348390,\t90347072,90347601,90348474,\t0"\
        "\tLYRM2\tcmpl\tcmpl\t0,0,0,"
    ]
    
    return ref_seq


@pytest.fixture(scope='function')
def gff_lines(request):
    """
    Create a small ccds file
    """
    gff = [
        "##gff-version 3",
        "##sequence-region 121CI_contig1 1 8876",
        "121CI_contig1\tProdigal:2.6\tCDS\t608\t748\t.\t-\t0\tID=121CI_00001;"\
        "inference=ab initio prediction:Prodigal:2.6;locus_tag=121CI_00001;"\
        "product=hypothetical protein",
        "121CI_contig1\tProdigal:2.6\tCDS\t832\t1593\t.\t-\t0\tID=121CI_00002;"\
        "inference=ab initio prediction:Prodigal:2.6,similar to AA sequence"\
        ":Protein_fasta_for_annot.faa:gi|15646195;locus_tag=121CI_00002"\
        ";note=protein of unknown function UPF0174 ;product=hypothetical protein"
    ]
    
    return gff


@pytest.fixture(scope='function')
def ccds_file(request, ccds_lines):
    """docstring for ccds_file"""
    annotations_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.txt')
    annotations_file.write('\n'.join(ccds_lines))
    annotations_file.seek(0)
    annotations_file.close()
    
    def teardown():
        os.remove(annotations_file.name)
    
    request.addfinalizer(teardown)
    
    return annotations_file.name
    
    