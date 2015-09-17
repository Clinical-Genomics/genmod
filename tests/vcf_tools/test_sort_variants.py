from tempfile import NamedTemporaryFile

from genmod.vcf_tools import sort_variants

def setup_csv_file():
    """
    Print some variants to a vcf file and return the filename
    """
    variant_lines = [
        '1\t1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/1:60\t1/1:60\n',
        '3\t3\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '1\t1\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '23\tX\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '1\t1\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '1\t1\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        ]
    csv_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
    csv_file.writelines(variant_lines)
    csv_file.seek(0)
    csv_file.close()
    
    return csv_file.name

def test_sort_variants():
    """
    Test to sort an unsorted file with variants
    """
    csv_file = setup_csv_file()
    
    sort_variants(infile=csv_file, mode='chromosome')
    
    variants = []
    with open(csv_file, 'r') as f:
        for line in f:
            variants.append(line.rstrip().split('\t'))

    assert variants[0][1] == '1'    
    assert variants[0][2] == '11900'
    
    assert variants[1][1] == '1'    
    assert variants[1][2] == '879585'
    
    assert variants[2][1] == '1'    
    assert variants[2][2] == '879586'

    assert variants[3][1] == '1'    
    assert variants[3][2] == '947378'
    
    assert variants[4][1] == '3'    
    assert variants[4][2] == '947378'
    
    assert variants[5][1] == 'X'    
    assert variants[5][2] == '879586'

