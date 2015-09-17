from tempfile import NamedTemporaryFile
from multiprocessing import Manager
from collections import OrderedDict

from vcf_parser import VCFParser

from genmod.utils import VariantPrinter

def setup_vcf_file():
    """
    Print some variants to a vcf file and return the filename
    """
    vcf_lines = [
        '##fileformat=VCFv4.1\n',
        '##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">\n',
        '##contig=<ID=1,length=249250621,assembly=b37>\n',
        '##reference=file:///humgen/gsa-hpprojects/GATK/bundle'\
        '/current/b37/human_g1k_v37.fasta\n',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'\
        'father\tmother\tproband\n',
        '1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/1:60\t1/1:60\n',
        '1\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '1\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '1\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '1\t973348\t.\tG\tA\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '3\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n',
        '3\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n',
        '3\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n',
        '3\t973348\t.\tG\tA\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n'
        ]
    vcf_file = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
    vcf_file.writelines(vcf_lines)
    vcf_file.seek(0)
    vcf_file.close()
    
    return vcf_file.name

def test_variant_printer():
    """Test the variant printer"""
    vcf_file = setup_vcf_file()
    variant_queue = Manager().Queue()
    variant_parser = VCFParser(infile=vcf_file)
    header = variant_parser.metadata
    
    outfile = NamedTemporaryFile(mode='w+t', delete=False, suffix='.vcf')
    outfile.close()
    
    
    variant_printer = VariantPrinter(
        task_queue=variant_queue, 
        head=header, 
        mode='chromosome', 
        outfile = outfile.name
    )
    
    variant_printer.start()
    
    batch = OrderedDict()
    
    for variant in variant_parser:
        variant_id = variant['variant_id']
        batch[variant_id] = variant
    
    variant_queue.put(batch)
    variant_queue.put(None)
    
    variant_printer.join()
    
    variants = []
    with open(outfile.name, 'r') as f:
        for line in f:
            variants.append(line.rstrip().split('\t'))
    
    assert variants[0][0] == '1'
    assert variants[0][2] == '11900'