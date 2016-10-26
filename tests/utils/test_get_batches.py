from genmod.utils import get_batches
from genmod.vcf_tools import HeaderParser
try:
    from Queue import Queue
except ImportError:
    from queue import Queue

HEADER = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"

def get_variant_line(chrom='1', pos='1', db_id='.', ref='A', alt='G', 
qual='100', filt='PASS', info="Annotation=ADK;Exonic"):
    """
    Return a variant dictionary
    """
    variant_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t".format(
        chrom, pos, db_id, ref, alt, qual, filt, info
    )
    
    return variant_line

def test_get_batches_one():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    first_variant = get_variant_line()
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    
    variants.append(first_variant)
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
    header=header)
    batch = batch_queue.get()
    
    assert chromosomes == ['1']
    assert len(batch) == 1

def test_get_batches_two():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    first_variant = get_variant_line()
    second_variant = get_variant_line(pos="2")
    variants.append(first_variant)
    variants.append(second_variant)
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
    header=header)
    batch = batch_queue.get()
    
    assert chromosomes == ['1']
    assert len(batch) == 2


def test_get_batches_two_regions():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    first_variant = get_variant_line()
    second_variant = get_variant_line(pos="2", info="Annotation=DDD;Exonic")
    variants.append(first_variant)
    variants.append(second_variant)
    
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
    header=header)
    batch_1 = batch_queue.get()
    batch_queue.task_done()
    
    batch_2 = batch_queue.get()
    batch_queue.task_done()
    
    assert chromosomes == ['1']
    assert len(batch_1) == 1
    assert len(batch_2) == 1

def test_get_batches_vep():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    
    first_variant = get_variant_line(info="MQ;CSQ=G|ADK")
    
    second_variant = get_variant_line(pos="2", info="MQ;CSQ=G|ADK")
    
    variants.append(first_variant)
    variants.append(second_variant)

    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    header.vep_columns = ['Allele', 'SYMBOL']
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
                header=header)
    
    batch_1 = batch_queue.get()
    batch_queue.task_done()
    
    batch_2 = batch_queue.get()
    batch_queue.task_done()
    
    assert chromosomes == ['1']
    assert len(batch_1) == 1
    assert len(batch_2) == 1

def test_get_batches_vep_no_allele():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    
    first_variant = get_variant_line(info="MQ;CSQ=ADK")
    
    second_variant = get_variant_line(pos="2", info="MQ;CSQ=ADK")
    
    variants.append(first_variant)
    variants.append(second_variant)

    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    header.vep_columns = ['SYMBOL']
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
                header=header)
    
    batch_1 = batch_queue.get()
    batch_queue.task_done()
    
    batch_2 = batch_queue.get()
    batch_queue.task_done()
    
    assert chromosomes == ['1']
    assert len(batch_1) == 1
    assert len(batch_2) == 1

def test_get_batches_no_regions():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    
    first_variant = get_variant_line(info="MQ")
    
    second_variant = get_variant_line(pos="2")
    
    variants.append(first_variant)
    variants.append(second_variant)

    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
    header=header)
    
    batch_1 = batch_queue.get()
    batch_queue.task_done()
    
    batch_2 = batch_queue.get()
    batch_queue.task_done()
    
    assert chromosomes == ['1']
    assert len(batch_1) == 1
    assert len(batch_2) == 1

def test_get_batches_new_chromosome():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    
    first_variant = get_variant_line()
    second_variant = get_variant_line(chrom="2")
    
    variants.append(first_variant)
    variants.append(second_variant)
    
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))
    
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, 
    header=header)
    
    batch_1 = batch_queue.get()
    batch_queue.task_done()
    
    batch_2 = batch_queue.get()
    batch_queue.task_done()
    
    assert chromosomes == ['1', '2']
    assert len(batch_1) == 1
    assert len(batch_2) == 1
    
