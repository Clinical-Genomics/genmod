from genmod.utils import get_batches
from Queue import Queue

def get_variant(chrom='1', pos='1', ref='A', alt='G', annotation=["ADK"]):
    """
    Return a variant dictionary
    """
    variant_id = '_'.join([chrom, pos, ref, alt])
    variant = {
        "CHROM":chrom,
        "POS":pos,
        "INFO":"Annotation={0}".format(annotation),
        'info_dict':{
            "Annotation":[region for region in annotation]
        },
        "variant_id": variant_id
    }
    return variant

def test_get_batches_one():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    first_variant = get_variant()
    
    variants.append(first_variant)
    
    chromosomes = get_batches(variants, batch_queue)
    batch = batch_queue.get()
    
    assert chromosomes == ['1']
    assert len(batch) == 1

def test_get_batches_two():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    first_variant = get_variant()
    second_variant = get_variant(pos="2")
    variants.append(first_variant)
    variants.append(second_variant)
    chromosomes = get_batches(variants, batch_queue)
    batch = batch_queue.get()
    
    assert chromosomes == ['1']
    assert len(batch) == 2


def test_get_batches_two_regions():
    """
    Test to get a batch
    """
    batch_queue = Queue()
    variants = []
    first_variant = get_variant()
    second_variant = get_variant(pos="2", annotation=["DDD"])
    variants.append(first_variant)
    variants.append(second_variant)
    chromosomes = get_batches(variants, batch_queue)
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
    
    first_variant = get_variant()
    first_variant['info_dict'].pop('Annotation')
    second_variant = get_variant(pos="2")
    
    variants.append(first_variant)
    variants.append(second_variant)
    chromosomes = get_batches(variants, batch_queue)
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
    
    first_variant = get_variant()
    second_variant = get_variant(chrom="2")
    
    variants.append(first_variant)
    variants.append(second_variant)
    chromosomes = get_batches(variants, batch_queue)
    batch_1 = batch_queue.get()
    batch_queue.task_done()
    
    batch_2 = batch_queue.get()
    batch_queue.task_done()
    
    assert chromosomes == ['1', '2']
    assert len(batch_1) == 1
    assert len(batch_2) == 1
    
