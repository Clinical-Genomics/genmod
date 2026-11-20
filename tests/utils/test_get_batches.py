from genmod.utils import SNPEFF, get_batches
from genmod.vcf_tools import HeaderParser

try:
    from Queue import Queue
except ImportError:
    from queue import Queue

HEADER = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def get_variant_line(
    chrom="1",
    pos="1",
    db_id=".",
    ref="A",
    alt="G",
    qual="100",
    filt="PASS",
    info="Annotation=ADK;Exonic",
):
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

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)
    batch = batch_queue.get()

    assert chromosomes == ["1"]
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

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)
    batch = batch_queue.get()

    assert chromosomes == ["1"]
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

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)
    batch_1 = batch_queue.get()
    batch_queue.task_done()

    batch_2 = batch_queue.get()
    batch_queue.task_done()

    assert chromosomes == ["1"]
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
    header.vep_columns = ["Allele", "SYMBOL"]

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)

    batch_1 = batch_queue.get()
    batch_queue.task_done()

    batch_2 = batch_queue.get()
    batch_queue.task_done()

    assert chromosomes == ["1"]
    assert len(batch_1) == 1
    assert len(batch_2) == 1


def test_get_batches_snpeff_shared_gene():
    # GIVEN two SnpEff annotated variants in the same gene and an intergenic variant
    ann_header = (
        '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: '
        "'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID' \">"
    )
    variants = [
        get_variant_line(info="ANN=G|missense_variant|MODERATE|ADK|ENSG01"),
        get_variant_line(pos="2", info="ANN=G|stop_gained|HIGH|ADK|ENSG01"),
        get_variant_line(pos="3", info="ANN=G|intergenic_region|MODIFIER|ADK-DDD|ENSG01-ENSG02"),
    ]
    header = HeaderParser()
    header.parse_meta_data(ann_header)
    header.parse_header_line("#{0}".format(HEADER))
    batch_queue = Queue()

    # WHEN batching with the SnpEff source
    get_batches(variants=variants, batch_queue=batch_queue, header=header, source=SNPEFF)
    batch_1 = batch_queue.get()
    batch_2 = batch_queue.get()

    # THEN the variants in the gene share a batch and the intergenic variant is alone
    assert list(batch_1) == ["1_1_A_G", "1_2_A_G"]
    assert batch_1["1_1_A_G"]["annotation"] == {"ENSG01"}
    assert list(batch_2) == ["1_3_A_G"]
    assert batch_2["1_3_A_G"]["annotation"] == set()


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
    header.vep_columns = ["SYMBOL"]

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)

    batch_1 = batch_queue.get()
    batch_queue.task_done()

    batch_2 = batch_queue.get()
    batch_queue.task_done()

    assert chromosomes == ["1"]
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

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)

    batch_1 = batch_queue.get()
    batch_queue.task_done()

    batch_2 = batch_queue.get()
    batch_queue.task_done()

    assert chromosomes == ["1"]
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

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)

    batch_1 = batch_queue.get()
    batch_queue.task_done()

    batch_2 = batch_queue.get()
    batch_queue.task_done()

    assert chromosomes == ["1", "2"]
    assert len(batch_1) == 1
    assert len(batch_2) == 1
