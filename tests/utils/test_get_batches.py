from __future__ import annotations

from queue import Queue

from genmod.utils import get_batches
from genmod.vcf_tools import HeaderParser

HEADER = "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"


def get_variant_line(
    chrom: str = "1",
    pos: str = "1",
    db_id: str = ".",
    ref: str = "A",
    alt: str = "G",
    qual: str = "100",
    filt: str = "PASS",
    info: str = "Annotation=ADK;Exonic",
) -> str:
    """
    Return a variant line
    """
    variant_line = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t".format(
        chrom, pos, db_id, ref, alt, qual, filt, info
    )

    return variant_line


def test_get_batches_one() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
    variants = []
    first_variant = get_variant_line()
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))

    variants.append(first_variant)

    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)
    batch = batch_queue.get()

    assert chromosomes == ["1"]
    assert len(batch) == 1


def test_get_batches_two() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
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


def test_get_batches_two_regions() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
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


def test_get_batches_vep() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
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


def test_get_batches_vep_no_allele() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
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


def test_get_batches_no_regions() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
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


def test_get_batches_new_chromosome() -> None:
    """
    Test to get a batch
    """
    batch_queue: Queue = Queue()
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


def test_get_batches_strips_chr_prefix_from_chromosomes() -> None:
    # GIVEN variants on chromosomes written with and without a 'chr' prefix
    batch_queue: Queue = Queue()
    variants = [get_variant_line(chrom="chr1"), get_variant_line(chrom="chr2")]
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))

    # WHEN creating batches
    chromosomes = get_batches(variants=variants, batch_queue=batch_queue, header=header)

    # THEN the chromosomes are returned without the prefix
    assert chromosomes == ["1", "2"]


def test_get_batches_without_variants_gives_no_chromosomes_and_batches() -> None:
    # GIVEN a header line and a comment line but no variants
    batch_queue: Queue = Queue()
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))

    # WHEN creating batches
    chromosomes = get_batches(variants=["##comment\n"], batch_queue=batch_queue, header=header)

    # THEN no chromosomes are found, and nothing is put in the queue
    assert chromosomes == []
    assert batch_queue.empty()


def test_get_batches_overlapping_annotations_are_kept_together() -> None:
    # GIVEN a gene with variants, and a third variant in another gene
    batch_queue: Queue = Queue()
    variants = [
        get_variant_line(pos="1", info="Annotation=ADK"),
        get_variant_line(pos="2", info="Annotation=ADK,DDD"),
        get_variant_line(pos="3", info="Annotation=DDD"),
        get_variant_line(pos="4", info="Annotation=OTHER"),
    ]
    header = HeaderParser()
    header.parse_header_line("#{0}".format(HEADER))

    # WHEN creating batches
    get_batches(variants=variants, batch_queue=batch_queue, header=header)
    batches = [batch_queue.get(), batch_queue.get()]

    # THEN variants linked through overlapping genes are in one batch, in order
    assert [len(batch) for batch in batches] == [3, 1]
    assert list(batches[0]) == ["1_1_A_G", "1_2_A_G", "1_3_A_G"]
