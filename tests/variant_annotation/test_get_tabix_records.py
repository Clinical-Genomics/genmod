from pathlib import Path

import pytest
from genmod.annotate_variants.read_tabix_files import get_contig_name, get_tabix_records


def test_get_tabix_record(thousand_g_handle):
    chrom = "1"
    start = "879537"
    i = None
    for i, row in enumerate(get_tabix_records(thousand_g_handle, chrom, start)):
        print(row)
    # Should find one row
    assert i == 0


def test_get_tabix_record_chr(thousand_g_chr_handle):
    chrom = "1"
    start = "879537"
    i = None
    for i, row in enumerate(get_tabix_records(thousand_g_chr_handle, chrom, start)):
        print(row)
    # Should find one row
    assert i == 0


def test_get_non_existing_tabix_record(thousand_g_handle):
    chrom = "1"
    start = "10"
    i = None
    for i, row in enumerate(get_tabix_records(thousand_g_handle, chrom, start)):
        print(row)
    # Should find one row
    assert i is None


def test_get_tabix_record_unknown_contig(thousand_g_path, thousand_g_handle):
    # GIVEN an indexed file that has no contig "nonexistent" (with or without "chr")
    assert Path(thousand_g_path).exists()

    # WHEN fetching records for that contig
    records = get_tabix_records(thousand_g_handle, "nonexistent", "100")

    # THEN no records are returned
    assert records == []


def test_get_tabix_record_returns_fields(thousand_g_path, thousand_g_handle):
    # GIVEN an indexed file with a record at 1:879537
    assert Path(thousand_g_path).exists()

    # WHEN fetching that position
    record = get_tabix_records(thousand_g_handle, "1", "879537")[0]

    # THEN the record can be indexed by column
    assert record[0] == "1"
    assert record[1] == "879537"


def test_get_tabix_record_chr_input_unprefixed_file(thousand_g_handle):
    # GIVEN a file with contigs without prefix ("1")
    assert "1" in thousand_g_handle.contigs

    # WHEN fetching with a prefixed chromosome ("chr1")
    records = get_tabix_records(thousand_g_handle, "chr1", "879537")

    # THEN the 1:879537 record is found
    assert len(records) == 1
    assert records[0][1] == "879537"


def test_get_tabix_record_chr_input_prefixed_file(thousand_g_chr_handle):
    # GIVEN a file with prefixed contigs ("chr1") and no "chrchr1"
    assert "chr1" in thousand_g_chr_handle.contigs
    assert "chrchr1" not in thousand_g_chr_handle.contigs

    # WHEN fetching with a prefixed chromosome ("chr1")
    records = get_tabix_records(thousand_g_chr_handle, "chr1", "879537")

    # THEN the 1:879537 record is found
    assert len(records) == 1
    assert records[0][1] == "879537"


def test_get_tabix_record_unprefixed_input_prefixed_file(thousand_g_chr_handle):
    # GIVEN a file with prefixed contigs ("chr1")
    assert "chr1" in thousand_g_chr_handle.contigs

    # WHEN fetching with an unprefixed chromosome ("1")
    records = get_tabix_records(thousand_g_chr_handle, "1", "879537")

    # THEN the 1:879537 record is found
    assert len(records) == 1
    assert records[0][1] == "879537"


@pytest.mark.parametrize(
    "chrom, expected",
    [("1", "chr1"), ("chr1", "chr1"), ("X", "chrX"), ("nonexistent", None)],
)
def test_get_contig_name_prefixed_file(thousand_g_chr_handle, chrom, expected):
    # GIVEN a file with prefixed contigs
    # WHEN looking up a chromosome, with or without prefix
    # THEN the contig as named in the file is returned, or None if it is not there
    assert get_contig_name(thousand_g_chr_handle, chrom) == expected


@pytest.mark.parametrize("chrom, expected", [("1", "1"), ("chr1", "1"), ("chrchr1", None)])
def test_get_contig_name_unprefixed_file(thousand_g_handle, chrom, expected):
    # GIVEN a file with unprefixed contigs
    # WHEN looking up a chromosome, with or without prefix
    # THEN the contig as named in the file is returned, or None if it is not there
    assert get_contig_name(thousand_g_handle, chrom) == expected
