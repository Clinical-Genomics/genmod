from pathlib import Path

from genmod.annotate_variants.read_tabix_files import get_tabix_records


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
