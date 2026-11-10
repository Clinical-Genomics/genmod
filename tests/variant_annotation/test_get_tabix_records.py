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


def test_get_tabix_record_unknown_contig(thousand_g_handle):
    assert get_tabix_records(thousand_g_handle, "nonexistent", "100") == []


def test_get_tabix_record_returns_fields(thousand_g_handle):
    record = get_tabix_records(thousand_g_handle, "1", "879537")[0]
    assert record[0] == "1"
    assert record[1] == "879537"
