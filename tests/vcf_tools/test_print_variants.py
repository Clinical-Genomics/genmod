from pathlib import Path

import pytest
from genmod.vcf_tools import print_variant, print_variant_dict, print_variant_for_sorting

LINE = "1\t100\t.\tA\tT\t50\tPASS\tMQ=1\n"
HEADER_LINE = ["CHROM", "POS", "ID", "REF", "ALT"]
VARIANT = {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "T"}


def test_print_variant_line_to_file(tmp_path: Path) -> None:
    # GIVEN a variant line and an outfile
    out = tmp_path / "out.vcf"

    # WHEN printing it
    with open(out, "w") as handle:
        print_variant(variant_line=LINE, outfile=handle)

    # THEN the line is written as it came in
    assert out.read_text() == LINE


def test_print_variant_line_to_stdout(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing without an outfile
    print_variant(variant_line=LINE)

    # THEN it is printed to stdout
    assert capsys.readouterr().out == LINE


def test_print_variant_silent(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing silently
    print_variant(variant_line=LINE, silent=True)

    # THEN nothing is printed
    assert capsys.readouterr().out == ""


def test_print_variant_priority_is_prepended(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing with a priority
    print_variant(variant_line=LINE, priority="5")

    # THEN the priority is the first column
    assert capsys.readouterr().out == "5\t" + LINE


def test_print_variant_modified_mode_drops_the_sort_column(
    capsys: pytest.CaptureFixture[str],
) -> None:
    # GIVEN a line where the first column was added for sorting
    sorted_line = "5\t" + LINE

    # WHEN printing in modified mode
    print_variant(variant_line=sorted_line, mode="modified")

    # THEN the sorting column is removed again (and priority is not added)
    assert capsys.readouterr().out == LINE


def test_print_variant_dict_uses_header_columns(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing a variant dict with a header that has a column the variant lacks
    print_variant(variant_dict=VARIANT, header_line=HEADER_LINE)

    # THEN the columns follow the header and missing values are '.'
    assert capsys.readouterr().out == "1\t100\t.\tA\tT\n"


def test_print_variant_dict_without_header_line_raises() -> None:
    # WHEN printing a variant dict without a header line
    # THEN an error is raised
    with pytest.raises(IOError, match="header_line"):
        print_variant(variant_dict=VARIANT)


def test_print_variant_without_line_or_dict_raises() -> None:
    # WHEN there is nothing to print
    # THEN a clear error is raised
    with pytest.raises(ValueError, match="variant_line or variant_dict"):
        print_variant()


def test_print_variant_dict_function(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing with print_variant_dict
    print_variant_dict(VARIANT, HEADER_LINE)

    # THEN the columns follow the header
    assert capsys.readouterr().out == "1\t100\t.\tA\tT\n"


def test_print_variant_for_sorting(tmp_path: Path) -> None:
    # GIVEN a temp file
    out = tmp_path / "tmp.vcf"

    # WHEN printing a variant for sorting
    with open(out, "w") as handle:
        print_variant_for_sorting(LINE, "7", handle)

    # THEN the priority is the first column
    assert out.read_text() == "7\t" + LINE
