from pathlib import Path

import pytest
from genmod.vcf_tools import HeaderParser, print_headers


def get_header() -> HeaderParser:
    head = HeaderParser()
    head.add_fileformat("VCFv4.2")
    head.parse_header_line("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")
    return head


def test_print_headers_to_file(tmp_path: Path) -> None:
    # GIVEN a header and an open file
    out = tmp_path / "out.vcf"

    # WHEN printing the header to the file
    with open(out, "w") as handle:
        print_headers(get_header(), outfile=handle)

    # THEN the file has the header lines
    lines = out.read_text().splitlines()
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[-1].startswith("#CHROM\tPOS")


def test_print_headers_to_stdout(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing without an outfile
    print_headers(get_header())

    # THEN the header is printed to stdout
    assert capsys.readouterr().out.startswith("##fileformat=VCFv4.2\n#CHROM")


def test_print_headers_silent_prints_nothing(capsys: pytest.CaptureFixture[str]) -> None:
    # WHEN printing silently without an outfile
    print_headers(get_header(), silent=True)

    # THEN nothing is printed
    assert capsys.readouterr().out == ""


def test_print_headers_silent_still_writes_to_file(tmp_path: Path) -> None:
    # GIVEN an outfile and silent=True
    out = tmp_path / "out.vcf"

    # WHEN printing
    with open(out, "w") as handle:
        print_headers(get_header(), outfile=handle, silent=True)

    # THEN silent only affects stdout, the file still gets the header
    assert out.read_text().startswith("##fileformat=VCFv4.2")
