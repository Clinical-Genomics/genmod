import gzip
import io
from pathlib import Path

import pytest
from genmod.commands.utils import get_file_handle

LINES = ["##fileformat=VCFv4.2\n", "#CHROM\tPOS\n", "1\t100\n"]


def test_get_file_handle_plain_file(tmp_path: Path) -> None:
    # GIVEN a plain text file
    path = tmp_path / "variants.vcf"
    path.write_text("".join(LINES))

    # WHEN opening it
    with get_file_handle(str(path)) as handle:
        # THEN the lines can be iterated as text
        assert list(handle) == LINES


def test_get_file_handle_gzipped_file(tmp_path: Path) -> None:
    # GIVEN a gzipped file
    path = tmp_path / "variants.vcf.gz"
    with gzip.open(path, "wt", encoding="utf-8") as out:
        out.writelines(LINES)

    # WHEN opening it
    with get_file_handle(str(path)) as handle:
        # THEN the lines are decoded to text
        assert list(handle) == LINES


def test_get_file_handle_gzipped_file_with_invalid_utf8_is_replaced(tmp_path: Path) -> None:
    # GIVEN a gzipped file with a byte that is not valid utf-8
    path = tmp_path / "bad.vcf.gz"
    with gzip.open(path, "wb") as out:
        out.write(b"1\t100\tcaf\xe9\n")

    # WHEN reading it
    with get_file_handle(str(path)) as handle:
        lines = list(handle)

    # THEN the byte is replaced instead of raising an error
    assert lines == ["1\t100\tcaf�\n"]


def test_get_file_handle_utf8_content(tmp_path: Path) -> None:
    # GIVEN a file with non-ascii characters, plain and gzipped
    text = "#comment Åström\n"
    plain = tmp_path / "u.vcf"
    plain.write_text(text, encoding="utf-8")
    zipped = tmp_path / "u.vcf.gz"
    with gzip.open(zipped, "wt", encoding="utf-8") as out:
        out.write(text)

    # THEN both are read correctly
    for path in (plain, zipped):
        with get_file_handle(str(path)) as handle:
            assert handle.read() == text


def test_get_file_handle_dash_is_stdin(monkeypatch: pytest.MonkeyPatch) -> None:
    # GIVEN something on stdin
    fake_stdin = io.StringIO("".join(LINES))
    monkeypatch.setattr("sys.stdin", fake_stdin)

    # WHEN opening '-'
    handle = get_file_handle("-")

    # THEN stdin itself is returned
    assert handle is fake_stdin
    assert list(handle) == LINES


def test_get_file_handle_missing_file_raises(tmp_path: Path) -> None:
    # WHEN opening a file that does not exist
    # THEN an error is raised
    with pytest.raises(FileNotFoundError):
        get_file_handle(str(tmp_path / "missing.vcf"))
    with pytest.raises(FileNotFoundError):
        get_file_handle(str(tmp_path / "missing.vcf.gz"))
