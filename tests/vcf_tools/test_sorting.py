from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest
from genmod.vcf_tools import sort_variants


def setup_csv_file():
    """
    Print some variants to a vcf file and return the filename
    """
    variant_lines = [
        "1\t1\t11900\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/1:60\t1/1:60\n",
        "3\t3\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n",
        "1\t1\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n",
        "23\tX\t879586\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/1:60\t0/1:60\n",
        "1\t1\t879585\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/1:60\t0/0:60\t0/1:60\n",
        "1\t1\t947378\t.\tA\tT\t100\tPASS\tMQ=1\tGT:GQ\t0/0:60\t0/0:60\t0/1:60\n",
    ]
    csv_file = NamedTemporaryFile(mode="w+t", delete=False, suffix=".vcf")
    csv_file.writelines(variant_lines)
    csv_file.seek(0)
    csv_file.close()

    return csv_file.name


def test_sort_variants():
    """
    Test to sort an unsorted file with variants
    """
    csv_file = setup_csv_file()

    sort_variants(infile=csv_file, mode="chromosome")

    variants = []
    with open(csv_file, "r") as f:
        for line in f:
            variants.append(line.rstrip().split("\t"))

    assert variants[0][1] == "1"
    assert variants[0][2] == "11900"

    assert variants[1][1] == "1"
    assert variants[1][2] == "879585"

    assert variants[2][1] == "1"
    assert variants[2][2] == "879586"

    assert variants[3][1] == "1"
    assert variants[3][2] == "947378"

    assert variants[4][1] == "3"
    assert variants[4][2] == "947378"

    assert variants[5][1] == "X"
    assert variants[5][2] == "879586"


def read_rows(path: str) -> list[list[str]]:
    with open(path, "r") as f:
        return [line.rstrip().split("\t") for line in f]


def test_sort_variants_rank_mode_sorts_on_first_column_descending(tmp_path: Path) -> None:
    # GIVEN lines with a numeric priority in the first column
    infile = tmp_path / "rank.txt"
    infile.write_text("2\ta\n10\tb\n5\tc\n")

    # WHEN sorting on rank
    sort_variants(infile=str(infile), mode="rank")

    # THEN the highest score is first
    assert [row[0] for row in read_rows(str(infile))] == ["10", "5", "2"]


def test_sort_variants_vcf_mode_sorts_on_chromosome_and_position(tmp_path: Path) -> None:
    # GIVEN vcf lines out of order (version sorting keeps 2 before 10)
    infile = tmp_path / "vcf.txt"
    infile.write_text("10\t5\n2\t100\n2\t20\n1\t7\n")

    # WHEN sorting in vcf mode
    sort_variants(infile=str(infile), mode="vcf")

    # THEN they are sorted on chromosome, then on position, both with version sorting
    assert read_rows(str(infile)) == [["1", "7"], ["2", "20"], ["2", "100"], ["10", "5"]]


def test_sort_variants_unknown_mode_raises(tmp_path: Path) -> None:
    # GIVEN a file
    infile = tmp_path / "x.txt"
    infile.write_text("b\na\n")

    # WHEN sorting with an unknown mode
    # THEN an error is raised and the file is not changed
    with pytest.raises(ValueError, match="sort mode"):
        sort_variants(infile=str(infile), mode="alphabet")  # type: ignore[arg-type]
    assert infile.read_text() == "b\na\n"
