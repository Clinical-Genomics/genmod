import pytest
from genmod import __version__
from genmod.vcf_tools import (
    HeaderParser,
    add_annotation_header,
    add_compounds_header,
    add_exonic_header,
    add_genetic_models_header,
    add_metadata,
    add_model_score_header,
    add_version_header,
)


def get_header() -> HeaderParser:
    head = HeaderParser()
    head.add_fileformat("VCFv4.2")
    return head


def test_add_metadata_info_adds_an_info_line() -> None:
    # GIVEN a header without INFO lines
    head = get_header()

    # WHEN adding an info entry
    add_metadata(
        head,
        "info",
        "MyKey",
        annotation_number="1",
        entry_type="Integer",
        description="A test key.",
    )

    # THEN the INFO line is in the header
    assert (
        head.info_dict["MyKey"]
        == '##INFO=<ID=MyKey,Number=1,Type=Integer,Description="A test key.">'
    )
    assert head.extra_info["MyKey"]["Type"] == "Integer"


def test_add_metadata_number_defaults_to_dot() -> None:
    # GIVEN a header
    head = get_header()

    # WHEN adding an info entry without a number
    add_metadata(head, "info", "MyKey", entry_type="String", description="A test key.")

    # THEN the number is '.'
    assert "Number=.," in head.info_dict["MyKey"]


def test_add_metadata_unknown_type_raises() -> None:
    # GIVEN a header
    head = get_header()

    # WHEN adding metadata of an unsupported type
    # THEN an error is raised instead of silently doing nothing
    with pytest.raises(ValueError, match="metadata type"):
        add_metadata(head, "format", "MyKey", entry_type="String", description="x")  # type: ignore[arg-type]
    assert head.info_dict == {}


@pytest.mark.parametrize(
    "add_header, key, number, entry_type",
    [
        (add_annotation_header, "Annotation", ".", "String"),
        (add_exonic_header, "Exonic", "0", "Flag"),
        (add_genetic_models_header, "GeneticModels", ".", "String"),
        (add_model_score_header, "ModelScore", "1", "Integer"),
        (add_compounds_header, "Compounds", ".", "String"),
    ],
)
def test_info_header_helpers(add_header, key: str, number: str, entry_type: str) -> None:
    # GIVEN a header
    head = get_header()

    # WHEN adding one of the genmod INFO headers
    add_header(head)

    # THEN the key is added with its number and type
    assert f"ID={key},Number={number},Type={entry_type}," in head.info_dict[key]


def test_add_version_header_tracks_genmod_version() -> None:
    # GIVEN a header
    head = get_header()

    # WHEN adding the version header with the command line
    add_version_header(head, "genmod models trio.vcf")

    # THEN the genmod version and command are recorded in the header lines
    printed = "\n".join(head.print_header())
    assert __version__ in printed
    assert "genmod models trio.vcf" in printed
