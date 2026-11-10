"""String rules: legacy substring matching (default) and opt-in exact matching.

Exact matching (`exact_match=True`) compares whole tokens instead of substrings,
see genmod issue #87.
"""

import pytest
from genmod.plugin import Plugin

CLNSIG_RULES = {
    "pathogenic": 5,
    "likely_pathogenic": 4,
    "uncertain_significance": 3,
    "likely_benign": 2,
    "benign": 1,
}


def get_clnsig_plugin(record_rule="max", separators=("|", "/", ","), **kwargs):
    return Plugin(
        name="CLNSIG",
        field="INFO",
        data_type="string",
        info_key="CLNSIG",
        separators=list(separators),
        record_rule=record_rule,
        string_rules=dict(CLNSIG_RULES),
        **kwargs,
    )


def get_value(clnsig, exact_match=True, **kwargs):
    variant = {"info_dict": {"CLNSIG": clnsig}}
    plugin = get_clnsig_plugin(exact_match=exact_match, **kwargs)
    return plugin.get_value(variant_dict=variant)


# Legacy behaviour: the default must not change


def test_exact_match_is_off_by_default():
    assert get_clnsig_plugin().exact_match is False


@pytest.mark.parametrize(
    "clnsig, expected",
    [
        ("Pathogenic", "pathogenic"),
        ("Likely_benign", "likely_benign"),
        # The substring behaviour reported in #87 is kept for legacy rank models
        ("Likely_pathogenic", "pathogenic"),
        ("Conflicting_classifications_of_pathogenicity", "pathogenic"),
        ("risk_factor", None),
    ],
)
def test_legacy_substring_matching_is_unchanged(clnsig, expected):
    variant = {"info_dict": {"CLNSIG": clnsig}}
    assert get_clnsig_plugin().get_value(variant_dict=variant) == expected
    assert get_value(clnsig, exact_match=False) == expected


def test_legacy_filter_substring_matching_is_unchanged():
    plugin = Plugin(
        name="Filter",
        field="FILTER",
        data_type="string",
        string_rules={"PASS": 2, "NOT_PASS": 1},
        record_rule="max",
    )
    assert plugin.get_value(variant_dict={"FILTER": "NOT_PASS"}) == "PASS"


# Exact matching


@pytest.mark.parametrize(
    "clnsig, expected",
    [
        ("Pathogenic", "pathogenic"),
        ("Likely_pathogenic", "likely_pathogenic"),
        ("Likely_benign", "likely_benign"),
        ("Benign", "benign"),
        ("Uncertain_significance", "uncertain_significance"),
    ],
)
def test_exact_single_classification(clnsig, expected):
    assert get_value(clnsig) == expected


@pytest.mark.parametrize(
    "clnsig",
    [
        "Conflicting_classifications_of_pathogenicity",
        "Conflicting_interpretations_of_pathogenicity",
    ],
)
def test_exact_conflicting_classification_does_not_match_pathogenic(clnsig):
    """'pathogenic' is a substring of 'pathogenicity' but must not match it."""
    assert get_value(clnsig) is None


def test_exact_multiple_submissions_priority_is_kept():
    """Priority policy is unchanged: the best-priority submission wins."""
    assert get_value("Benign|Pathogenic") == "pathogenic"
    assert get_value("Benign|Pathogenic", record_rule="min") == "benign"
    assert get_value("Likely_benign|Likely_pathogenic") == "likely_pathogenic"


def test_exact_likely_pathogenic_alongside_benign_is_not_pathogenic():
    assert get_value("Benign|Likely_pathogenic") == "likely_pathogenic"


def test_exact_slash_separated_classifications():
    assert get_value("Pathogenic/Likely_pathogenic") == "pathogenic"
    assert get_value("Likely_pathogenic/Likely_benign") == "likely_pathogenic"


def test_exact_conflicting_with_real_submission_uses_submission():
    assert get_value("Benign|Conflicting_classifications_of_pathogenicity") == "benign"


def test_exact_no_rule_matches():
    assert get_value("risk_factor") is None


def test_exact_matching_is_case_insensitive():
    assert get_value("LIKELY_PATHOGENIC") == "likely_pathogenic"


def test_exact_filter_field_matches_whole_tokens():
    """PASS;NOT_PASS: each token is compared exactly."""
    plugin = Plugin(
        name="Filter",
        field="FILTER",
        data_type="string",
        string_rules={"PASS": 2, "NOT_PASS": 1},
        record_rule="max",
        exact_match=True,
    )
    assert plugin.get_value(variant_dict={"FILTER": "PASS;NOT_PASS"}) == "PASS"
    assert plugin.get_value(variant_dict={"FILTER": "NOT_PASS"}) == "NOT_PASS"


def test_exact_vep_score_suffix_is_ignored():
    """Older VEP prints e.g. 'probably_damaging(0.95)'; the rule is still 'probably_damaging'."""
    plugin = Plugin(
        name="PolyPhen",
        field="INFO",
        data_type="string",
        info_key="PolyPhen",
        separators=[",", ":", "|"],
        record_rule="max",
        string_rules={"probably_damaging": 3, "possibly_damaging": 2, "benign": 1},
        exact_match=True,
    )
    variant = {"info_dict": {"PolyPhen": "possibly_damaging(0.5)|probably_damaging(0.95)"}}
    assert plugin.get_value(variant_dict=variant) == "probably_damaging"
