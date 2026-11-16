"""String rules: exact matching (default) and legacy substring matching (`exact_match=False`).

Exact matching compares whole tokens instead of substrings, see genmod issue #87.
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


def test_exact_match_is_on_by_default():
    # GIVEN a plugin created without the exact_match argument
    plugin = get_clnsig_plugin()

    # THEN exact matching is on
    assert plugin.exact_match is True


def test_default_matching_is_exact():
    # GIVEN a plugin without the exact_match argument and a Likely_pathogenic variant
    variant = {"info_dict": {"CLNSIG": "Likely_pathogenic"}}

    # WHEN getting the plugin value
    value = get_clnsig_plugin().get_value(variant_dict=variant)

    # THEN 'pathogenic' does not match 'Likely_pathogenic'
    assert value == "likely_pathogenic"


# Legacy behaviour: exact_match=False keeps the substring matching


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
    # GIVEN a variant with a CLNSIG value and a plugin with exact_match=False
    # WHEN getting the plugin value
    value = get_value(clnsig, exact_match=False)

    # THEN rules are matched as substrings, as before exact_match existed
    assert value == expected


def test_legacy_filter_substring_matching_is_unchanged():
    # GIVEN a FILTER plugin with exact_match=False, where PASS is a substring of NOT_PASS
    plugin = Plugin(
        name="Filter",
        field="FILTER",
        data_type="string",
        string_rules={"PASS": 2, "NOT_PASS": 1},
        record_rule="max",
        exact_match=False,
    )

    # WHEN getting the value of a NOT_PASS variant
    value = plugin.get_value(variant_dict={"FILTER": "NOT_PASS"})

    # THEN the substring still matches the PASS rule (legacy behaviour)
    assert value == "PASS"


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
    # GIVEN a single classification and a plugin with exact_match
    # WHEN getting the plugin value
    value = get_value(clnsig)

    # THEN the rule of exactly that classification is returned
    assert value == expected


@pytest.mark.parametrize(
    "clnsig",
    [
        "Conflicting_classifications_of_pathogenicity",
        "Conflicting_interpretations_of_pathogenicity",
    ],
)
def test_exact_conflicting_classification_does_not_match_pathogenic(clnsig):
    """'pathogenic' is a substring of 'pathogenicity' but must not match it."""
    # GIVEN a conflicting classifications value and a plugin with exact_match
    # WHEN getting the plugin value
    value = get_value(clnsig)

    # THEN no rule matches, so the variant is not scored as pathogenic
    assert value is None


def test_exact_multiple_submissions_priority_is_kept():
    """Priority policy is unchanged: the best-priority submission wins."""
    # GIVEN values with several submissions and a plugin with exact_match
    # WHEN getting the value with record_rule max and min
    # THEN the highest (max) or lowest (min) priority rule wins
    assert get_value("Benign|Pathogenic") == "pathogenic"
    assert get_value("Benign|Pathogenic", record_rule="min") == "benign"
    assert get_value("Likely_benign|Likely_pathogenic") == "likely_pathogenic"


def test_exact_likely_pathogenic_alongside_benign_is_not_pathogenic():
    # GIVEN a Benign and a Likely_pathogenic submission and a plugin with exact_match
    # WHEN getting the plugin value
    value = get_value("Benign|Likely_pathogenic")

    # THEN it is Likely_pathogenic, not Pathogenic
    assert value == "likely_pathogenic"


def test_exact_slash_separated_classifications():
    # GIVEN ClinVar values combined with '/' and '/' as a separator
    # WHEN getting the plugin value
    # THEN each part is matched on its own and the best priority wins
    assert get_value("Pathogenic/Likely_pathogenic") == "pathogenic"
    assert get_value("Likely_pathogenic/Likely_benign") == "likely_pathogenic"


def test_exact_conflicting_with_real_submission_uses_submission():
    # GIVEN a real Benign submission next to a conflicting classifications value
    # WHEN getting the plugin value
    value = get_value("Benign|Conflicting_classifications_of_pathogenicity")

    # THEN only the real submission matches a rule
    assert value == "benign"


def test_exact_no_rule_matches():
    # GIVEN a value that none of the rules know and a plugin with exact_match
    # WHEN getting the plugin value
    # THEN there is no value, so the not_reported score applies
    assert get_value("risk_factor") is None


def test_exact_matching_is_case_insensitive():
    # GIVEN an upper case value and a plugin with exact_match
    # WHEN getting the plugin value
    # THEN it matches the lower case rule
    assert get_value("LIKELY_PATHOGENIC") == "likely_pathogenic"


def test_exact_filter_field_matches_whole_tokens():
    """PASS;NOT_PASS: each token is compared exactly."""
    # GIVEN a FILTER plugin with exact_match, where PASS is a substring of NOT_PASS
    plugin = Plugin(
        name="Filter",
        field="FILTER",
        data_type="string",
        string_rules={"PASS": 2, "NOT_PASS": 1},
        record_rule="max",
        exact_match=True,
    )

    # WHEN getting the value of variants with PASS;NOT_PASS and with NOT_PASS
    # THEN NOT_PASS no longer matches the PASS rule
    assert plugin.get_value(variant_dict={"FILTER": "PASS;NOT_PASS"}) == "PASS"
    assert plugin.get_value(variant_dict={"FILTER": "NOT_PASS"}) == "NOT_PASS"


def test_exact_vep_score_suffix_is_ignored():
    """Older VEP prints e.g. 'probably_damaging(0.95)'; the rule is still 'probably_damaging'."""
    # GIVEN a PolyPhen plugin with exact_match
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
    # AND a variant with VEP scores after the predictions
    variant = {"info_dict": {"PolyPhen": "possibly_damaging(0.5)|probably_damaging(0.95)"}}

    # WHEN getting the plugin value
    value = plugin.get_value(variant_dict=variant)

    # THEN the scores are ignored and the best priority prediction wins
    assert value == "probably_damaging"


# CSQ column (real VEP annotation of ClinVar variation 198801, trimmed to relevant columns)

CSQ_FORMAT = [
    "Allele",
    "Consequence",
    "SYMBOL",
    "Feature",
    "SIFT",
    "PolyPhen",
    "CLNVID",
    "CLNSIG",
    "CLNVID_2",
    "CLNREVSTAT",
]

CONFLICTING_CSQ = ",".join(
    [
        "G|missense_variant|VCAN|ENST00000265077|tolerated|benign|198801|"
        "Conflicting_classifications_of_pathogenicity|198801|"
        "criteria_provided&_conflicting_classifications",
        "G|intron_variant|VCAN|ENST00000342785|||198801|"
        "Conflicting_classifications_of_pathogenicity|198801|"
        "criteria_provided&_conflicting_classifications",
    ]
)


def get_csq_plugin(**kwargs):
    return Plugin(
        name="CLNSIG",
        field="INFO",
        data_type="string",
        info_key="CSQ",
        csq_key="CLNSIG",
        separators=["&"],
        record_rule="max",
        string_rules=dict(CLNSIG_RULES),
        **kwargs,
    )


def get_csq_value(csq, **kwargs):
    variant = {"info_dict": {"CSQ": csq}}
    return get_csq_plugin(**kwargs).get_value(variant_dict=variant, csq_format=CSQ_FORMAT)


def test_exact_csq_conflicting_classification_does_not_match_pathogenic():
    # GIVEN the CSQ of ClinVar variation 198801 (conflicting classifications)
    # WHEN getting the CLNSIG value with exact_match
    value = get_csq_value(CONFLICTING_CSQ, exact_match=True)

    # THEN 'pathogenic' does not match 'Conflicting_..._pathogenicity'
    assert value is None


def test_exact_csq_only_uses_the_csq_key_column():
    """'benign' is in the PolyPhen column but only the CLNSIG column is searched."""
    # GIVEN a CSQ with 'benign' in the PolyPhen column but not in the CLNSIG column
    # WHEN getting the CLNSIG value with exact_match
    value = get_csq_value(CONFLICTING_CSQ, exact_match=True)

    # THEN the PolyPhen value is not matched
    assert value != "benign"


def test_exact_csq_classification_across_transcripts_and_separators():
    csq = (
        "G|missense_variant|VCAN|ENST1|tolerated|benign|1|Benign&Likely_benign|1|x,"
        "G|intron_variant|VCAN|ENST2|||1|Likely_pathogenic|1|x"
    )
    # GIVEN a CSQ with two transcripts and an '&' separated CLNSIG value
    # WHEN getting the CLNSIG value with exact_match
    value = get_csq_value(csq, exact_match=True)

    # THEN the values of all transcripts are considered and the best priority wins
    assert value == "likely_pathogenic"


def test_legacy_csq_searches_the_whole_raw_value():
    """Legacy behaviour, kept as is: substring of the whole CSQ string."""
    # GIVEN the CSQ of ClinVar variation 198801 and a plugin with exact_match=False
    # WHEN getting the CLNSIG value
    value = get_csq_value(CONFLICTING_CSQ, exact_match=False)

    # THEN the whole raw CSQ is searched and 'pathogenic' matches 'pathogenicity'
    assert value == "pathogenic"
