"""Multiple VEP consequences are joined by '&', e.g. 'intron_variant&non_coding_transcript_variant'.

With exact matching, '&' has to be one of the plugin separators, see genmod issue #205.
"""

from genmod.plugin import Plugin
from genmod.plugin.plugin import Value

CSQ_FORMAT = ["Allele", "Consequence", "SYMBOL"]
CONSEQUENCE_RULES = {"intron_variant": 1, "non_coding_transcript_variant": 2}
CSQ = "G|intron_variant&non_coding_transcript_variant|GENE1"


def get_consequence_value(separators: list[str], exact_match: bool = True) -> Value:
    plugin = Plugin(
        name="most_severe_consequence",
        field="INFO",
        data_type="string",
        info_key="CSQ",
        csq_key="Consequence",
        separators=separators,
        record_rule="max",
        string_rules=dict(CONSEQUENCE_RULES),
        exact_match=exact_match,
    )
    return plugin.get_value(variant_dict={"info_dict": {"CSQ": CSQ}}, csq_format=CSQ_FORMAT)


def test_ampersand_separator_picks_highest_priority_consequence():
    # GIVEN a rank model with '&' among the separators
    # WHEN scoring a variant with two '&' joined consequences
    value = get_consequence_value([",", ":", "|", "&"])

    # THEN the consequence with the highest priority is chosen, not just the first one
    assert value == "non_coding_transcript_variant"


def test_without_ampersand_separator_exact_match_finds_no_consequence():
    # GIVEN a rank model without '&' among the separators
    # WHEN scoring a variant with two '&' joined consequences using exact matching
    value = get_consequence_value([",", ":", "|"])

    # THEN the joined value is a single token that matches no rule
    assert value is None


def test_without_ampersand_separator_legacy_substring_match_still_works():
    # GIVEN a rank model without '&' among the separators and exact_match = false
    # WHEN scoring a variant with two '&' joined consequences
    value = get_consequence_value([",", ":", "|"], exact_match=False)

    # THEN the substring match finds the highest priority consequence
    assert value == "non_coding_transcript_variant"
