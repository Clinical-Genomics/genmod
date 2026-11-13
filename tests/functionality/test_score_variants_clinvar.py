"""End-to-end tests of scoring ClinVar classifications with string rules (genmod issue #87).

The VCF has a real VEP annotated variant, ClinVar variation 198801, classified as
"Conflicting_classifications_of_pathogenicity", together with control variants. The ClinVar
classification is read from a CSQ column (``csq_key = CLNSIG``).

By default whole tokens in the CLNSIG column are compared. With ``exact_match = false`` rules are
matched as substrings of the raw CSQ value (legacy behaviour), so the conflicting variant is
scored as pathogenic.
"""

import operator
from pathlib import Path
from typing import Dict, Tuple

import pytest
from click.testing import CliRunner
from genmod.commands import score_command
from genmod.plugin import Plugin
from genmod.score_variants import config_parser
from test_utils import generate_variants_from_file

CLINVAR_VCF = "tests/fixtures/score_variants/clinvar_csq.vcf"
EXACT_CONFIG = "tests/fixtures/score_variants/clinvar_string_exact.ini"
LEGACY_CONFIG = "tests/fixtures/score_variants/clinvar_string_legacy.ini"
EXAMPLE_CONFIG = "tests/fixtures/score_variants/genmod_example.ini"
ANNOTATED_VCF = "tests/fixtures/test_vcf_annotated.vcf"

# Positions in clinvar_csq.vcf, with the CLNSIG column of the variant
CONFLICTING = 83539905  # Conflicting_classifications_of_pathogenicity (ClinVar 198801)
PATHOGENIC = 83540001  # Pathogenic
LIKELY_PATHOGENIC = 83540002  # Likely_pathogenic
LIKELY_BENIGN = 83540003  # Likely_benign
BENIGN = 83540004  # Benign
PATHOGENIC_MIX = 83540005  # Pathogenic&Likely_pathogenic
BENIGN_MIX = 83540006  # Benign&Likely_benign
NO_CLINVAR = 83540007  # No ClinVar annotation


def score_vcf(vcf: str, config: str, tmp_path: Path) -> Dict[str, float]:
    """Run `genmod score` and return the rank score of each variant keyed by CHROM:POS."""
    result = CliRunner().invoke(score_command, [vcf, "-c", config])
    assert result.exit_code == 0, result.output

    scored_vcf = tmp_path / "scored.vcf"
    scored_vcf.write_bytes(result.stdout_bytes)

    scores = {}
    for variant in generate_variants_from_file(str(scored_vcf)):
        family_score = variant["info_dict"]["RankScore"]  # family_id:score
        scores[f"{variant['CHROM']}:{variant['POS']}"] = float(family_score.split(":")[1])
    return scores


def clinvar_scores(config: str, tmp_path: Path) -> Dict[int, float]:
    """Score the ClinVar VCF and return the rank score keyed by position."""
    scores = score_vcf(CLINVAR_VCF, config, tmp_path)
    return {int(key.split(":")[1]): score for key, score in scores.items()}


def test_exact_match_conflicting_classification_is_not_scored_as_pathogenic(tmp_path: Path):
    """The reported bug: 'pathogenic' must not match 'Conflicting_..._pathogenicity'."""
    # GIVEN a rank model with the default matching and VCAN variant 198801, which has
    # conflicting classifications and benign SIFT/PolyPhen predictions
    # WHEN scoring the variants
    scores = clinvar_scores(EXACT_CONFIG, tmp_path)

    # THEN the conflicting variant is scored as a variant without ClinVar annotation
    assert scores[CONFLICTING] == scores[NO_CLINVAR] == 0.0
    # THEN it scores clearly below a pathogenic variant
    assert scores[CONFLICTING] < scores[PATHOGENIC]


def test_exact_match_scores_each_classification_by_its_own_rule(tmp_path):
    # GIVEN a rank model with the default matching
    # WHEN scoring the control variants
    scores = clinvar_scores(EXACT_CONFIG, tmp_path)

    # THEN each classification gets the score of its own rule
    assert scores[PATHOGENIC] == 5.0
    # 'pathogenic' must not match 'Likely_pathogenic'
    assert scores[LIKELY_PATHOGENIC] == 2.0
    # benign and likely benign variants are not scored high
    assert scores[LIKELY_BENIGN] == -1.0
    assert scores[BENIGN] == -2.0
    # The highest priority classification wins if there are several (unchanged policy)
    assert scores[PATHOGENIC_MIX] == 5.0
    assert scores[BENIGN_MIX] == -1.0


def test_exact_match_false_keeps_the_legacy_scores(tmp_path: Path):
    """With exact_match = false the scores are what they always were, including #87."""
    # GIVEN a rank model with exact_match = false
    # WHEN scoring the variants
    scores = clinvar_scores(LEGACY_CONFIG, tmp_path)

    # THEN rules are matched as substrings of the raw CSQ value
    assert scores == {
        CONFLICTING: 5.0,  # 'Pathogenic' in '..._of_pathogenicity'
        PATHOGENIC: 5.0,
        LIKELY_PATHOGENIC: 5.0,  # 'Pathogenic' in 'Likely_pathogenic'
        LIKELY_BENIGN: -1.0,
        BENIGN: -2.0,
        PATHOGENIC_MIX: 5.0,
        BENIGN_MIX: -1.0,
        NO_CLINVAR: -2.0,  # 'Benign' in the PolyPhen value, the whole CSQ is searched
    }


def test_exact_match_true_is_the_same_as_the_default(tmp_path: Path):
    # GIVEN the default rank model with an explicit exact_match = true
    config = tmp_path / "explicit_true.ini"
    config.write_text(
        Path(EXACT_CONFIG)
        .read_text()
        .replace("  csq_key = CLNSIG\n", "  csq_key = CLNSIG\n  exact_match = true\n")
    )
    assert "exact_match = true" in config.read_text()

    # WHEN scoring the variants
    explicit = clinvar_scores(str(config), tmp_path)

    # THEN the scores are those of the default
    assert explicit == clinvar_scores(EXACT_CONFIG, tmp_path)


def test_exact_match_only_differs_for_the_string_plugin(tmp_path):
    """The two rank models are identical except for exact_match = false in the legacy one."""
    # GIVEN the default and the legacy ClinVar rank models
    exact = Path(EXACT_CONFIG).read_text()
    legacy = Path(LEGACY_CONFIG).read_text()

    # WHEN comparing them line by line
    differing_lines = set(exact.splitlines()) ^ set(legacy.splitlines())
    # THEN only the exact_match setting and the model name differ, so score differences
    # between them are caused by exact_match alone
    assert {line.strip() for line in differing_lines if not line.startswith("#")} == {
        "exact_match = false",
        "name = clinvar string exact",
        "name = clinvar string legacy",
    }


class FrozenLegacyPlugin(Plugin):
    """A copy of the string matching genmod had before `exact_match` existed.

    String rules are matched as substrings of the raw entry, ignoring `exact_match`. It is kept
    here, as it was, to check that `exact_match = false` still scores as before.
    """

    def get_value(self, **kwargs):
        if self.data_type != "string" or not self.record_rule:
            return super().get_value(**kwargs)

        raw_entry = self.get_raw_entry(
            variant_line=kwargs.get("variant_line"),
            variant_dict=kwargs.get("variant_dict"),
            vcf_header=kwargs.get("vcf_header"),
            individual_id=kwargs.get("individual_id"),
            dict_key=kwargs.get("dict_key"),
        )
        if not raw_entry:
            return None

        sorted_rules = sorted(
            self.string_rules.items(), key=operator.itemgetter(1), reverse=self.record_rule == "max"
        )
        for rule, _ in sorted_rules:
            if rule.lower() in raw_entry.lower():
                return rule
        return None


def score_with_rank_results(vcf: str, config: str, tmp_path: Path) -> Dict[str, Tuple[str, str]]:
    """Run `genmod score -r` and return (rank score, score of each category) per CHROM:POS."""
    result = CliRunner().invoke(score_command, [vcf, "-c", config, "-r"])
    assert result.exit_code == 0, result.output

    scored_vcf = tmp_path / "scored.vcf"
    scored_vcf.write_bytes(result.stdout_bytes)

    scores = {}
    for variant in generate_variants_from_file(str(scored_vcf)):
        info = variant["info_dict"]
        scores[f"{variant['CHROM']}:{variant['POS']}"] = (info["RankScore"], info["RankResult"])
    return scores


def opt_out_of_exact_match(config: str, tmp_path: Path) -> str:
    """Return a copy of a rank model with exact_match = false on all string plugins."""
    text = Path(config).read_text()
    if "exact_match" not in text:
        text = text.replace("data_type = string\n", "data_type = string\n  exact_match = false\n")
    opted_out = tmp_path / "exact_match_false.ini"
    opted_out.write_text(text)
    return str(opted_out)


@pytest.mark.parametrize(
    "vcf, config",
    [(ANNOTATED_VCF, EXAMPLE_CONFIG), (CLINVAR_VCF, LEGACY_CONFIG)],
    ids=["example_rank_model", "clinvar_legacy_rank_model"],
)
def test_exact_match_false_scores_like_the_legacy_matching(vcf, config, tmp_path, monkeypatch):
    """Opting out of exact_match scores exactly as before exact_match existed."""
    # GIVEN rank models with exact_match = false, and the scores genmod gives with the legacy
    # substring matching (a frozen copy of what it did before)
    config = opt_out_of_exact_match(config, tmp_path)
    with monkeypatch.context() as patch:
        patch.setattr(config_parser, "Plugin", FrozenLegacyPlugin)
        legacy_scores = score_with_rank_results(vcf, config, tmp_path)
    assert legacy_scores

    # WHEN scoring the same VCF with the current matching
    scores = score_with_rank_results(vcf, config, tmp_path)

    # THEN every variant gets the same rank score and category scores
    assert scores == legacy_scores


def test_legacy_comparison_detects_a_change_in_matching(tmp_path, monkeypatch):
    """The comparison above has teeth: it fails for a rank model with the default matching."""
    # GIVEN the legacy scores of the ClinVar VCF, scored with the default rank model
    with monkeypatch.context() as patch:
        patch.setattr(config_parser, "Plugin", FrozenLegacyPlugin)
        legacy_scores = score_with_rank_results(CLINVAR_VCF, EXACT_CONFIG, tmp_path)

    # WHEN scoring it with the current code, where the default matching is applied
    scores = score_with_rank_results(CLINVAR_VCF, EXACT_CONFIG, tmp_path)

    # THEN the conflicting variant is scored differently
    assert scores != legacy_scores
    assert scores[f"chr5:{CONFLICTING}"] != legacy_scores[f"chr5:{CONFLICTING}"]
