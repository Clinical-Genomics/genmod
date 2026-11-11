import pytest
from genmod.score_variants import ConfigParser
from validate import ValidateError

CONFIG = "tests/fixtures/score_variants/genmod_example.ini"

# from genmod import logger
# from genmod.log import init_log
# init_log(logger, loglevel='DEBUG')


def test_config_parser():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)

    assert set(config_reader.plugins) == set(["1000G", "CADD", "GeneticModels", "CLNSIG"])
    assert set(config_reader.categories.keys()) == set(
        ["allele_frequency", "deleteriousness", "inheritance", "clinical_significance"]
    )
    assert set(config_reader.categories["allele_frequency"].keys()) == set(
        ["category_aggregation", "plugins"]
    )
    assert set(config_reader.categories["allele_frequency"]["plugins"]) == set(["1000G"])


def test_get_score():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)

    variant = {"info_dict": {"1000GAF": "0.1", "CADD": "12"}}

    assert config_reader.plugins["1000G"].get_value(variant_dict=variant) == 0.1
    assert config_reader.plugins["CADD"].get_value(variant_dict=variant) == 12

    assert config_reader.score_functions["1000G"].get_score(0.01) == 1.0
    assert config_reader.score_functions["1000G"].get_score(0.001) == 2.0
    assert config_reader.score_functions["1000G"].get_score(None) == 3.0


def test_get_score_string():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)

    variant = {
        "info_dict": {
            "1000GAF": "0.1",
            "CADD": "12",
            "GeneticModels": "1:AD|AD_dn",
        }
    }

    assert config_reader.plugins["GeneticModels"].get_value(variant_dict=variant) == "AD"

    assert config_reader.score_functions["GeneticModels"].get_score("AD") == 3.0
    assert config_reader.score_functions["GeneticModels"].get_score("AD_dn") == 2.0
    assert config_reader.score_functions["GeneticModels"].get_score(None) == -12


def test_get_score_value():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)

    variant = {
        "info_dict": {
            "1000GAF": "0.1",
            "CADD": "12",
            "CLNSIG": "2",
        }
    }

    assert config_reader.plugins["CLNSIG"].get_value(variant_dict=variant) == 2

    assert config_reader.score_functions["CLNSIG"].get_score(2) == -1.0


def test_get_score_value_multiple_values():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)

    variant = {
        "info_dict": {
            "1000GAF": "0.1",
            "CADD": "12",
            "CLNSIG": "2|5",
        }
    }

    assert config_reader.plugins["CLNSIG"].get_raw_entry(variant_dict=variant) == "2|5"
    assert config_reader.plugins["CLNSIG"].get_entry(variant_dict=variant) == ["2", "5"]

    assert config_reader.plugins["CLNSIG"].get_value(variant_dict=variant) == 5

    assert config_reader.score_functions["CLNSIG"].get_score(5) == 2


STRING_CLNSIG_INI = (
    "[Version]\n  version = 1.0\n  name = clnsig_string\n"
    "[Categories]\n  [[Clinical_significance]]\n    category_aggregation = max\n"
    "[CLNSIG]\n  field = INFO\n  data_type = string\n"
    "  category = Clinical_significance\n  record_rule = max\n"
    "  separators = '|', '/', ',',\n  info_key = CLNSIG\n{extra}"
    "  [[not_reported]]\n    score = 0\n"
    "  [[pathogenic]]\n    score = 5\n    priority = 5\n    string = 'Pathogenic'\n"
    "  [[likely_pathogenic]]\n    score = 2\n    priority = 4\n    string = 'Likely_pathogenic'\n"
)


def get_clnsig_scorer(tmp_path, extra=""):
    ini = tmp_path / "clnsig_string.ini"
    ini.write_text(STRING_CLNSIG_INI.format(extra=extra))
    config_reader = ConfigParser(str(ini))
    plugin = config_reader.plugins["CLNSIG"]
    score_function = config_reader.score_functions["CLNSIG"]

    def score(clnsig):
        value = plugin.get_value(variant_dict={"info_dict": {"CLNSIG": clnsig}})
        return score_function.get_score(value)

    return score


def test_string_clnsig_exact_match_scores_conflicting_below_pathogenic(tmp_path):
    """Issue #87: with exact_match, conflicting classifications must not score as pathogenic."""
    # GIVEN a string CLNSIG rank model with exact_match = true
    score = get_clnsig_scorer(tmp_path, extra="  exact_match = true\n")

    # WHEN scoring pathogenic, likely pathogenic and conflicting classifications
    # THEN each gets its own score and the conflicting one is not scored as pathogenic
    assert score("Pathogenic") == 5
    assert score("Likely_pathogenic") == 2
    assert score("Conflicting_classifications_of_pathogenicity") == 0


@pytest.mark.parametrize("extra", ["", "  exact_match = false\n"])
def test_string_clnsig_default_is_legacy_substring_matching(tmp_path, extra):
    """Without exact_match the scores of existing rank models do not change."""
    # GIVEN a string CLNSIG rank model without exact_match (or with exact_match = false)
    score = get_clnsig_scorer(tmp_path, extra=extra)

    # WHEN scoring likely pathogenic and conflicting classifications
    # THEN substring matching still scores both as pathogenic
    assert score("Likely_pathogenic") == 5
    assert score("Conflicting_classifications_of_pathogenicity") == 5


def test_exact_match_must_be_a_boolean(tmp_path):
    # GIVEN a rank model where exact_match is not true or false
    # WHEN parsing the config
    # THEN a ValidateError is raised
    with pytest.raises(ValidateError):
        get_clnsig_scorer(tmp_path, extra="  exact_match = maybe\n")
