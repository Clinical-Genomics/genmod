import pytest
from typing import Union, Dict, Any
from tempfile import NamedTemporaryFile
from configobj import ConfigObj

from genmod.score_variants import ConfigParser
from genmod.vcf_tools import get_info_dict
from genmod.score_variants.score_variant import get_category_score


class ConfigObjWithNamedTemporaryFile(ConfigObj):
    """
    Class that wraps a NamedTemporaryFile inside a ConfigObject
    """
    def __init__(self, named_temporary_file: NamedTemporaryFile, *args, **kwargs):
        """
        Args:
            named_temporary_file: A tmp file that will eventually contain the written config
        """
        super().__init__(*args, **kwargs)
        self._file_pointer: NamedTemporaryFile = named_temporary_file
        self.filename: str = named_temporary_file.name

    def add_category(self,
                     name: str,
                     aggregation_mode: str):
        """
        Add a category to config file
        Args:
            name: Category name
            aggregation_mode: mode of score aggregation for this category
        Returns:
            self for chaining
        """
        self['Categories'][name]: Dict = {}
        self['Categories'][name]['category_aggregation'] = aggregation_mode
        self.write()
        return self

    def add_plugin(self,
                   category: str,
                   plugin_name: str,
                   score_categories: Dict[str, Dict[str, float]],
                   info_key: str,
                   not_reported_score: float = 0.0,
                   field: str = 'INFO',
                   data_type: str = 'float',
                   record_rule: str = 'max',
                   separators: str = ',',
                   ):
        """
        Add plugin definition to scoring config.
        Arguments is 1:1 mapped to scoring config file.

        For the remaining arguments, please refer to the scoring config spec.
        Args:
            category: Plugin category this plugin belong sto
            plugin_name: Name of plugin
            score_categories: A map containing a range, example:
              {'rare': {'score': 2.0, 'lower': 0.01, 'upper': 0.05}}
            which maps to:
              [[rare]]
                score = 1
                lower = 0.01
                upper = 1.0
            info_key: Mapping key that holds the plugin input value used for scoring

        Returns:
            self for chaining
        """
        self[plugin_name]: Dict = {}
        self[plugin_name]['field'] = field
        self[plugin_name]['data_type'] = data_type
        self[plugin_name]['category'] = category
        self[plugin_name]['record_rule'] = record_rule
        self[plugin_name]['separators'] = separators
        self[plugin_name]['info_key'] = info_key

        self[plugin_name]['not_reported']: Dict = {}
        self[plugin_name]['not_reported']['score'] = not_reported_score

        for score_category in score_categories.keys():
            self[plugin_name][score_category]: Dict = {}

        for score_category, score_range in score_categories.items():
            for score_range_key, score_range_value in score_range.items():
                self[plugin_name][score_category][score_range_key] = score_range_value
        self.write()
        return self


def set_info_dict_in_variant(variant: Dict[str, str]) -> Dict[str, Union[str, Any]]:
    """
    Create info_dict attribute required by get_category_score.
    Args:
        variant: A variant dict, example:
            variant = {'CHROM': '1',
                       'POS': '1',
                       'INFO': 'plugin_a=0.01;plugin_b=1.0'
                       }
    Returns:
        variant with ['info_dict'] attribute
    """
    variant['info_dict']: Dict[str, str] = get_info_dict(variant['INFO'])
    return variant


@pytest.fixture
def config_file() -> ConfigObjWithNamedTemporaryFile:
    """
    Provide a scoring config file for tests
    Returns:
        A minimal scoring config file as ConfigObjWithNamedTemporaryFile
    """
    named_temporary_file: NamedTemporaryFile = NamedTemporaryFile(dir='/tmp')
    config: ConfigObjWithNamedTemporaryFile = \
        ConfigObjWithNamedTemporaryFile(named_temporary_file=named_temporary_file)
    config['Version'] = {}
    config['Version']['version'] = 0.1
    config['Version']['name'] = 'genmod example'
    config['Categories'] = {}
    return config


@pytest.fixture
def variant_0():
    variant = {
        'CHROM': '1',
        'POS': '1',
        'INFO': 'plugin_a=0.01;plugin_b=5.0'
    }
    variant = set_info_dict_in_variant(variant)
    return variant


@pytest.mark.parametrize('aggregation_mode,expected',
                         [('sum', (2, 0, 2)),
                          ('max', (2, 0, 2)),
                          ('min', (2, 0, 2))])
def test_single_category_single_plugin(config_file, variant_0, aggregation_mode, expected):
    """
    Test with a single category.
    """
    # GIVEN a scoring config in mode aggregation_mode
    config_file.add_category('CategoryA', aggregation_mode)
    config_file.add_plugin('CategoryA',
                           'PluginA',
                           {'rare': {'score': 2.0, 'lower': 0.01, 'upper': 0.05}},
                           'plugin_a')
    config_parser: ConfigParser = ConfigParser(config_file.filename)
    # WHEN scoring a variant
    score, min, max = get_category_score(variant=variant_0,
                                         category='CategoryA',
                                         config_parser=config_parser)
    # THEN expect proper score and min max bounds
    assert score == expected[0]
    assert min == expected[1]
    assert max == expected[2]


@pytest.mark.parametrize('aggregation_mode,expected',
                         [('sum', (2, -5, 2)),
                          ('max', (2, -5, 2)),
                          ('min', (2, -5, 2))])
def test_single_category_single_plugin_notreported(config_file, variant_0, aggregation_mode, expected):
    """
    Test with a single category, with a custom not reported score.
    """
    # GIVEN a scoring config in mode aggregation_mode
    config_file.add_category('CategoryA', aggregation_mode)
    config_file.add_plugin('CategoryA',
                           'PluginA',
                           {'rare': {'score': 2.0, 'lower': 0.01, 'upper': 0.05}},
                           'plugin_a',
                           not_reported_score=-5)
    config_parser: ConfigParser = ConfigParser(config_file.filename)
    # WHEN scoring a variant
    score, min, max = get_category_score(variant=variant_0,
                                         category='CategoryA',
                                         config_parser=config_parser)
    # THEN expect proper score and min max bounds
    assert score == expected[0]
    assert min == expected[1]
    assert max == expected[2]


@pytest.mark.parametrize('aggregation_mode,expected',
                         [('sum', (3, 0, 3)),
                          ('max', (2, 0, 2)),
                          ('min', (1, 0, 1))])
def test_single_category_two_plugin(config_file, variant_0, aggregation_mode, expected):
    """
    Test two categories.
    """
    # GIVEN a scoring config in mode aggregation_mode
    config_file.add_category('CategoryA', aggregation_mode)
    config_file.add_plugin('CategoryA',
                           'PluginA',
                           {'rare': {'score': 2.0, 'lower': 0.01, 'upper': 0.05}},
                           'plugin_a')
    config_file.add_plugin('CategoryA',
                           'PluginB',
                           {'rare': {'score': 1.0, 'lower': 4.0, 'upper': 6.0}},
                           'plugin_b')
    config_parser: ConfigParser = ConfigParser(config_file.filename)
    # WHEN scoring a variant
    score, min, max = get_category_score(variant=variant_0,
                                         category='CategoryA',
                                         config_parser=config_parser)
    # THEN expect proper score and min max bounds
    assert score == expected[0]
    assert min == expected[1]
    assert max == expected[2]


@pytest.mark.parametrize('aggregation_mode,expected',
                         [('sum', (-3, -5, 2)),
                          ('max', (2, 0, 2)),
                          ('min', (-5, -5, 0))])
def test_single_category_two_plugin(config_file, variant_0, aggregation_mode, expected):
    """
    Test two categories, negative score range.
    """
    # GIVEN a scoring config in mode aggregation_mode
    config_file.add_category('CategoryA', aggregation_mode)
    config_file.add_plugin('CategoryA',
                           'PluginA',
                           {'rare': {'score': 2, 'lower': 0.01, 'upper': 0.05}},
                           'plugin_a')
    config_file.add_plugin('CategoryA',
                           'PluginB',
                           {'rare': {'score': -5, 'lower': 4.0, 'upper': 6.0}},
                           'plugin_b')
    config_parser: ConfigParser = ConfigParser(config_file.filename)
    # WHEN scoring a variant
    score, min, max = get_category_score(variant=variant_0,
                                         category='CategoryA',
                                         config_parser=config_parser)
    # THEN expect proper score and min max bounds
    assert score == expected[0]
    assert min == expected[1]
    assert max == expected[2]


@pytest.mark.parametrize('aggregation_mode,expected',
                         [('sum', (2, 0, 2)),
                          ('max', (2, 0, 2)),
                          ('min', (2, 0, 2))])
def test_multi_category_two_plugins(config_file, variant_0, aggregation_mode, expected):
    """
    Test with a single category, shall ignore other category
    """
    # GIVEN a scoring config in mode aggregation_mode
    config_file.add_category('CategoryA', aggregation_mode)
    config_file.add_category('CategoryB', aggregation_mode)
    config_file.add_plugin('CategoryA',
                           'PluginA',
                           {'rare': {'score': 2.0, 'lower': 0.01, 'upper': 0.05}},
                           'plugin_a')
    config_file.add_plugin('CategoryB',
                           'PluginB',
                           {'rare': {'score': 1.0, 'lower': 4.0, 'upper': 6.0}},
                           'plugin_b')
    config_parser: ConfigParser = ConfigParser(config_file.filename)
    # WHEN scoring a variant, expect only from category A
    score, min, max = get_category_score(variant=variant_0,
                                         category='CategoryA',
                                         config_parser=config_parser)
    # THEN expect proper score and min max bounds
    assert score == expected[0]
    assert min == expected[1]
    assert max == expected[2]

@pytest.mark.parametrize('aggregation_mode,expected',
                         [('sum', (10, 0, 10)),
                          ('max', (10, 0, 10)),
                          ('min', (10, 0, 10))])
def test_single_category_single_plugin_multirange(config_file, variant_0, aggregation_mode, expected):
    """
    Test with a single category, having multiple score ranges.
    """
    # GIVEN a scoring config in mode aggregation_mode
    config_file.add_category('CategoryA', aggregation_mode)
    config_file.add_plugin('CategoryA',
                           'PluginA',
                           {'rare': {'score': 10.0, 'lower': 0.01, 'upper': 0.05},
                            'common': {'score': 1.0, 'lower': 0.06, 'upper': 1.0}},
                           'plugin_a')
    config_parser: ConfigParser = ConfigParser(config_file.filename)
    # WHEN scoring a variant
    score, min, max = get_category_score(variant=variant_0,
                                         category='CategoryA',
                                         config_parser=config_parser)
    # THEN expect proper score and min max bounds
    assert score == expected[0]
    assert min == expected[1]
    assert max == expected[2]