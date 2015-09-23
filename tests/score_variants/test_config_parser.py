from genmod.score_variants import ConfigParser

CONFIG = "tests/fixtures/score_variants/genmod_example.ini"

# from genmod import logger
# from genmod.log import init_log
# init_log(logger, loglevel='DEBUG')


def test_config_parser():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)
    
    assert set(config_reader.plugins) == set(['1000G', 'CADD', 'GeneticModels','CLNSIG'])
    assert set(config_reader.categories.keys()) == set(['allele_frequency', 
                'deleteriousness', 'inheritance', 'clinical_significance'])
    assert set(config_reader.categories['allele_frequency'].keys()) == set(['category_aggregation', 'plugins'])
    assert set(config_reader.categories['allele_frequency']['plugins']) == set(['1000G'])

def test_get_score():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)
    
    variant = {'info_dict':{
        '1000GAF': '0.1',
        'CADD': '12'
    }}
    
    assert config_reader.plugins['1000G'].get_value(variant_dict=variant) == 0.1
    assert config_reader.plugins['CADD'].get_value(variant_dict=variant) == 12

    assert config_reader.score_functions['1000G'].get_score(0.01) == 1.0
    assert config_reader.score_functions['1000G'].get_score(0.001) == 2.0
    assert config_reader.score_functions['1000G'].get_score(None) == 3.0

def test_get_score_string():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)
    
    variant = {'info_dict':{
        '1000GAF': '0.1',
        'CADD': '12',
        'GeneticModels': '1:AD|AD_dn',
    }}
    
    assert config_reader.plugins['GeneticModels'].get_value(variant_dict=variant) == "AD"

    assert config_reader.score_functions['GeneticModels'].get_score("AD") == 3.0
    assert config_reader.score_functions['GeneticModels'].get_score("AD_dn") == 2.0
    assert config_reader.score_functions['GeneticModels'].get_score(None) == -12

def test_get_score_value():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)
    
    variant = {'info_dict':{
        '1000GAF': '0.1',
        'CADD': '12',
        'CLNSIG': '2',
    }}
    
    assert config_reader.plugins['CLNSIG'].get_value(variant_dict=variant) == 2

    assert config_reader.score_functions['CLNSIG'].get_score(2) == -1.0

def test_get_score_value_multiple_values():
    """Test the config parser"""
    config_reader = ConfigParser(CONFIG)
    
    variant = {'info_dict':{
        '1000GAF': '0.1',
        'CADD': '12',
        'CLNSIG': '2|5',
    }}

    assert config_reader.plugins['CLNSIG'].get_raw_entry(variant_dict=variant) == '2|5'
    assert config_reader.plugins['CLNSIG'].get_entry(variant_dict=variant) == ['2','5']
    
    assert config_reader.plugins['CLNSIG'].get_value(variant_dict=variant) == 5

    assert config_reader.score_functions['CLNSIG'].get_score(5) == 2
