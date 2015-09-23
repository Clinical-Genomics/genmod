from genmod.score_variants import ScoreFunction

def test_string_score():
    """Test the score function with a string function"""
    
    not_reported_score = 4
    string_dict = {
        'hello': 1,
        'world': 2
    }

    score_function = ScoreFunction(match_type = 'string')
    
    for key in string_dict:
        score_function.add_string_rule(key, string_dict[key])

    score_function.set_not_reported(not_reported_score)
    
    assert score_function.get_score('hello') == 1
    assert score_function.get_score('world') == 2
    assert score_function.get_score(None) == not_reported_score
    assert score_function.get_score('non_existing') == 0


def test_int_score():
    """Test the score function with a integer"""
    
    not_reported_score = 1
    
    score_function = ScoreFunction(match_type = 'integer')
    score_function.add_interval(lower=0, upper=10, score=1)
    score_function.add_interval(lower=10, upper=15, score=2)
    score_function.add_interval(lower=15, upper=20, score=3)

    score_function.set_not_reported(not_reported_score)
    
    assert score_function.get_score(3) == 1
    assert score_function.get_score(12) == 2
    assert score_function.get_score(15) == 3
    assert score_function.get_score(19) == 3
    assert score_function.get_score(None) == not_reported_score
    assert score_function.get_score(-3) == 0

def test_eq_score():
    """Test the score function when the score found should be returned"""
    
    not_reported_score = 1
    
    score_function = ScoreFunction(match_type='integer', equal=True)
    
    score_function.set_not_reported(not_reported_score)
    
    assert score_function.get_score(3) == 3
    assert score_function.get_score(12) == 12
    assert score_function.get_score(None) == not_reported_score
