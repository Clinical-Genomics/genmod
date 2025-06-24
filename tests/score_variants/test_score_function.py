import pytest
from genmod.score_variants import ScoreFunction


def test_string_score():
    """Test the score function with a string function"""

    not_reported_score = 4
    string_dict = {"hello": 1, "world": 2}

    score_function = ScoreFunction(match_type="string")

    for key in string_dict:
        score_function.add_string_rule(key, string_dict[key])

    score_function.set_not_reported(not_reported_score)

    assert score_function.get_score("hello") == 1
    assert score_function.get_score("world") == 2
    assert score_function.get_score(None) == not_reported_score
    assert score_function.get_score("non_existing") == 0


def test_int_score():
    """Test the score function with a integer"""

    not_reported_score = 1

    score_function = ScoreFunction(match_type="integer")
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

    score_function = ScoreFunction(match_type="integer", equal=True)

    score_function.set_not_reported(not_reported_score)

    assert score_function.get_score(3) == 3
    assert score_function.get_score(12) == 12
    assert score_function.get_score(None) == not_reported_score


def test_score_mode_user_defined_range():
    """Test score mode with user defined override."""

    # GIVEN a score function with user defined rank score override (missing plugin defined min-max scores)
    score_function = ScoreFunction(match_type="integer", equal=True)

    with pytest.raises(ValueError) as error:
        # WHEN trying to get the score range
        # THEN expect this to trigger an ValueError
        _ = score_function.score_range
    assert "User supplied score values does not have a known score range" in str(error.value)

    # GIVEN a score function with user defined rank score override (missing plugin defined min-max scores)
    score_function = ScoreFunction(match_type="integer")
    score_function.set_equal()
    with pytest.raises(ValueError) as error:
        # WHEN trying to get the score range
        # THEN expect this to trigger a ValueError
        _ = score_function.score_range
    assert "User supplied score values does not have a known score range" in str(error.value)


def test_score_mode_invalid_double_lookup():
    """Test ScoreFunction sanity check when using multiple scoring maps."""
    # GIVEN a score function with multiple maps
    score_function = ScoreFunction(match_type="integer")
    score_function.set_not_reported(-100)
    score_function.add_interval(-10, -1, 0.1)
    score_function.add_value(0, 0)
    with pytest.raises(ValueError) as error:
        # WHEN trying to get a plugin min or max bounds
        # THEN expect this to trigger a ValueError
        _ = score_function.score_min
    assert "Unable to accurately determine what mapping to use for determining score range" in str(
        error.value
    )


def test_score_mode_tree_lookup():
    """Test ScoreFunctions min max bounds property."""
    # GIVEN a score function with a tree (range) map
    score_function = ScoreFunction(match_type="integer")
    score_function.set_not_reported(-100)
    score_function.add_interval(-10, -1, 0.1)
    score_function.add_interval(0, 10, 0.5)
    score_function.add_interval(11, 15, 0.9)
    # WHEN accessing min max plugin score
    # THEN expect the proper min max values
    assert score_function.score_min == -100.0
    assert score_function.score_max == 0.9


def test_score_mode_value_lookup():
    """Test ScoreFunctions min max bounds property."""
    # GIVEN a score function with value dict map
    score_function = ScoreFunction(match_type="value")
    score_function.set_not_reported(-100.0)
    score_function.add_value(0, 0)
    score_function.add_value(1, 1)
    score_function.add_value(2, 3)
    score_function.add_value(3, 5)
    # WHEN accessing min max plugin score
    # THEN expect the proper min max values
    assert score_function.score_min == -100
    assert score_function.score_max == 5


def test_score_mode_string_lookup():
    """Test ScoreFunctions min max bounds property."""
    # GIVEN a score function with string dict map
    score_function = ScoreFunction(match_type="string")
    score_function.set_not_reported(-100)
    score_function.add_string_rule("foo", 0)
    score_function.add_string_rule("bar", 1)
    score_function.add_string_rule("0xdead", 2)
    # WHEN accessing min max plugin score
    # THEN expect the proper min max values
    assert score_function.score_min == -100
    assert score_function.score_max == 2

def test_score_mode_flag_lookup():
    """Test ScoreFunctions min max bounds property."""
    # GIVEN a score function with a flag reported/not_reported score
    score_function = ScoreFunction(match_type="flag")
    score_function.set_not_reported(-100)
    score_function.set_reported(10)

    # WHEN accessing min max plugin score
    # THEN expect the proper min max values
    assert score_function.score_min == -100
    assert score_function.score_max == 10
