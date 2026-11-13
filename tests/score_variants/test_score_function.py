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
    not_reported_score = -100
    reported_score = 10
    score_function.set_not_reported(not_reported_score)
    score_function.set_reported(reported_score)

    # WHEN accessing min max plugin score
    # THEN expect the proper min max values
    assert score_function.score_min == not_reported_score
    assert score_function.score_max == reported_score

    # WHEN asking for a missing value, THEN not reported is flagged and scored
    assert score_function.get_score(None) == not_reported_score

    # WHEN asking for a set value, THEN reported is flagged and scored
    assert score_function.get_score("c'est une pipe") == reported_score


def test_interval_lower_bound_is_inclusive_and_upper_bound_is_exclusive():
    """Intervals are [lower, upper): see genmod issue #206."""
    # GIVEN a score function with two adjacent intervals sharing the boundary 0.02
    score_function = ScoreFunction(match_type="float")
    score_function.add_interval(lower=0.02, upper=1.1, score=-12)
    score_function.add_interval(lower=0.005, upper=0.02, score=1)

    # WHEN scoring the shared boundary value
    # THEN the interval with lower = 0.02 is chosen, not the one with upper = 0.02
    assert score_function.get_score("0.02") == -12

    # AND the value just below the boundary is in the other interval
    assert score_function.get_score("0.0199") == 1

    # AND the lowest value of an interval is included
    assert score_function.get_score("0.005") == 1

    # AND the upper bound of the highest interval is excluded
    assert score_function.get_score("1.1") == 0


def test_intervals_can_be_added_in_any_order():
    # GIVEN a score function with intervals added out of order
    score_function = ScoreFunction(match_type="float")
    score_function.add_interval(lower=10, upper=20, score=3)
    score_function.add_interval(lower=0, upper=5, score=1)
    score_function.add_interval(lower=5, upper=10, score=2)

    # WHEN scoring values in each interval
    # THEN the score of the interval containing the value is returned
    assert score_function.get_score("0") == 1
    assert score_function.get_score("4.99") == 1
    assert score_function.get_score("5") == 2
    assert score_function.get_score("10") == 3
    assert score_function.get_score("19.99") == 3


def test_value_in_gap_between_intervals_gets_no_score():
    # GIVEN two intervals with a gap between them
    score_function = ScoreFunction(match_type="float")
    score_function.add_interval(lower=0, upper=5, score=1)
    score_function.add_interval(lower=10, upper=15, score=2)

    # WHEN scoring values in the gap and outside the intervals
    # THEN they score 0
    assert score_function.get_score("7") == 0
    assert score_function.get_score("-1") == 0
    assert score_function.get_score("100") == 0


def test_single_interval():
    # GIVEN a score function with one interval
    score_function = ScoreFunction(match_type="integer")
    score_function.add_interval(lower=1, upper=3, score=7)

    # THEN values inside score, the upper bound does not
    assert score_function.get_score("1") == 7
    assert score_function.get_score("2") == 7
    assert score_function.get_score("3") == 0
    assert score_function.get_score("0") == 0


def test_negative_and_float_intervals():
    # GIVEN intervals with negative and fractional bounds
    score_function = ScoreFunction(match_type="float")
    score_function.add_interval(lower=-1.5, upper=-0.5, score=-4)
    score_function.add_interval(lower=-0.5, upper=0.25, score=2)

    # THEN the shared bound belongs to the interval where it is the lower bound
    assert score_function.get_score("-1.5") == -4
    assert score_function.get_score("-0.5") == 2
    assert score_function.get_score("0.25") == 0


@pytest.mark.parametrize("lower, upper", [(1, 2), (0, 10), (-5, 1.5), (4, 6)])
def test_overlapping_intervals_are_rejected(lower, upper):
    # GIVEN a score function with the interval [1, 5)
    score_function = ScoreFunction(match_type="float")
    score_function.add_interval(lower=1, upper=5, score=1)

    # WHEN adding an interval that overlaps it
    # THEN an error is raised instead of silently picking one of the scores
    with pytest.raises(ValueError, match="overlap"):
        score_function.add_interval(lower=lower, upper=upper, score=2)


def test_adjacent_intervals_are_allowed():
    # GIVEN the interval [1, 5)
    score_function = ScoreFunction(match_type="float")
    score_function.add_interval(lower=1, upper=5, score=1)

    # WHEN adding intervals that only touch it
    score_function.add_interval(lower=5, upper=6, score=2)
    score_function.add_interval(lower=0, upper=1, score=3)

    # THEN no error is raised and the intervals are scored separately
    assert score_function.get_score("5") == 2
    assert score_function.get_score("0") == 3


@pytest.mark.parametrize("lower, upper", [(5, 5), (6, 5)])
def test_interval_without_length_is_rejected(lower, upper):
    # GIVEN a score function
    score_function = ScoreFunction(match_type="float")

    # WHEN adding an interval where lower >= upper
    # THEN an error is raised
    with pytest.raises(ValueError, match="lower < upper"):
        score_function.add_interval(lower=lower, upper=upper, score=1)


def test_score_range_includes_all_interval_scores_and_not_reported():
    # GIVEN a score function with intervals and a not reported score
    score_function = ScoreFunction(match_type="float")
    score_function.set_not_reported(-3)
    score_function.add_interval(lower=0, upper=1, score=0.5)
    score_function.add_interval(lower=1, upper=2, score=2)

    # THEN the range holds every interval score and the not reported score
    assert sorted(score_function.score_range) == [-3.0, 0.5, 2.0]


@pytest.mark.parametrize(
    "match_type, value", [("integer", "abc"), ("float", "abc"), ("integer", "1.5")]
)
def test_non_numeric_value_raises(match_type, value):
    # GIVEN a numeric score function
    score_function = ScoreFunction(match_type=match_type)
    score_function.add_interval(lower=0, upper=10, score=1)

    # WHEN scoring a value that is not a number of that type
    # THEN a clear error is raised
    with pytest.raises(ValueError, match="Value has to be a number"):
        score_function.get_score(value)


def test_value_lookup_scores_matching_value_and_zero_otherwise():
    # GIVEN a score function with value rules, given as strings like in a config file
    score_function = ScoreFunction(match_type="integer")
    score_function.add_value(value="1", score="2.5")
    score_function.add_value(value="3", score="7")

    # THEN matching values get their (truncated) score and others get 0
    assert score_function.get_score("1") == 2
    assert score_function.get_score("3") == 7
    assert score_function.get_score("2") == 0


def test_flag_and_missing_value_scores():
    # GIVEN a flag score function
    score_function = ScoreFunction(match_type="flag")
    score_function.set_reported(5)
    score_function.set_not_reported(-2)

    # THEN a present flag gets the reported score and a missing value the not reported score
    assert score_function.get_score(True) == 5
    assert score_function.get_score(None) == -2
