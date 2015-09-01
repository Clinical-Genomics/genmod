from genmod.utils import generate_pairs

import pytest

def test_generate_pairs():
    """Test if generate pairs behave as suspected"""
    objects = [1,2]
    pairs = []
    for pair in generate_pairs(objects):
        pairs.append(pair)
    
    assert pairs == [(1,2)]

def test_non_iterator():
    """Test if generate pairs behave as suspected"""
    objects = 1
    pairs = []
    with pytest.raises(SyntaxError):
        for pair in generate_pairs(objects):
            pairs.append(pair)

def test_one_object():
    """Test if generate pairs behave as suspected"""
    objects = [1]
    pairs = []
    with pytest.raises(SyntaxError):
        for pair in generate_pairs(objects):
            pairs.append(pair)

def test_generate_multiple_pairs():
    """Test if generate pairs behave as suspected"""
    objects = [1,2,3,4]
    pairs = []
    for pair in generate_pairs(objects):
        pairs.append(pair)
    
    assert len(pairs) == 6
    assert pairs[0] == (1,2)
    assert pairs[1] == (1,3)
    assert pairs[-1] == (3,4)
