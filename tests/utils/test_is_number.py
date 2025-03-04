from genmod.utils import is_number


def test_int():
    """Test if is_number behave as suspected"""
    obj = 2
    assert is_number(obj) is True


def test_float():
    """Test if is_number behave as suspected"""
    obj = 2.5
    assert is_number(obj) is True


def test_non_number():
    """Test if is_number behave as suspected"""
    obj = "a"
    assert is_number(obj) is False


def test_str_int():
    """Test if is_number behave as suspected"""
    obj = "1"
    assert is_number(obj) is True


def test_str_float():
    """Test if is_number behave as suspected"""
    obj = "1.3"
    assert is_number(obj) is True
