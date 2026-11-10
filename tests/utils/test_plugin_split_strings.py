from genmod.plugin.split_strings import split_strings

def test_simple_string():
    """Test how split_strings splits a simple string"""
    string = "1,2"
    splitters = [',']
    
    assert set(split_strings(string, splitters)) == set(['1','2'])

def test_simple_string_no_splitters():
    """Test how split_strings splits a simple string"""
    string = "1,2"
    splitters = []
    
    assert set(split_strings(string, splitters)) == set(['1,2'])

    
def test_two_splitters():
    """
    Test how split_strings splits a string with two splitters
    """
    
    string = "1:a,2:b"
    splitters = [',',':']
    
    assert set(split_strings(string, splitters)) == set(['1','a','2','b'])
    
def test_complex():
    """
    Test how split_strings splits a complex string
    """
    
    string = "a:1,b:2|3"
    splitters = [',',':','|']
    
    assert set(split_strings(string, splitters)) == set(['1','a','2','b','3'])