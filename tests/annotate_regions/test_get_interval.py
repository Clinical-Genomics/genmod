from genmod.annotate_regions.parse_annotations import get_interval

def test_get_interval():
    # GIVEN some coordinates and a symbol
    start = 1
    stop = 10
    symbol = 'first'
    
    # WHEN building an interval
    interval = get_interval(start, stop, symbol)
    
    # THEN the interval should have the right properties
    
    assert interval.begin == start
    assert interval.end == stop
    assert interval.data == symbol