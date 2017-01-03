from intervaltree import IntervalTree
from genmod.annotate_regions.parse_annotations import build_region_trees

def test_build_region_trees(bed_lines):
    # GIVEN some lines in the bed format
    # WHEN building region trees
    region_trees = build_region_trees(bed_lines, padding=0)
    # THEN assert there are interval trees in the result
    assert len(region_trees) > 0
    for chrom in region_trees:
        assert isinstance(region_trees[chrom], IntervalTree)


def test_query_region_trees():
    # GIVEN some bed formated lines
    lines = [
        "13\t1\t100\tHMGA1P6\tENSG00000233440\n",
        "13\t50\t200\tRNY3P4\tENSG00000207157\n",
        "13\t300\t1000\tLINC00362\tENSG00000229483\n"
    ]
    # WHEN building region trees
    region_trees = build_region_trees(lines, padding=0)
    
    # THEN make sure that we get the correct info when querying
    interval_tree = region_trees['13']
    
    # There should be two results here
    assert len(interval_tree[50]) == 2
    
    result = interval_tree[100]
    #Intervals are half opened
    assert len(result) == 1
    for interval in result:
        assert interval.data == 'RNY3P4'

    result = interval_tree[100:100]
    #Intervals without lenght does not match anythin
    assert len(result) == 0
        
    result = interval_tree[250:600]
    for interval in result:
        assert interval.data == 'LINC00362'
    
    # Test none overlapping
    result = interval_tree[10000:11000]
    assert not result
    assert isinstance(result, set)


# def test_build_all(ensembl_file):
#     with open(ensembl_file, 'r') as ensembl_handle:
#         region_trees = build_region_trees(ensembl_handle, padding=4000)
