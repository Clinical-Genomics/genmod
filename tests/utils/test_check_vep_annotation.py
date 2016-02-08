from genmod.utils import check_vep_annotation

def test_get_none():
    """
    Test to get wrong annotation terms
    """
    
    vep_variant = {'vep_info':{ 
        'A':[
            {
                'Consequence': 'Two&Different',
                'Gene': 'ADK'
            }
            ]
        }
    }
    
    # The result should be empty since the terms do not exist
    assert check_vep_annotation(vep_variant) == set()

def test_get_annotation():
    """
    Test to get true annotation terms
    """
    
    vep_variant = {'vep_info':{ 
        'A':[
            {
                'Consequence': 'transcript_ablation&splice_donor_variant',
                'Gene': 'ADK'
            }
            ]
        }
    }
    
    # The result should be empty since the terms do not exist
    assert check_vep_annotation(vep_variant) == set(['ADK'])