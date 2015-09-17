from . import get_model_score

def make_print_version(variant, families):
    """
    Get the variants ready for printing
    
    This function collects the annotations added and merge them in the INFO
    dict.
    
    Arguments:
        variant (dict): A variant dictionary
    
    """
    
    variant_id = variant['variant_id']
    vcf_info = variant['INFO'].split(';')
    
    feature_list = variant.get('annotation', set())
    
    # variant[compounds] is a dictionary with family id as keys and a set of compounds as values
    compounds = variant.get('compounds', dict())
    # Here we store the compound strings that should be added to the variant:
    family_compound_strings = []

    # We need to check if compounds have already been annotated.
    if 'Compounds' not in variant['info_dict']:
    
        for family_id in compounds:
            compound_string = ''
            compound_set = compounds[family_id]
            #We do not want reference to itself as a compound:
            compound_set.discard(variant_id)
            # If there are any compounds for the family:
            if compounds[family_id]:
                compound_string = '|'.join(compound_set)
                family_compound_strings.append(':'.join([family_id, compound_string]))

            if len(family_compound_strings) > 0:
                vcf_info.append('Compounds=' + ','.join(family_compound_strings))
    
    # Check if any genetic models are followed
    if 'GeneticModels' not in variant['info_dict']:
        # Here we store the compound strings that should be added to the variant:
        family_model_strings = []
        model_scores = {}
        genetic_models = variant.get('inheritance_models', {})
        for family_id in genetic_models:
            model_string = ''
            model_list = []
            for model in genetic_models[family_id]:
                if genetic_models[family_id][model]:
                    model_list.append(model)
                model_string = '|'.join(model_list)
            if len(model_list) > 0:
                family_model_strings.append(':'.join(
                                            [family_id, model_string]))
                
                model_scores[family_id] = str(
                    get_model_score(families[family_id].individuals, variant))
        
        if len(family_model_strings) > 0:
            vcf_info.append(
                        'GeneticModels={0}'.format( 
                        ','.join(family_model_strings)))

            model_score_list = []
            for family_id in model_scores:
                if model_scores[family_id]:
                    if float(model_scores[family_id]) > 0:
                        model_score_list.append(
                            ':'.join(
                                        [
                                            family_id, 
                                            model_scores[family_id]
                                        ]
                                    )
                                )
            if len(model_score_list) > 0:
                vcf_info.append(
                            'ModelScore=' +  
                            ','.join(model_score_list)
                            )
    
    
    
    variant['INFO'] = ';'.join(vcf_info)

    return variant
