from genmod.annotate_models.models import check_X_recessive
from genmod.vcf_tools import Genotype

from ped_parser import FamilyParser

FAMILY_FILE = "tests/fixtures/recessive_trio.ped"

def get_family(family_file = None, family_lines = None):
    """Return a family object
    
    """
    family = None
    if family_file:
        family = FamilyParser(open(family_file, 'r'))
    elif family_lines:
        family = FamilyParser(family_lines)
        
    return family



################# Test affected ###############
def test_x_affected_recessive_male():
    """Test a sick male
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    
    assert check_X_recessive(
        variant = recessive_variant,
        family = family
    ) == True

def test_x_affected_recessive_female():
    """Test a sick heterozygote female
    
    Females needs to bo hom alt to follow pattern
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t2\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    
    assert check_X_recessive(
        variant = recessive_variant,
        family = family
    ) == False

def test_x_affected_homozygote_male():
    """Test an affected homozygote male"""
    
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    homozygote_variant = {'genotypes': {}}
    homozygote_variant['genotypes']['proband'] = Genotype(**{'GT':'1/1'})
    
    assert check_X_recessive(
        variant = homozygote_variant,
        family = family
    ) == True

def test_x_affected_homozygote_female():
    """Test an affected homozygote male"""
    
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t2\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    homozygote_variant = {'genotypes': {}}
    homozygote_variant['genotypes']['proband'] = Genotype(**{'GT':'1/1'})
    
    assert check_X_recessive(
        variant = homozygote_variant,
        family = family
    ) == True

def test_x_affected_male_ref_call():
    """Test an affected ref call male"""
    
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    homozygote_variant = {'genotypes': {}}
    homozygote_variant['genotypes']['proband'] = Genotype(**{'GT':'0/0'})
    
    assert check_X_recessive(
        variant = homozygote_variant,
        family = family
    ) == False

def test_x_affected_female_ref_call():
    """Test an affected ref call male"""
    
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t2\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    homozygote_variant = {'genotypes': {}}
    homozygote_variant['genotypes']['proband'] = Genotype(**{'GT':'0/0'})
    
    assert check_X_recessive(
        variant = homozygote_variant,
        family = family
    ) == False
    

def test_x_affected_no_call_male():
    """Test a sick male with no gt call
    
    This should be true since there is no information that contradicts the model
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    no_call_variant = {'genotypes': {}}
    no_call_variant['genotypes']['proband'] = Genotype(**{'GT':'./.'})
    
    assert check_X_recessive(
        variant = no_call_variant,
        family = family
    ) == True

def test_x_affected_no_call_male_strict():
    """Test a sick male with no gt call
    
    This should not be true since we allways need 'proof'
    for an inheritance pattern if strict mode.
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t2\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    no_call_variant = {'genotypes': {}}
    no_call_variant['genotypes']['proband'] = Genotype(**{'GT':'./.'})
    
    assert check_X_recessive(
        variant = no_call_variant,
        family = family,
        strict = True
    ) == False

############### Test healthy ##############

def test_x_healthy_recessive_male():
    """Test a healthy recessive male
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t1\t1\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    
    assert check_X_recessive(
        variant = recessive_variant,
        family = family
    ) == False

def test_x_healthy_recessive_female():
    """Test a healthy heterozygote female
    
    Females needs to bo hom alt to follow pattern
    """
    family_lines = [
        "#FamilyID\tSampleID\tFather\tMother\tSex\tPhenotype\n",
        "1\tproband\t0\t0\t2\t1\n"
    ]
    
    family = get_family(family_lines=family_lines)
    
    recessive_variant = {'genotypes': {}}
    recessive_variant['genotypes']['proband'] = Genotype(**{'GT':'0/1'})
    
    assert check_X_recessive(
        variant = recessive_variant,
        family = family
    ) == True
