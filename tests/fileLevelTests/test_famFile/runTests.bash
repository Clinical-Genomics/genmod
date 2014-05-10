#!/bin/bash


# Testing whether non-standard FAM files are handled appropriately

# Call with FAM file
function run(){
echo "######################## ${1}"
cat ${1}
echo -e "\n"
../../../scripts/run_genmod.py -phased ${1} oneARvariant.vcf
echo -e "\n\n\n"
}


echo "FOLLOWING SHOULD RUN"
run standardTrio.fam

echo "FOLLOWING SHOULD THROW EXCEPTION"
run standardTrio_extraHealthySister.fam

echo "FOLLOWING SHOULD THROW EXCEPTION"
run standardTrio_missingFather.fam

echo "FOLLOWING SHOULD THROW EXCEPTION"
run standardTrio_noRelationships.fam

echo "FOLLOWING SHOULD THROW EXCEPTION"
run standardTrio_missingColForSon.fam

echo "FOLLOWING SHOULD THROW EXCEPTION"
run standardTrio_extraColForSon.fam

echo "FOLLOWING SHOULD THROW EXCEPTION"
run standardTrio_sexInversionOnParents.fam





