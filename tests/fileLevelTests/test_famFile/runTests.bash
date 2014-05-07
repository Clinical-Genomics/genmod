#!/bin/bash

# Testing whether non-standard FAM files are handled appropriately

# Call with FAM file
function run(){
echo "######################## ${1}"
cat ${1}
echo -e "\n\n"
../../../scripts/run_genmod.py -phased ${1} oneARvariant.vcf
echo -e "\n\n\n\n"
}


echo "SHOULD RUN THROUGH"
run standardTrio.fam


echo "SHOULD THROW EXCEPTION"
run standardTrio_extraHealthySister.fam
run standardTrio_missingFather.fam
run standardTrio_noRelationships.fam
run standardTrio_missingColForSon.fam
run standardTrio_extraColForSon.fam
run standardTrio_sexInversionOnParents.fam





