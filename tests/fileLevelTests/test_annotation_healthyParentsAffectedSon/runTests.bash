#!/bin/bash

#set -x

# Annotation of VCF file with genotypes conforming to different models

echo "############# One call to update the gene and exon db"
sudo ../../../scripts/run_genmod.py -s -at gtf -an ../Homo_sapiens.GRCh37.75_partChr1andX.gtf healthyParentsAffectedSon.fam AD_TP.vcf | grep -v "##"




# Call with VCF file
# And any extra options eg "-phased" or "-phased -v" 
function run(){
echo "######################## $*"
# Running with GTF file so that testing of AR_comp is possible (default is just exonic)
# FAM file is fixed
../../../scripts/run_genmod.py $2 healthyParentsAffectedSon.fam $1 | grep -v "##"
}

<<COMMENT
echo "#######################################################################################"
echo "################### True positives - should see genmod add an annotation ##############"
echo "#######################################################################################"
run AD_TP.vcf
run AR_TP.vcf
run XD_TP.vcf
run XR_TP.vcf
COMMENT

run ARcomp.vcf | column -t
run ARcomp.vcf "-phased"  | column -t # calling the phased option is more restrictive on TP

<<COMMENT
echo "#######################################################################################"
echo "################### True negatives - should NOT see genmod add an annotation ##########"
echo "#######################################################################################"
run AD_TN.vcf
run AR_TN.vcf
run XD_TN.vcf
run XR_TN.vcf
run ARcomp_TN_ifRunningWiPhaseOpt.vcf "-phased"
run ARcomp_TN_ifRunningWiPhaseOpt.vcf # not calling -phase is more restrictive for TNs
<<COMMENT
