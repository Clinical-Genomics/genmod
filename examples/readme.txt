This is a small test family with a couple of variants.

run with:

python scripts/run_genmod.py test_data/small_pedigree.ped test_data/test_vcf.vcf test_data/small_geneset.gtf -at gtf

the vcf file have 4 variants

1. Should follow autosomal recessive
2 and 3. Should be a compound pair
4 This should follow adtosomal dominant de novo

Output:

##fileformat=VCFv4.1
##INFO=<ID=ANN,Number=.,Type=String,Description="Annotates what feature(s) this variant belongs to.">
##INFO=<ID=Comp,Number=.,Type=String,Description="':'-separated list of compound pairs for this variant.">
##INFO=<ID=GM,Number=.,Type=String,Description="':'-separated list of genetic models for this variant.">
##INFO=<ID=MS,Number=1,Type=Integer,Description="PHRED score for genotype models.">
##contig=<ID=1,length=249250621,assembly=b37>
##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband
1	11900	.	A	T	100	PASS	MQ=1;ANN=ENSG00000223972;Comp=-;GM=AR_hom;MS=55	GT:GQ	0/1:60	0/1:60	1/1:60
1	879585	.	A	T	100	PASS	MQ=1;ANN=ENSG00000187634:ENSG00000188976;Comp=1_879586_A_T;GM=AR_comp;MS=55	GT:GQ	0/1:60	0/0:60	0/1:60
1	879586	.	A	T	100	PASS	MQ=1;ANN=ENSG00000187634:ENSG00000188976;Comp=1_879585_A_T;GM=AR_comp;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
1	947378	.	A	T	100	PASS	MQ=1;ANN=ENSG00000224969;Comp=-;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60



Please start an issue on https://github.com/moonso/genmod if any problems.