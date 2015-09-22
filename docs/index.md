# Genmod #

<p align="center">
  <a href="https://github.com/moonso/genmod">
    <img src="https://github.com/moonso/genmod/raw/master/artwork/tree_man.JPG"/>
  </a>
</p>

Annotating and analyzing primarily genetic models in the [VCF](http://samtools.github.io/hts-specs/VCFv4.1.pdf) file format. 

## Overview ##


The main function is to annotate **accurate patterns of inheritance** for each 
variant in a vcf file including one or several families of arbitrary size. 

By default variants will be annotated from the human refGene database but 
alternative annotations can be built with ``genmod build``. 
**genmod** can also be used for annotating variants with their frequencies in 1000G, 
ExAC or their predicted CADD scores, see **genmod annotate**.
genmod uses multithreading to annotate variants and operates fast for both 
whole exome data and whole genome data.

## Installation ##

**GENMOD** works with Python 2.7 and Python 3.2 and above

    pip install genmod

or

	git clone https://github.com/moonso/genmod.git
	cd genmod
	python setup.py install



### Example: ###

The basic idea with genmod is to make fast and easy analysis of vcf variants 
for rare disease.
It can still be interesting to use in other cases, such as annotating what 
genetic regions the variants in a bacteria belongs to.
**genmod** can annotate accurate patterns of inheritance in arbitrary sized families.
The genetic models checked are the basic mendelian ones, these are:

* Autsomal Recessive, denoted 'AR_hom'
* Autsomal Recessive denovo, denoted 'AR\_hom\_dn'
* Autsomal Dominant, 'AD'
* Autsomal Dominant denovo, 'AD_dn'
* Autosomal Compound Heterozygote, 'AR_comp'
* X-linked dominant, 'XD'
* X-linked dominant de novo, 'XD_dn'
* X-linked Recessive, 'XR'
* X-linked Recessive de novo, 'XR_dn'

The following command should work when installed successfully. The files are distributed with the package.

```bash
$ cat examples/test_vcf.vcf
##fileformat=VCFv4.1
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##contig=<ID=1,length=249250621,assembly=b37>
##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband	father_2	mother_2	proband_2
1	879537	.	T	C	100	PASS	MQ=1	GT:AD:GQ	0/1:10,10:60	0/1:10,10:60	1/1:10,10:60	0/0:10,10:60	0/1:10,10:60	1/1:10,10:60
1	879541	.	G	A	100	PASS	MQ=1	GT:AD:GQ	./.	0/1:10,10:60	1/1:10,10:60	./.	0/1:10,10:60	0/1:10,10:60
1	879595	.	C	T	100	PASS	MQ=1	GT:AD:GQ	0/1:10,10:60	0/0:10,10:60	1/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60
1	879676	.	G	A	100	PASS	MQ=1	GT:AD:GQ	0/1:10,10:60	1/1:10,10:60	1/1:10,10:60	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60
1	879911	.	G	A	100	PASS	MQ=1	GT:AD:GQ	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60
1	880012	.	A	G	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60
1	880086	.	T	C	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
1	880199	.	G	A	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
1	880217	.	T	G	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
10	76154051	.	A	G	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60
10	76154073	.	T	G	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
10	76154074	.	C	G	100	PASS	MQ=1	GT:AD:GQ	./.	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60
10	76154076	.	G	C	100	PASS	MQ=1	GT:AD:GQ	./.	0/0:10,10:60	0/1:10,10:60	./.	0/0:10,10:60	0/1:10,10:60
X	302253	.	CCCTCCTGCCCCT	C	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	1/1:10,10:60	0/0:10,10:60	1/1:10,10:60	1/1:10,10:60
MT	302253	.	CCCTCCTGCCCCT	C	100	PASS	MQ=1	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	1/1:10,10:60	0/0:10,10:60	1/1:10,10:60	1/1:10,10:60

$ cat examples/test_vcf.vcf |\
>genmod annotate - --annotate_regions |\ 
>genmod models - --family_file examples/recessive_trio.ped > test_vcf_models_annotated.vcf

$ cat test_vcf_models_annotated.vcf
##fileformat=VCFv4.1
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS Mapping Quality">
##INFO=<ID=Annotation,Number=.,Type=String,Description="Annotates what feature(s) this variant belongs to.">
##INFO=<ID=Exonic,Number=0,Type=Flag,Description="Indicates if the variant is exonic.">
##INFO=<ID=GeneticModels,Number=.,Type=String,Description="':'-separated list of genetic models for this variant.">
##INFO=<ID=ModelScore,Number=.,Type=String,Description="PHRED score for genotype models.">
##INFO=<ID=Compounds,Number=.,Type=String,Description="List of compound pairs for this variant.The list is splitted on ',' family id is separated with compoundswith ':'. Compounds are separated with '|'.">
##contig=<ID=1,length=249250621,assembly=b37>
##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
##Software=<ID=genmod,Version=3.0.1,Date="2015-09-22 08:40",CommandLineOptions="processes=4 keyword=Annotation family_type=ped family_file=<open file 'examples/recessive_trio.ped', mode 'r' at 0x102d3a780> variant_file=<_io.TextIOWrapper name='<stdin>' encoding='utf-8'> logger=<logging.Logger object at 0x102d64250>">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband	father_2	mother_2	proband_2
1	879537	.	T	C	100	PASS	MQ=1;Exonic;Annotation=SAMD11;GeneticModels=1:AR_hom;ModelScore=1:55.0	GT:AD:GQ	0/1:10,10:60	0/1:10,10:60	1/1:10,10:60	0/0:10,10:60	0/1:10,10:60	1/1:10,10:60
1	879541	.	G	A	100	PASS	MQ=1;Exonic;Annotation=SAMD11;GeneticModels=1:AR_hom_dn|AR_hom;ModelScore=1:57.0	GT:AD:GQ	./.	0/1:10,10:60	1/1:10,10:60	./.	0/1:10,10:60	0/1:10,10:60
1	879595	.	C	T	100	PASS	MQ=1;Exonic;Annotation=NOC2L,SAMD11;GeneticModels=1:AR_hom_dn;ModelScore=1:55.0	GT:AD:GQ	0/1:10,10:60	0/0:10,10:60	1/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60
1	879676	.	G	A	100	PASS	MQ=1;Exonic;Annotation=NOC2L,SAMD11	GT:AD:GQ	0/1:10,10:60	1/1:10,10:60	1/1:10,10:60	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60
1	879911	.	G	A	100	PASS	MQ=1;Exonic;Annotation=NOC2L,SAMD11;Compounds=1:1_880086_T_C|1_880012_A_G;GeneticModels=1:AR_comp|AR_comp_dn;ModelScore=1:55.0	GT:AD:GQ	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60
1	880012	.	A	G	100	PASS	MQ=1;Exonic;Annotation=NOC2L;Compounds=1:1_879911_G_A|1_880086_T_C;GeneticModels=1:AR_comp|AR_comp_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60
1	880086	.	T	C	100	PASS	MQ=1;Exonic;Annotation=NOC2L;Compounds=1:1_879911_G_A|1_880012_A_G;GeneticModels=1:AD_dn|AR_comp_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
1	880199	.	G	A	100	PASS	MQ=1;Annotation=NOC2L;GeneticModels=1:AD_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
1	880217	.	T	G	100	PASS	MQ=1;Annotation=NOC2L;GeneticModels=1:AD_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
10	76154051	.	A	G	100	PASS	MQ=1;Exonic;Annotation=ADK;Compounds=1:10_76154073_T_G;GeneticModels=1:AR_comp_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60	0/0:10,10:60	0/1:10,10:60	0/1:10,10:60
10	76154073	.	T	G	100	PASS	MQ=1;Exonic;Annotation=ADK;Compounds=1:10_76154051_A_G;GeneticModels=1:AD_dn|AR_comp_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60	0/0:10,10:60	0/0:10,10:60	0/1:10,10:60
10	76154074	.	C	G	100	PASS	MQ=1;Annotation=ADK	GT:AD:GQ	./.	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60	0/1:10,10:60
10	76154076	.	G	C	100	PASS	MQ=1;Annotation=ADK;GeneticModels=1:AD_dn|AD;ModelScore=1:57.0	GT:AD:GQ	./.	0/0:10,10:60	0/1:10,10:60	./.	0/0:10,10:60	0/1:10,10:60
X	302253	.	CCCTCCTGCCCCT	C	100	PASS	MQ=1;Annotation=PPP2R3B;GeneticModels=1:XD|XR;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	1/1:10,10:60	0/0:10,10:60	1/1:10,10:60	1/1:10,10:60
MT	302253	.	CCCTCCTGCCCCT	C	100	PASS	MQ=1;GeneticModels=1:AR_hom_dn;ModelScore=1:55.0	GT:AD:GQ	0/0:10,10:60	0/1:10,10:60	1/1:10,10:60	0/0:10,10:60	1/1:10,10:60	1/1:10,10:60
```

**genmod** is made for working on any type of annotated vcf.
To get relevant Autosomal Compound Heterozygotes we need to know what genetic regions that the variants belong to.
We can use annotations from the [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html) 
or let **genmod** do the annotation.

**genmod** comes with a prebuilt annotation data base that is made from the latest refSeq dataset. 
We can also build new annotation sets with **genmod build**, please see wiki for mor info.

(There are files for testing the following commands in genmod/examples)

To annotate the variants with regions use

```bash
$genmod annotate <vcf_file> -r/--annotate_regions (-a/--annotation_dir)

```

Now the variants are ready to get their models annotated:

```bash
$genmod models <vcf_file> -f/--family_file <family.ped>

```