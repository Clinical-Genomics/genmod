# GENMOD #



[![DOI](https://zenodo.org/badge/5735/moonso/genmod.png)](http://dx.doi.org/10.5281/zenodo.11424)
[![Build Status](https://travis-ci.org/moonso/vcf_parser.svg)](https://travis-ci.org/moonso/genmod)


**GENMOD** is a simple to use command line tool for annotating and analyzing genomic variations in the [VCF](http://samtools.github.io/hts-specs/VCFv4.1.pdf) file format.
**GENMOD** can annotate genetic patterns of inheritance in vcf:s with single or multiple families of arbitrary size.

The tools in the genmod suite are:

- **genmod build**, for building new annotation sets from different sources
- **genmod annotate**, for annotating regions, frequencies, cadd scores etc.
- **genmod models**, For annotating patterns of inheritance
- **genmod sort**, To sort the variants of a vcf file, either on rank score or position
- **genmod score**, Score the variants of a vcf based on their annotation
- **genmod filter**, Filter the variants of a vcf based on their annotation

##Installation:##

**GENMOD** works with Python 2.7 and Python 3.2 and above

    pip install genmod

or

	git clone https://github.com/moonso/genmod.git
	cd genmod
	python setup.py install


## USAGE: ##

<!-- TODO change documentation link -->
*This is an overview, for more in depth documentation see [documentation](http://moonso.github.io/genmod/)*


### Example: ###


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

<!-- ###genmod annotate###


    genmod annotate variant_file.vcf

This will print a new vcf to standard out with all variants annotated according to the statements below.
All individuals described in the [ped](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped) file must be present in the vcf file

See examples in the folder ```genmod/examples```.

**From version 1.9 genmod can split multiallelic calls in vcf:s, use flag -split/--split_variants.**

To get an example of how splitting variants work, run genmod on the file ```examples/multi_allele_example.vcf``` with the dominant trio.
That is:
    ```genmod annotate examples/multi_allele_example.vcf -f examples/dominant_trio.ped -split```

Compare the result when not using the ```-split``` flag.

Genmod is distributed with a annotation database that is built from the refGene data.
If the user wants to build a new annotation set use the command below:

	genmod build_annotation [--type] annotation_file


Each variant in the VCF-file will be annotated with which genetic models that are followed in the family if a family file
(ped file) is provided.

The genetic models that are checked are the following:

* Autsomal Recessive, denoted 'AR_hom'
* Autsomal Recessive denovo, denoted 'AR\_hom\_dn'
* Autsomal Dominant, 'AD'
* Autsomal Dominant denovo, 'AD_dn'
* Autosomal Compound Heterozygote, 'AR_comp'
* X-linked dominant, 'XD'
* X-linked dominant de novo, 'XD_dn'
* X-linked Recessive, 'XR'
* X-linked Recessive de novo, 'XR_dn'

Se description of how genetic models are annotated in the section **Conditions for genetic models** below.

It is possible to run without a family file, in this case all variants will be annotated with which region(s) they belong to, and if other annotation files are provided(1000G, CADD scores etc.) the variants will get the proper values from these.

[Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html)(vep) annotations are supported, use the ```--vep```-flag if variants are already annotated with vep.

**GENMOD** will add entrys to the INFO column for the given VCF file depending on what information is given.

If ```--vep``` is NOT provided:

- **Annotation** Comma separated list with features overlapped in the annotation file

If ```--vep``` is used **Annotation** will not be annotated since all information is in the vep entry.

If a pedigree file is provided the following will be added:

- **GeneticModels** A comma separated list with which genetic models that are followed in each family described in the ped file. Annotation are separated with pipes on the form ```GeneticModels=fam_id_1:AR_hom, fam_id_2:AR_comp|AD_dn``` etc..
- **Compounds** Comma separated list with compound pairs(if any) for each family. These are described like 'CHR\_POS\_REF\_ALT'
- **ModelScore** Model Score, a phred-score based on the genotype qualities to describe the uncertainty of the genetic model in each family

Also a line for logging is added in the vcf header with the id **genmod**, here the date of run, version and command line arguments are printed.

- Compound heterozygote inheritance pattern will be checked if two variants are exonic (or in canonical splice sites) and if they reside in the same gene.

- If compounds should be checked in the whole gene (including introns) use ```--whole_gene```
- GENMOD supports phased data, use the ```-phased``` flag. Data should follow the [GATK way](http://gatkforums.broadinstitute.org/discussion/45/read-backed-phasing) of phasing.

All annotations will be present only if they have a value.

- **GENMOD** can annotate the variants with 1000 genome frequencies. Use the flag `-kg/--thousand_g path/to/bgzipped/thousand_genomes.vcf.gz`
- **GENMOD** also supports annotation of frequencies from the [ExAC](http://exac.broadinstitute.org/). Use the flag `--exac path/to/bgzipped/ExAC_file.vcf.gz`
- Annotate with [CADD scores](http://cadd.gs.washington.edu/), use `-cadd/--cadd_file path/to/huge_cadd_file.tsv.gz`.
- There several cadd files with different variant sets to cover as much as possible.
	- One with all 1000 genomes positions (this one include some indels), if annotation with this one use `-c1kg/--cadd_1000_g path/to/CADD_1000g.txt.gz`.
	- One with all variants from the ESP6500 dataset. If annotation with this one use `--cadd_esp path/to/CADD_ESP.tsv.gz`.
	- One with all variants from the ExAC dataset. If annotation with this one use `--cadd_exac path/to/CADD_ExAC.tsv.gz`.
	- One with 12.3M InDels from the CADD resources. If annotation with this one use `--cadd_indels path/to/CADD_InDels.txt.gz`.
- By default the relative cadd scores is annotated with 'CADD=score', there is also an alternative to annotate with the raw cadd scores using the `--cadd_raw` flag. In this case a info field 'CADD_raw=score'.
- If your VCF is already annotated with VEP, use `-vep/--vep`
- If data is phased use `-phased/--phased`
- If you want to allow compound pairs in intronic regions to use `-gene/--whole_gene`
- If you want canonical splice site region to be bigger than 2 base pairs on each side of the exons, use `-splice/--splice_padding <integer>`
- The `-strict/--strict` flag tells **genmod** to only annotate genetic models if they are proved by the data. If a variant is not called in a family member it will not be annotated.


###genmod build_annotation###

	genmod build_annotation [--type] [-o/--outdir] annotation_file

The following file formats are supported for building new annotations:

- bed
- ccds
- gtf
- gene_pred

The user can also specify the amount of positions around exon boundaries that should be considered as splice sites. Use

```--splice_padding INTEGER```

###genmod analyze###

From version 1.6 there is also a tool for analyzing the variants annotated by **genmod**. This tool will look at all variants in a vcf and do an analysis based on which inheritance patterns they follow. The variants are then ranked based on the cadd scores, the highest ranked variants for each category is printed to screen and the full list for each category is printed to new vcf files.
Run with:

	genmod analyze path/to/file.vcf

For more information do

	genmod analyze --help


### genmod sort ###


Sort a VCF file based on Rank Score.

```
Usage: genmod sort [OPTIONS] <vcf_file> or -

  Sort a VCF file based on rank score.

Options:
  -o, --outfile PATH    Specify the path to a file where results should be
                        stored.
  -f, --family_id TEXT  Specify the family id for sorting. If no family id the
                        first family found in annotation will be used.
  -v, --verbose         Increase output verbosity.
  --help                Show this message and exit.
```

###genmod summarize###

Tool to get basic statistics of the annotated in a vcf file.
Run

	genmod summarize --help

for more information.

## Conditions for Genetic Models ##

### Short explanation of genotype calls in VCF format:###

Since we only look at humans, that are diploid, the genotypes represent what we see on both alleles in a single position.
0 represents the reference sequence, 1 is the first of the alternative alleles, 2 second alternative and so on.
If no phasing has been done the genotype is an unordered pair on the form x/x, so 0/1 means that the individual is heterozygote in this given position with the reference base on one of the alleles and the first of the alternatives on the other.
2/2 means that we see the second of the alternatives on both alleles.
Some chromosomes are only present in one copy in humans, here it is allowed to only use a single digit to show the genotype. A 0 would mean reference and 1 first of alternatives.

If phasing has been done the pairs are not unordered anymore and the delimiter is then changed to '|', so one can be heterozygote in two ways; 0|1 or 1|0.


### Autosomal Recessive ###

For this model individuals can be carriers so healthy individuals can be heterozygous. Both alleles need to have the variant for an individual to be sick so a healthy individual can not be homozygous alternative and a sick individual *has* to be homozygous alternative.

* Affected individuals have to be homozygous alternative (hom. alt.)
* Healthy individuals cannot be hom. alt.
* Variant is considered _de novo_ if both parents are genotyped and do not carry the variant


### Autosomal Dominant ###

* Affected individuals have to be heterozygous (het.)
* Healthy individuals cannot have the alternative variant
* Variant is considered _de novo_ if both parents are genotyped and do not carry the variant


### Autosomal Compound Heterozygote ###

This model includes pairs of exonic variants that are present within the same gene.
**The default behaviour of GENMOD is to look for compounds only in exonic/canonical splice sites**.
The reason for this is that since some genes have huge intronic regions the data will be drowned in compound pairs.
If the user wants all variants in genes checked use the flag -gene/--whole_gene.

1. Non-phased data:
	* Affected individuals have to be het. for both variants
	* Healthy individuals can be het. for one of the variants but cannot have both variants
	* Variant is considered _de novo_ if only one or no variant is found in the parents


2. Phased data:
	* All affected individuals have to be het. for both variants **and** the variants has to be on two different alleles
	* Healthy individuals can be heterozygous for one but cannot have both variants
	* If only one or no variant is found in parents it is considered _de novo_


### X-Linked Dominant###

These traits are inherited on the x-chromosome, of which men have one allele and women have two.

* Variant has to be on chromosome X
* Affected individuals have to be het. or hom. alt.
* Healthy males cannot carry the variant
* Healthy females can carry the variant (because of X inactivation)
* If sex is male the variant is considered _de novo_ if mother is genotyped and does not carry the variant
* If sex is female variant is considered _de novo_ if none of the parents carry the variant


### X Linked Recessive ###

* Variant has to be on chromosome X
* Affected males have to be het. or hom. alt. (het is theoretically not possible in males, but can occur due to Pseudo Autosomal Regions).
* Affected females have to be hom. alt.
* Healthy females cannot be hom. alt.
* Healthy males cannot carry the variant
* If sex is male the variant is considered _de novo_ if mother is genotyped and does not carry the variant
* If sex is female variant is considered _de novo_ if not both parents carry the variant
 -->
