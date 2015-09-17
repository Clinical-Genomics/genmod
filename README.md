# GENMOD #



[![DOI](https://zenodo.org/badge/5735/moonso/genmod.png)](http://dx.doi.org/10.5281/zenodo.11424)


**GENMOD** is a simple to use command line tool for annotating and analyzing genomic variations in the [VCF](http://samtools.github.io/hts-specs/VCFv4.1.pdf) file format.
**GENMOD** can annotate genetic patterns of inheritance in vcf:s with single or multiple families of arbitrary size.

The tools in the genmod suite are:

- **genmod annotate**, for annotating regions, frequencies, cadd scores etc.
- **genmod build_annotation**, for building new annotation sets from different sources
- **genmod models**, For annotating patterns of inheritance
- **genmod sort**, To sort the variants of a vcf file

##Installation:##

**GENMOD** works with Python 2.7 and Python v3.2 and above

    pip install genmod

or

	git clone https://github.com/moonso/genmod.git
	cd genmod
	python setup.py install



# WARNING: Deprecated after version 2.6. This will be updated soon!! Please see command line help for now #

##USAGE:##

###genmod annotate###


    genmod annotate variant_file.vcf --family_file ped_file

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

