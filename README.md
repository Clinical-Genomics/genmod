# GENMOD #



[![DOI](https://zenodo.org/badge/5735/moonso/genmod.png)](http://dx.doi.org/10.5281/zenodo.11424)


##USAGE:##

###Basic functions###

    genmod annotate variant_file.vcf --family_file ped_file

This will print a new vcf to standard out with all variants annotated according to the statements below.

Genmod is distributed with a annotation database that are built from the refGene data.
If the user wants to build a new annotation set use the command below:

	genmod build_annotation [--type] annotation_file


##General##

Tool for annotating patterns of inheritance in Variant Call Format (VCF) files with arbitrary pedigrees.

Each variant in the VCF-file will be annotated with which genetic models that are followed in the family if a family file
(ped file) is provided.
It is possible to run without a family file, in this case all variants will be annotated with which region(s) they belong to, and if other annotation files are provided(1000G, CADD scores etc.) the variants will get the proper values from these.

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

**GENMOD** will add entrys to the INFO column for the given VCF file depending on what information is given.

If ```-vep/--vep``` is NOT provided:
- **Annotation** Comma separated list with features overlapped in the annotation file

If a pedigree file is provided:

- **GeneticModels** A comma separated list with genetic models followed
- **Compounds** Comma separated list with compound pairs(if any). These are described like 'CHR\_POS\_REF\_ALT'.
- **ModelScore** Model Score, a phred-score based on the genotype qualities to describe the uncertainty of the genetic model.

Also a line for logging is added in the vcf header with the id **genmod**, here the date of run, version and command line arguments are printed.

All annotations will be present only if they have a value.

If ```-vep/--vep``` is used **Annotation** will not be annotated since all information is in the vep entry.

##Installation:##

**GENMOD** works with Python 2.7 and Python v3.2 and above

    pip install genmod

or

    git clone https://github.com/moonso/genmod.git
    cd genmod
    python setup.py install

###Alternatives###

- **GENMOD** can annotate the variants with 1000 genome frequencies. Use the flag `-kg/--thousand_g path/to/bgzipped/thousand_genomes.vcf.gz`
- Annotate with [CADD scores](http://cadd.gs.washington.edu/), use `-cadd/--cadd\_file path/to/huge_cadd_file.txt.gz`. 
- There several cadd files with different variant sets to cover as much as possible. 
	- One with all 1000 genomes positions (this one include some indels), if annotation with this one use `-c1kg/--cadd_1000_g path/to/CADD_1000g.txt.gz`. 
	- One with all variants from the ESP6500 dataset. If annotation with this one use `--cadd_esp path/to/CADD_ESP.txt.gz`.
	- One with 12.3M InDels from the CADD resources. If annotation with this one use `--cadd_indels path/to/CADD_InDels.txt.gz`.
- By default the relative cadd scores is annotated with 'CADD=score', there is also an alternative to annotate with the raw cadd scores using the `--cadd_raw` flag. In this case a info field 'CADD_raw=score'.
- If your VCF is already annotated with VEP, use `-vep/--vep`
- If data is phased use `-phased/--phased`
- If you want to allow compound pairs in intronic regions to use `-gene/--whole_gene`
- If you want canonical splice site region to be bigger than 2 base pairs on each side of the exons, use `-splice/--splice_padding _integer_`
- The `-strict/--strict` flag tells **genmod** to only annotate genetic models if they are proved by the data. If a variant is not called in a family member it will not be annotated.

###Distribution###

- GENMOD includes db like files in the genmod/annotations folder, this is the exon and gene definitions from ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz.

If the user wants to use another annotation:

    genmod build_annotation [--type] "annotation type" path/to/annotation_file -o path/to/outdir

In this case the new annotation will be built into the outdir specified (default is the genmods annotation dir).
When the user want to annotate a vcf with this new annotation set use 

    genmod annotate variant_file -fam ped_file [-a/--annotation_dir] /path/to/new/annotation_dir


- Compound heterozygote inheritance pattern will be checked if two variants are exonic (or in canonical splice sites) and if they reside in the same gene.

- GENMOD supports phased data use the -phased flag. Data should follow the [GATK way](http://gatkforums.broadinstitute.org/discussion/45/read-backed-phasing) of phasing.

- GENMOD support VCF files annotated with [VEP](http://www.ensembl.org/info/docs/tools/vep/index.html), use -vep flag. This means that GENMOD will use the **VEP** annotation for checking if variants are in the same gene. 

- GENMOD can annotate variants with their [CADD](http://cadd.gs.washington.edu/) score. This is done by adding the flag -cadd "path/to/cadd_file".


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




##Details##

###Exonic variants###

Variants are defined as exonic if they are within an interval that is defined as an exon in the annotation file, and if they are in the canonical splice sites. The size of the canonical splice sites can be altered with -splice 'integer', in this case the annotation needs to be rebuilded. Example:

    run_genmod ped_file variant_file -an refGene.txt -at gene_pred -splice 6

In this case the exons will be padded by 6 bases on each side.

If VEP annotation is used the following SO-terms is counted as compound candidate sites:

**transcript ablation, splice donor variant, splice acceptor variant, stop gained, frameshift variant, stop lost,
initiator codon variant, inframe insertion, inframe deletion, missense variant, transcript amplification, splice region variant,incomplete terminal codon variant, synonymous variant, stop retained variant, coding sequence variant**.



###Annotation Formats###

The following formats are supported:

- [gene pred](http://genome.ucsc.edu/FAQ/FAQformat.html#format9). This is the format of the refSeq genes
- [gtf](http://www.ensembl.org/info/website/upload/gff.html). This is the format use by ensembl
- CCDS format
- [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1)format.

When BED format is used all entrys will both count as exons and genes(for compounds).


<!-- ## Detailed Structure ##

Here all attributes and methods of the classes will be showed:

### Genotype ###

Store the genotype information of a variant that is specific for an individual

**Attributes:**

* genotype STRING
* allele_1 STRING
* allele_2 STRING
* nocall BOOL
* heterozygote BOOL
* homo_alt BOOL
* homo_ref BOOL
* has_variant BOOL
* filter STRING
* ref_depth INT
* alt_depth INT
* phred_likelihoods TUPLE with INT
* depth_of_coverage INT
* genotype_quality FLOAT


### Variant ###

Holds the info of a variant and it's specific behaviour in a family.

**Attributes**

* chr STRING Have to be string since X, Y
* start INT
* stop INT 
* ref STRING Reference nucleotide(s)
* alt STRING Alternative sequence
* identity STRING dbSNP-id
* var_info DICT A dictionary with all the info from the variant file
* qual STRING A value for the score of the base call
* filter STRING The filter status
* genotypes LIST A list with the genotypes found for this variants
* gene STRING Semicolon separated string with ensemble gene names
* ad BOOL If following Autosomal Dominant pattern
* ad_dn BOOL If following Autosomal Dominant De novo pattern
* ar BOOL If following Autosomal Recessive pattern
* ar_dn BOOL If following Autosomal Recessive De nove pattern
* ar_comp BOOL If following Autosomal Recessive compound pattern
* ar_comp_dn BOOL If following Autosomal Recessive Compound De Novo pattern

**Methods**

* get_variant(self):
    Returns a dictionary with basic info to stdout
* print_model_info(self):
    Print for each variant which patterns of inheritance they follow.    
* print_vcf_variant(self):
    Print the variant in vcf-format to stdout
* print_original_version(self, header_columns):
    Prints the variant in its original format.
* check_noncomplete_call(self):
    Check if GATK have missed to report some info.
* get_genotype(self):
    Returns the list with genotypes for this variant.


### Individual ###

Holds the information about an individual and the individual specific genotypes.

**Attributes**

* ind STRING Can be any id unique within the family
* family STRING Can be any unique id within the cohort
* mother STRING The ind_id of the mother or [0,-9] if info is missing
* father STRING ---------||------ father --------------||---------------
* sex INT 1=male 2=female 0=unknown
* phenotype INT 1=unaffected, 2=affected, missing = [0,-9]
* genotypes DICT Container with genotype information on the form {<variant_id>: <Genotype>}
* phasing BOOL If the genotype information includes phasing for this individual

### Family ###

**Attributes**

* individuals DICT dictionary with family members on the form {<ind_id>:<Individual>}
* variants DICT dictionary with all the variants that exists in the family on the form {<var_id>:<Variant>} -->


<!-- [![Bitdeli Badge](https://d2weczhvl823v0.cloudfront.net/moonso/genmod/trend.png)](https://bitdeli.com/free "Bitdeli Badge") -->

