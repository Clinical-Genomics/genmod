# Annotating Patterns of Inheritance #


### Basic command ###

```bash

$ genmod models --help
Usage: genmod models [OPTIONS] <vcf_file> or -

  Annotate genetic models for vcf variants.

  Checks what patterns of inheritance that are followed in a VCF file. The
  analysis is family based so each family that are specified in the family
  file and exists in the variant file will get it's own annotation.

Options:
 Options:
   -f, --family_file <ped_file>
   -t, --family_type [ped|alt|cmms|mip]
                                   If the analysis use one of the known setups,
                                   please specify which one.
   -r, --reduced_penetrance <tsv_file>
                                   File with gene ids that have reduced
                                   penetrance.
   --vep                           If variants are annotated with the Variant
                                   Effect Predictor.
   --phased                        If data is phased use this flag.
   -s, --strict                    If strict model annotations should be
                                   used(see documentation).
   -p, --processes INTEGER         Define how many processes that should be use
                                   for annotation.
   --silent                        Do not print the variants.
   -w, --whole_gene                If compounds should be checked over the
                                   whole gene.
   -k, --keyword TEXT              What annotation keyword that should be used
                                   when
                                   searching for features.
   -o, --outfile FILENAME          Specify the path to a file where results
                                   should be stored.
   --help                          Show this message and exit.

```

## Overview ##

The ``genmod models`` command is used to annotate patterns of inheritance in vcf files. 
**genmod** can handle multiple families of arbitrary sizes in the same vcf.  
The individuals and families to be included in the analysis is specified in 
the family file that is given with the ``-f/--family_file`` option.

## Options and Arguments ##

### vcf_file ###

This is the only mandatory argument. A proper vcf file or a stream with a vcf file is allways input

### family_file ###

This is a file in ped like format that describes which individuals that should be considered in the analysis and how they are related.
The default file format is [ped](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped) but there are some alternatives.
Please see ```-t/--family_type```.

### family_type ###

Family files are usually in the ped format but there are possibilities to use ped like files.
If a ped file includes more columns than the first six mandatory (described [here](http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped)) one may use the ```alt``` option for this parameter.

### reduced_penetrance ###

Some phenotypes are known to have reduced penetrance. This means that healthy 
individuals can carry a heterozygous variant without being affected, while others
that have the same variant are affected.
When annotating inheritance patterns with **genmod** the user can provide a tab separated 
file where the first column describes a gene id per row.
Use the flag ```-r/--reduced_penetrance``` to provide this file to genmod.
Genmod will then allow healthy carriers for variants that resides in these genes when
 annotating the Autosomal Dominant inheritance pattern.

### vep ###

This flag tell genmod that the variants are annotated with the [Variant Effect Predictor](http://www.ensembl.org/info/docs/tools/vep/index.html).
In this case the vep annotations will be used to determine if a compound pair is in the same gene.

### strict ###

This flag tell genmod that the variant calls are phased. This will make it possible for more accurate compound checks.
For all autosomal chromosomes humans have two copies. If two variants resides on the same copy there will still be one copy that should produce "healthy" transcripts. 

With phased variants **genmod** can take this information into account.

### processes ###

How many processes should be used during the analysis. If ```-p/--processes 1``` the streaming functionality will work better for genmod since no intermediate file needs to be used, the variants will be processes line by line and printed directly.

If there are several processes running at the same time the variants will be printed unordered to a intermediate file that will be sorted when all variants are analyzed.

### silent ###
 
If no variants or headers should be printed to screen. This is mainly for testing.

### whole_gene ###
 
If variants in introns should be considered when checking compound pairs. This will lead to a huge number of compound candidates/pairs when dealing with whole genome data.

### keyword ###

This is the keyword in the vcf info field that genmod will look for when trying to find what genes the variant is annotated with.
Default is ```Annotation``` since this is the keyword used by ```genmod annotate```.



