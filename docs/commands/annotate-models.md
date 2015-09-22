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
  -f, --family_file <ped_file>
  -t, --family_type [ped|alt|cmms|mip]
                                  If the analysis use one of the known setups,
                                  please specify which one.
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
 