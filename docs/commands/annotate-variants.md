# Annotate Variants #

## Overview ##


To annotate the compound heterozygote model of inheritance we need to know if variants reside in the same gene.
With ``genmod annotate <file.vcf> -r/--annotate_regions`` each variant will be annotate with what genes they overlap and if they are exonic variants. By default the annotations that follow the package is used, these are built from the latest refSeq dataset.
It is possible to build a new annotation database with the ``genmod build`` command if the user prefere other annotations.   

``genmod annotate`` can also be used to annotate variants from different sources such as frequency databases etc.

### Command ###

```bash
$ genmod annotate --help
Usage: genmod annotate [OPTIONS] <vcf_file> or -

  Annotate vcf variants.

  Annotate variants with a number of different sources. Please use --help
  for more info.

Options:
  -r, --annotate_regions     Increase output verbosity.
  -c, --cadd_file PATH       Specify the path to a bgzipped cadd file (with
                             index) with variant scores. This command can be
                             used multiple times if multiple cadd files.
  --thousand_g PATH          Specify the path to a bgzipped vcf file (with
                             index) with 1000g variants
  --exac PATH                Specify the path to a bgzipped vcf file (with
                             index) with exac variants.
  -a, --annotation_dir PATH  Specify the path to the directory where the
                             annotation
                             databases are.
                             Default is the gene
                             pred files that comes with the distribution.
  -o, --outfile FILENAME     Specify the path to a file where results should
                             be stored.
  -s, --silent               Do not print the variants.
  --cadd_raw                 If the raw cadd scores should be annotated.
  -p, --processes INTEGER    Define how many processes that should be use for
                             annotation.
  --help                     Show this message and exit.
```