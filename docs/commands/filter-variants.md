# Filter Variants #

Filter variants based on some annotation from the INFO field of a vcf file.

### Examples ###

```bash
genmod -v filter examples/test_vcf_annotated.vcf -a 1000G_freq -t 0.01
```

This command will skip all variants that have a ```1000G_freq``` < 0.01.
By default all variants that lack ```1000G_freq``` will be kept, this can be modified by the ```-d/--discard``` flag.

```bash
genmod -v filter examples/test_vcf_annotated.vcf -a CADD -t 7 --greater --discard
```

Here we will filter out all variants that have annotation ```CADD``` less than 7 and all variants that miss the ```CADD``` flag.

## Command ##

```bash
$ genmod filter --help
Usage: genmod filter [OPTIONS] <vcf_file> or -

  Filter vcf variants.

  Filter variants based on their annotation

Options:
  -a, --annotation TEXT   Specify the info annotation to search for. Default
                          1000GAF
  -t, --threshold FLOAT   Threshold for filter variants. Default 0.05
  -d, --discard           If variants without the annotation should be
                          discarded
  -g, --greater           If greater than threshold should be used instead of
                          less thatn threshold.
  -s, --silent            Do not print the variants.
  -o, --outfile FILENAME  Specify the path to a file where results should be
                          stored.
  --help                  Show this message and exit.
```

