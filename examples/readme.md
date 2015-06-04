# Examples for GENMOD #

These are some example files for getting to know **genmod** and the possible alternatives.
There are three test families

- A "recessive" family in recessive_trio.ped, with two healthy parents and one affected child. In this case we would be interested in variants that follow the **autosomal recessive** or **autosomal recessive compound** inheritance pattern.
- A "dominant" family in dominant_trio.ped, with one affected parent, one healthy parent and one affected child.
- A multi family file with two families.

There are annotation files included in the distribution of **genmod** that are made from the [latest](ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz) refGene dataset. These will be used as default if no other annotation is given.
If the user want to build own annotations please use **genmod build_annotation**.

##Annotate variants for Recessive Family##


	genmod annotate examples/test_vcf.vcf -f examples/recessive_trio.ped -o examples/test_vcf_recessive_annotated.vcf

The vcf file have a couple of variants made up so it will be easy to understand how the genetic inheritance patterns are annotated.

With the basic command listed above the output should look like the variants in ```examples/test_vcf_recessive_annotated.vcf```

1. First variant is a classic Autosomal Recessive case, each parent are carriers and the affected child is homozygous alternative.
2. This variant is annotated as ``AR_hom:AR_hom_dn`` since we miss information from one parent.
3. Here we have a de novo case, notice that we will only annotate this with ``AR_hom_dn`` for more effective sorting in a downstream analysis.
4. Next two variants form a compound pair, that is each parent carry one variant each while the affected child have both. This implies that both variants are exonic(or in canonical splice region) and belong to the same gene.
5. At position 947378 we see an example of a Autosoma Dominant de novo variant.
6. The two following variants could form a compound pair since they are heterozygous in the affected child but they are not exonic. If the flag ``-g/--whole_gene`` is used these two will be annotated as compounds.

The following variants are to show how the ``-strict`` flag affects the analysis. When in strict mode we will only annotate a variant to follow a pattern if there is *proof* in the data. So if there are no calls in some individuals it will not follow any pattern. 


##Annotate variants for Dominant Family##

	genmod annotate test_data/test_vcf.vcf -f test_data/dominant_trio.ped -o examples/test_vcf_dominant_annotated.vcf

We can now see how the conditions change when one of the parents are affected. For example the recessive pattern for the first variant is not followed since all affected needs to be homozygote alternative if the variant should follow the Autosomal Recessive pattern.


##Annotate variants for Multiple Families##

	genmod annotate test_data/test_vcf.vcf -f test_data/multi_family.ped -o examples/test_vcf_multi_annotated.vcf

We can now see how the conditions change when one of the parents are affected. For example the recessive pattern for the first variant is not followed since all affected needs to be homozygote alternative if the variant should follow the Autosomal Recessive pattern.


##Annotate variants with CADD scores and population frequencies##

This is another example of how one can annotate with genmod:

	genmod annotate examples/test_vcf.vcf --cadd_file examples/small_CADD.tsv.gz --thousand_g examples/small_1000G.vcf.gz


Please post issues on http://github.com/moonso/genmod if any problems.