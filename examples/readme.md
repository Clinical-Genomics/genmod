# Examples for GENMOD #

These are some example files for getting to know **genmod** and the possible alternatives.
There are two small test families

- A "recessive" family in recessive_trio.ped, with two healthy parents and one affected child. In this case we would be interested in variants that follow the **autosomal recessive** or **autosomal recessive compound** inheritance pattern.
- A "dominant" family in dominant_trio.ped, with one affected parent, one healthy parent and one affected child.

There are annotation files included in the distribution of **genmod** that are made from the [latest](ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz) ref gene dataset. These will be used as default if no other annotation is given.
If the user want to build own annotations please use **genmod build_annotation**.

##Annotate variants for Recessive Family##


	genmod annotate test_data/test_vcf.vcf -fam test_data/recessive_trio.ped

The vcf file have a couple of variants made up so it will be easy to understand how the genetic inheritance patterns are annotated.

With the basic command listed above the output should look like

	##fileformat=VCFv4.1
	##INFO=<ID=ANN,Number=.,Type=String,Description="Annotates what feature(s) this variant belongs to.">
	##INFO=<ID=Comp,Number=.,Type=String,Description="':'-separated list of compound pairs for this variant.">
	##INFO=<ID=GM,Number=.,Type=String,Description="':'-separated list of genetic models for this variant.">
	##INFO=<ID=MS,Number=1,Type=Integer,Description="PHRED score for genotype models.">
	##contig=<ID=1,length=249250621,assembly=b37>
	##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband
	1	11900	.	C	A	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom;MS=55	GT:GQ	0/1:60	0/1:60	1/1:60
	1	11901	.	T	C	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom:AR_hom_dn;MS=57	GT:GQ	./.	0/1:60	1/1:60
	1	11902	.	T	G	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom_dn;MS=55	GT:GQ	0/1:60	0/0:60	1/1:60
	1	11903	.	C	G	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	0/1:60	1/1:60	1/1:60
	1	879585	.	A	T	100	PASS	MQ=1;ANN=NOC2L:SAMD11;Comp=1_879586_A_T;GM=AR_comp;MS=55	GT:GQ	0/1:60	0/0:60	0/1:60
	1	879586	.	A	T	100	PASS	MQ=1;ANN=NOC2L:SAMD11;Comp=1_879585_A_T;GM=AR_comp;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
	1	947378	.	T	G	100	PASS	MQ=1;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973348	.	G	A	100	PASS	MQ=1;ANN=AGRN;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973349	.	T	G	100	PASS	MQ=1;ANN=AGRN;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	3	879586	.	G	T	100	PASS	MQ=1;ANN=LOC101927215	GT:GQ	0/0:60	0/1:60	0/1:60
	3	947378	.	A	T	100	PASS	MQ=1;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	3	973348	.	G	A	100	PASS	MQ=1	GT:GQ	./.	0/1:60	0/1:60
	3	973349	.	A	T	100	PASS	MQ=1;GM=AD_dn:AD;MS=57	GT:GQ	./.	0/0:60	0/1:60

1. First variant is a classic Autosomal Recessive case, each parent are carriers and the affected child is homozygous alternative.
2. This variant is annotated as ``AR_hom:AR_hom_dn`` since we miss information from one parent.
3. Here we have a de novo case, notice that we will only annotate this with ``AR_hom_dn`` for more effective sorting in a downstream analysis.
4. Next two variants form a compound pair, that is each parent carry one variant each while the affected child have both. This implies that both variants are exonic(or in canonical splice region) and belong to the same gene.
5. At position 947378 we see an example of a Autosoma Dominant de novo variant.
6. The two following variants could form a compound pair since they are heterozygous in the affected child but they are not exonic. If the flag ``-g/--whole_gene`` is used these two will be annotated as compounds.

The following variants are to show how the ``-strict`` flag affects the analysis. When in strict mode we will only annotate a variant to follow a pattern if there is *proof* in the data. So if there are no calls in some individuals it will not follow any pattern. 


##Annotate variants for Dominant Family##

	genmod annotate test_data/test_vcf.vcf -fam test_data/dominant_trio.ped

Output:

	##fileformat=VCFv4.1
	##INFO=<ID=ANN,Number=.,Type=String,Description="Annotates what feature(s) this variant belongs to.">
	##INFO=<ID=Comp,Number=.,Type=String,Description="':'-separated list of compound pairs for this variant.">
	##INFO=<ID=GM,Number=.,Type=String,Description="':'-separated list of genetic models for this variant.">
	##INFO=<ID=MS,Number=1,Type=Integer,Description="PHRED score for genotype models.">
	##contig=<ID=1,length=249250621,assembly=b37>
	##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband
	1	11900	.	C	A	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	0/1:60	0/1:60	1/1:60
	1	11901	.	T	C	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	./.	0/1:60	1/1:60
	1	11902	.	T	G	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	0/1:60	0/0:60	1/1:60
	1	11903	.	C	G	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom;MS=55	GT:GQ	0/1:60	1/1:60	1/1:60
	1	879585	.	A	T	100	PASS	MQ=1;ANN=SAMD11:NOC2L	GT:GQ	0/1:60	0/0:60	0/1:60
	1	879586	.	A	T	100	PASS	MQ=1;ANN=SAMD11:NOC2L;GM=AD;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
	1	947378	.	T	G	100	PASS	MQ=1	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973348	.	G	A	100	PASS	MQ=1;ANN=AGRN	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973349	.	T	G	100	PASS	MQ=1;ANN=AGRN	GT:GQ	0/0:60	0/0:60	0/1:60
	3	879586	.	G	T	100	PASS	MQ=1;ANN=LOC101927215;GM=AD;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
	3	947378	.	A	T	100	PASS	MQ=1	GT:GQ	0/0:60	0/0:60	0/1:60
	3	973348	.	G	A	100	PASS	MQ=1;GM=AD;MS=57	GT:GQ	./.	0/1:60	0/1:60
	3	973349	.	A	T	100	PASS	MQ=1	GT:GQ	./.	0/0:60	0/1:60


We can now see how the conditions change when one of the parents are affected. For example the recessive pattern for the first variant is not followed since all affected needs to be homozygote alternative if the variant should follow the Autosomal Recessive pattern.

##Annotate variants with CADD scores and population frequencies##

This is another example of how one can annotate with genmod:

	genmod annotate examples/test_vcf.vcf --cadd_file examples/small_CADD.tsv.gz --thousand_g examples/small_1000G.vcf.gz

Output:

	##fileformat=VCFv4.1
	##INFO=<ID=Annotation,Number=.,Type=String,Description="Annotates what feature(s) this variant belongs to.">
	##INFO=<ID=Compounds,Number=.,Type=String,Description="':'-separated list of compound pairs for this variant.">
	##INFO=<ID=GeneModels,Number=.,Type=String,Description="':'-separated list of genetic models for this variant.">
	##INFO=<ID=ModelScore,Number=1,Type=Integer,Description="PHRED score for genotype models.">
	##INFO=<ID=CADD,Number=1,Type=Float,Description="The CADD relative score for this alternative.">
	##INFO=<ID=1000GMAF,Number=1,Type=Float,Description="Frequency in the 1000G database.">
	##contig=<ID=1,length=249250621,assembly=b37>
	##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
	##Software=<ID=genmod,Version=1.5.3,Date=2014-08-27 16:46:33.556358, CommandLineOptions="variant_file=examples/test_vcf.vcf thousand_g=examples/small_1000G.vcf.gz cadd_file=examples/small_CADD.tsv.gz outfile=res.tmp family_type=ped">
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband
	1	11900	.	C	A	100	PASS	MQ=1;Annotation=DDX11L1;CADD=0.106	GT:GQ	0/1:60	0/1:60	1/1:60
	1	11901	.	T	C	100	PASS	MQ=1;Annotation=DDX11L1;CADD=3.199	GT:GQ	./.	0/1:60	1/1:60
	1	11902	.	T	G	100	PASS	MQ=1;Annotation=DDX11L1;CADD=2.539	GT:GQ	0/1:60	0/0:60	1/1:60
	1	11903	.	C	G	100	PASS	MQ=1;Annotation=DDX11L1;CADD=1.459	GT:GQ	0/1:60	1/1:60	1/1:60
	1	879585	.	A	T	100	PASS	MQ=1;Annotation=SAMD11,NOC2L;CADD=5.268	GT:GQ	0/1:60	0/0:60	0/1:60
	1	879586	.	A	T	100	PASS	MQ=1;Annotation=SAMD11,NOC2L;CADD=4.221	GT:GQ	0/0:60	0/1:60	0/1:60
	1	947378	.	T	G	100	PASS	MQ=1;CADD=4.432;1000GMAF=0.01	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973348	.	G	A	100	PASS	MQ=1;Annotation=AGRN;CADD=5.587;1000GMAF=0.0018	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973349	.	T	G	100	PASS	MQ=1;Annotation=AGRN;CADD=3.073	GT:GQ	0/0:60	0/0:60	0/1:60
	3	879586	.	G	T	100	PASS	MQ=1;Annotation=LOC101927215;CADD=1.744;1000GMAF=0.0005	GT:GQ	0/0:60	0/1:60	0/1:60
	3	947378	.	A	T	100	PASS	MQ=1;CADD=2.345	GT:GQ	0/0:60	0/0:60	0/1:60
	3	973348	.	G	A	100	PASS	MQ=1;CADD=3.706	GT:GQ	./.	0/1:60	0/1:60
	3	973349	.	A	T	100	PASS	MQ=1;CADD=7.008	GT:GQ	./.	0/0:60	0/1:60
	


Please post issues on https://github.com/moonso/genmod if any problems.