# Examples for GENMOD #

These are some example files for getting to know **genmod** and the possible alternatives.
There are two small test families

- A "recessive" family in recessive_trio.ped, with two healthy parents and one affected child. In this case we would be interested in variants that follow the **autosomal recessive** or **autosomal recessive compound** inheritance pattern.
- A "dominant" family in dominant_trio.ped, with one affected parent, one healthy parent and one affected child.

There are annotation files included in the distribution of **genmod** that are made from the [latest](ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz) ref gene dataset. These will be used as default if no other annotation is given.
If the user want to build own annotations please use **genmod build_annotation**.

##Annotate variants for Recessive Family##


	genmod annotate test_data/recessive_trio.ped test_data/test_vcf.vcf

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
	1	11900	.	A	T	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom;MS=55	GT:GQ	0/1:60	0/1:60	1/1:60
	1	11901	.	C	T	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom:AR_hom_dn;MS=57	GT:GQ	./.	0/1:60	1/1:60
	1	11902	.	A	T	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom_dn;MS=55	GT:GQ	0/1:60	0/0:60	1/1:60
	1	11903	.	T	C	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	0/1:60	1/1:60	1/1:60
	1	879585	.	A	T	100	PASS	MQ=1;ANN=NOC2L:SAMD11;Comp=1_879586_A_T;GM=AR_comp;MS=55	GT:GQ	0/1:60	0/0:60	0/1:60
	1	879586	.	A	T	100	PASS	MQ=1;ANN=NOC2L:SAMD11;Comp=1_879585_A_T;GM=AR_comp;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
	1	947378	.	A	T	100	PASS	MQ=1;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973348	.	G	A	100	PASS	MQ=1;ANN=AGRN;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973349	.	T	G	100	PASS	MQ=1;ANN=AGRN;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	3	879586	.	A	T	100	PASS	MQ=1;ANN=LOC101927215	GT:GQ	0/0:60	0/1:60	0/1:60
	3	947378	.	A	T	100	PASS	MQ=1;GM=AD_dn;MS=55	GT:GQ	0/0:60	0/0:60	0/1:60
	3	973348	.	G	A	100	PASS	MQ=1	GT:GQ	./.	0/1:60	0/1:60
	3	973349	.	T	A	100	PASS	MQ=1;GM=AD_dn:AD;MS=57	GT:GQ	./.	0/0:60	0/1:60

1. First variant is a classic Autosomal Recessive case, each parent are carriers and the affected child is homozygous alternative.
2. This variant is annotated as ``AR_hom:AR_hom_dn`` since we miss information from one parent.
3. Here we have a de novo case, notice that we will only annotate this with ``AR_hom_dn`` for more effective sorting in a downstream analysis.
4. Next two variants form a compound pair, that is each parent carry one variant each while the affected child have both. This implies that both variants are exonic(or in canonical splice region) and belong to the same gene.
5. At position 947378 we see an example of a Autosoma Dominant de novo variant.
6. The two following variants could form a compound pair since they are heterozygous in the affected child but they are not exonic. If the flag ``-g/--whole_gene`` is used these two will be annotated as compounds.

The following variants are to show how the ``-strict`` flag affects the analysis. When in strict mode we will only annotate a variant to follow a pattern if there is *proof* in the data. So if there are no calls in some individuals it will not follow any pattern. 


##Annotate variants for Dominant Family##

	genmod annotate test_data/dominant_trio.ped test_data/test_vcf.vcf

Output:

	##fileformat=VCFv4.1
	##INFO=<ID=ANN,Number=.,Type=String,Description="Annotates what feature(s) this variant belongs to.">
	##INFO=<ID=Comp,Number=.,Type=String,Description="':'-separated list of compound pairs for this variant.">
	##INFO=<ID=GM,Number=.,Type=String,Description="':'-separated list of genetic models for this variant.">
	##INFO=<ID=MS,Number=1,Type=Integer,Description="PHRED score for genotype models.">
	##contig=<ID=1,length=249250621,assembly=b37>
	##reference=file:///humgen/gsa-hpprojects/GATK/bundle/current/b37/human_g1k_v37.fasta
	#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	father	mother	proband
	1	11900	.	A	T	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	0/1:60	0/1:60	1/1:60
	1	11901	.	C	T	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	./.	0/1:60	1/1:60
	1	11902	.	A	T	100	PASS	MQ=1;ANN=DDX11L1	GT:GQ	0/1:60	0/0:60	1/1:60
	1	11903	.	T	C	100	PASS	MQ=1;ANN=DDX11L1;GM=AR_hom;MS=55	GT:GQ	0/1:60	1/1:60	1/1:60
	1	879585	.	A	T	100	PASS	MQ=1;ANN=SAMD11:NOC2L	GT:GQ	0/1:60	0/0:60	0/1:60
	1	879586	.	A	T	100	PASS	MQ=1;ANN=SAMD11:NOC2L;GM=AD;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
	1	947378	.	A	T	100	PASS	MQ=1	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973348	.	G	A	100	PASS	MQ=1;ANN=AGRN	GT:GQ	0/0:60	0/0:60	0/1:60
	1	973349	.	T	G	100	PASS	MQ=1;ANN=AGRN	GT:GQ	0/0:60	0/0:60	0/1:60
	3	879586	.	A	T	100	PASS	MQ=1;ANN=LOC101927215;GM=AD;MS=55	GT:GQ	0/0:60	0/1:60	0/1:60
	3	947378	.	A	T	100	PASS	MQ=1	GT:GQ	0/0:60	0/0:60	0/1:60
	3	973348	.	G	A	100	PASS	MQ=1;GM=AD;MS=57	GT:GQ	./.	0/1:60	0/1:60
	3	973349	.	T	A	100	PASS	MQ=1	GT:GQ	./.	0/0:60	0/1:60


We can now see how the conditions change when one of the parents are affected. For example the recessive pattern for the first variant is not followed since all affected needs to be homozygote alternative if the variant should follow the Autosomal Recessive pattern.



Please post issues on https://github.com/moonso/genmod if any problems.