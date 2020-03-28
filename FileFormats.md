## Input files for NPS
To run NPS, you need the following set of input files: 

1. **GWAS summary statistics.** NPS supports two summary statistics formats: *MINIMAL* and *PREFORMATTED*. Internally, NPS uses PREFORMATTED summary statistics. With real data, however, we generally recommend preparing summary statistics in the MINIMAL format and harmonize them with training genotype data using provided NPS scripts. This will automatically convert the summary statistics file into the PREFORMATTED format ([step-by-step instruction](https://github.com/sgchun/nps#how-to-prepare-training-and-validation-cohorts-for-nps)). 
   - The MINIMAL format is a *tab-delimited* text file with the following seven or eight columns: 
     - **chr**: chromosome number. NPS expects only chromosomes 1-22. *Only chromosome numbers are expected.*
     - **pos**: base position of SNP.
     - **a1** and **a2**: alleles at each SNP in any order.
     - **effal**: effect allele. It should match either a1 or a2 allele. 
     - **pval**: p-value of association. 
     - **effbeta**: estimated *per-allele* effect size of *the effect allele*. For case/control GWAS, log(OR) should be used. *DO NOT pre-convert them to effect sizes relative to standardized genotypes. NPS will handle this automatically.*
     - **effaf**: (*Optional*) allele frequency of *the effect allele* in the discovery GWAS cohort. If this column is missing, NPS will use the allele frequencies of training cohort instead. Although this is optional, we **strongly recommend including effaf when it is available from a GWAS study.** When this column is provided, NPS will run a QC check for the consistency of the effect allele frequencies  between GWAS and training cohort data. 
     ```
     chr	pos	a1	a2	effal	pval	effbeta	effaf
     1	569406	G	A	G	0.8494	0.05191	0.99858
     1	751756	C	T	C	0.6996	0.00546	0.14418
     1	753405	C	A	C	0.8189	0.00316	0.17332
     1	753541	A	G	A	0.8945	0.00184	0.16054
     1	754182	A	G	A	0.7920	0.00361	0.18067
     1	754192	A	G	A	0.7853	0.00373	0.1809
     1	754334	T	C	T	0.7179	0.00500	0.18554
     1	755890	A	T	A	0.7516	0.00441	0.17327
     1	756604	A	G	A	0.9064	0.00162	0.18202
     ...
     ```
     
   - The *PREFORMATTED* format is the native format for NPS. We provide summary statistics of [our test cases](https://github.com/sgchun/nps#test-cases) in this format so that NPS can run on them directly without conversion. This is a *tab-delimited* text file format, and rows must be sorted by chromosome numbers and positions of SNPs. The following seven columns are required: 
     - **chr**: chromosome name starting with "chr." NPS expects only chromosomes 1-22. *Chromosomes should be designated by "chr1", ..., "chr22".*
     - **pos**: base position of SNP.
     - **ref** and **alt**: reference and alternative alleles of SNP, respectively.
     - **reffreq**: allele frequency of reference allele in the discovery GWAS cohort. 
     - **pval**: p-value of association. 
     - **effalt**: estimated *per-allele* effect size of *the alternative allele*. For case/control GWAS, log(OR) should be used. NPS will convert **effalt** to effect sizes relative to *the standardized genotype* internally using **reffreq**.  
     ```
     chr	pos	ref	alt	reffreq	pval	effalt
     chr1	676118	G	A	0.91584	0.7908	0.0012
     chr1	734349	G	A	0.90222	0.6989	0.001636
     chr1	770886	G	A	0.91708	0.721	0.001627
     chr1	785050	G	A	0.1139	0.3353	-0.00381
     chr1	798400	G	A	0.8032	0.03301	0.006736
     chr1	804759	G	A	0.8837	0.7324	-0.00134
     chr1	831489	G	A	0.2797	0.1287	0.004252
     chr1	832318	G	A	0.2797	0.4102	0.002304
     chr1	836924	G	A	0.7958	0.6591	-0.001374
     ...
     ```
     When summary statistics are provided in this format, NPS will not run its automatic data harmonization procedures. The user needs to ensure that the summary statistics file does not include InDels, tri-allelic SNPs or duplicated markers and also that no SNP in training genotypes files has missing summary statistics. If these requirements are violated, NPS will report an error and terminate. 

2. **Genotypes in the qctool dosage format.** NPS expects individual genotype data are provided in the dosage file format. We use [qctool](https://www.well.ox.ac.uk/~gav/qctool/) to generate these files (See [instructions](https://github.com/sgchun/nps#how-to-prepare-training-and-validation-cohorts-for-nps)). The genotype files need to be split by chromosomes for parallelization, and for each chromosome, the file should be named as "chrom*N*.*DatasetID*.dosage.gz." 
   
3. **Sample IDs in the PLINK .fam format.** The samples in the .fam file should appear in the exactly same order as the samples in the  genotype dosage files. This is a space- or tab-separated six-column text file without a column header. The phenotype information in this file will be ignored. See [PLINK documentation](https://www.cog-genomics.org/plink2/formats#fam) for the details on the format. 
   ```
   trainF2  trainI2  0  0  0 -9
   trainF3  trainI3  0  0  0 -9
   trainF39 trainI39 0  0  0 -9
   trainF41 trainI41 0  0  0 -9
   trainF58 trainI58 0  0  0 -9
   ```
4. **Phenotypes in the PLINK phenotype format.** NPS looks up phenotypes in a separately prepared phenotype file. The phenotype name has to be "**Outcome**" with cases and controls encoded by **1** and **0**, respectively. **FID** and **IID** are used together to match samples to .fam file. This file is tab-delimited, and samples can appear in any order. Missing phenotypes (e.g. missing entry of samples in .fam file or phenotypes encoded by **-9**) are not allowed.
   ```
   FID   IID    Outcome
   trainF2  trainI2  0
   trainF39 trainI39 0
   trainF3  trainI3  1
   trainF41 trainI41 1
   trainF58 trainI58 0
   ```

