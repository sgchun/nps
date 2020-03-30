## Input files for NPS
To run NPS, you need the following set of input files: 

1. **GWAS summary statistics.** This is a *tab-delimited* text file format, and rows must be sorted by chromosome numbers and positions of SNPs. The following seven columns are required: 
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
     The user needs to ensure that the summary statistics file does not include InDels, tri-allelic SNPs or duplicated markers and also that no SNP in training genotypes files has missing summary statistics. If these requirements are violated, NPS will report an error and terminate. 

2. **Genotype files in the qctool dosage format.** NPS expects individual genotype data are provided in the dosage file format. We use [qctool](https://www.well.ox.ac.uk/~gav/qctool/) to generate these files (See [instructions](https://github.com/sgchun/nps#how-to-prepare-training-and-validation-cohorts-for-nps)). The genotype files need to be split by chromosomes for parallelization, and for each chromosome, the file should be named as "chrom*N*.*DatasetID*.dosage.gz." 
   
3. **Sample IDs in the PLINK .fam format.** The samples in the .fam file should appear in the exactly same order as the samples in the  genotype dosage files. This is a space- or tab-separated six-column text file without a column header. The phenotype information in this file will be ignored. See [PLINK documentation](https://www.cog-genomics.org/plink2/formats#fam) for the details on the format. 
   ```
   trainF2  trainI2  0  0  1 -9
   trainF3  trainI3  0  0  2 -9
   trainF39 trainI39 0  0  1 -9
   trainF41 trainI41 0  0  2 -9
   trainF58 trainI58 0  0  1 -9
   ```
4. **Phenotypes in the PLINK phenotype format.** NPS looks up phenotypes in a separately prepared phenotype file. The phenotype name has to be "**Outcome**" with cases and controls encoded by "2" and **1**, respectively. **FID** and **IID** are used together to match samples to .fam file. This file is tab-delimited, and samples can appear in any order. Missing phenotypes (e.g. missing entry of samples in .fam file or phenotypes encoded by **-9**) are not allowed.
   ```
   FID   IID    Outcome
   trainF2  trainI2  1
   trainF39 trainI39 1
   trainF3  trainI3  2
   trainF41 trainI41 2
   trainF58 trainI58 1
   ```

