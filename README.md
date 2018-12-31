﻿
# Non-Parametric Shrinkage (NPS)
NPS is a non-parametric polygenic risk prediction algorithm described in Chun et al. (2018) BioRxiv [(preprint)](https://www.biorxiv.org/content/early/2018/07/16/370064). 

## How to install

1. The core NPS module is implemented in R. R can be downloaded from [here](https://www.r-project.org/). R-3.0 or later is required to run NPS. Although NPS can run on a plain-vanilla version of R, we strongly recommend to use R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [R open](https://mran.microsoft.com/open). These libraries reduce the running time of NPS substantially by speeding up matrix manipulation steps in NPS.  

2. (Optional) NPS relies on R modules, `pROC` and `DescTools`, to calculate the AUC and Nagelkerke's *R2* statistics. These modules are optional; if they are not installed, AUC and Nagelkerke's *R2* will not be reported. To enable this feature, please install these packages by running the following on command line: 

```bash
$ Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
$ Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
```

In case that it is preferred to install these R extensions in your home directory (e.g. ~/R) instead of the default system path, please do the following instead:
```bash
$ Rscript -e 'install.packages("pROC", "~/R", repos="http://cran.r-project.org")' 
$ Rscript -e 'install.packages("DescTools", "~/R", repos="http://cran.r-project.org")' 

# Add "~/R" to your local R library path in the login shell start up file 
# For examaple in case of bash, add this: 
$ export R_LIBS=~/R:$R_LIBS
```

3. Download and unpack NPS package as below. Some of NPS codes are optimized in C++ and need to be compiled. For this, you need GNU C++ compiler with C++0x support (GCC 4.4 or later). This step will create two executable binaries, `stdgt` and `grs`, in the top-level NPS directory. `stdgt` converts an imputed dosage file to standardized genotypes with the mean of 0 and variance of 1. `grs` calculates genetic risk scores for all indiviudals in an imputed dosage file using per-SNP genetic effects computed by NPS.

```bash
$ tar -zxvf nps-1.0.0.tar.gz
$ cd nps-1.0.0/
$ make
```

Note on cluster use: If you loaded a GCC module to compile NPS binaries, you also need to load the GCC module in `nps_stdgt.job` and `nps_score.job` in order to use shared C++ libraries in run time. 

4. We recommend to run NPS on computer clusters, processing all chromosomes in parallel. To make this easier, we provide job scripts for SGE and LSF clusters. Please see `sge` and `lsf` directories along with provided examples [below](https://github.com/sgchun/nps#test-cases). You may need to modify the provided job scripts to load necessary modules if they are not loaded by default. For example, you may need to do the following lines in the job scripts (Do not blindly add these lines. The details depend on individual system configurations): 

```bash 
###
# ADD CODES TO LOAD R MODULE HERE
#
# On clusters running environment modules and providing R-mkl
module add gcc/5.3.0 
module add R-mkl/3.3.2

# On clusters running DotKit instead and supporting OpenblasR
use GCC-5.2 
use OpenblasR
```

5. (Optional) We provide job scripts to prepare UK Biobank data for NPS training and validation cohorts. To gain access to UK Biobank data, please see [UK Biobank data access application procedure](https://www.ukbiobank.ac.uk/). Our scripts rely on [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) and [QCTOOL v2](https://www.well.ox.ac.uk/~gav/qctool/). If you use a different cohort for training dataset, you can still use these scripts as far as the genotype data are provided in .bgen files [**LINK TO INSTRUCTIONS HERE**]. 

## Input files for NPS
To run NPS, you need the following set of input files: 

1. **GWAS summary statistics.** Currently, NPS supports two summary statistics formats: *preformated* and *minimal*. 
   - The *preformatted* summary statistics format is the native format immediately ready for NPS. We provide summary statistics for [our test cases](https://github.com/sgchun/nps#test-cases) in this format. This is a tab-delimited text file, sorted by chromsome numbers and positions and with the following seven essential columns: 
     - `chr`: chromosome name starting with "chr." NPS expects only chromosomes chr1-chr22.
     - `pos`: base positions of SNP.
     - `ref` and `alt`: reference and alternative alleles of SNP. NPS does not allow InDels or tri-allelic SNPs. There should not be duplicated SNPs in the file. 
     - `reffreq`: allele frequency of reference allele in the discovery GWAS cohort. 
     - `pval`: p-value of association. 
     - `effalt`: estimated *per-allele* effect size of alternative allele. For case/control GWAS, log(OR) should be used. NPS will convert `effalt` to effect sizes relative to the standardized genotype using `reffreq` values.  
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
     ```
   - Because sometimes it is cumbersome to convert publicly available summary statistics into the preformatted summary statistics, we support the *minimal* summary statistics format, which can be automatically converted into the preformatted format and harmonized with training genotype data (See [**LINK HERE**]). Currently, this feature is supported only when UK Biobank is used as a training cohort. The minimal format is a tab-delimited text file with the seven or eight columns: 
     - `chr`: chromosome number. NPS expects only chromosomes 1-22.
     - `pos`: base positions of SNP.
     - `a1` and `a2`: Alleles at each SNP in any order. There should not be duplicated SNPs in the file. 
     - `effal`: the effect allele. It should be either a1 or a2 allele. 
     - `pval`: p-value of association. 
     - `effbeta`: estimated *per-allele* effect size of the effect allele. For case/control GWAS, log(OR) should be used. 
     - `effaf`: (Optional) allele frequency of effect allele in the discovery GWAS cohort. If this column is provided, the allele frequency will be compared between GWAS and training cohort, and markers with too divergent allele frequencies will be filtered out. If this column is not provided, allele frequencies of training cohort will be copied to summary statistics. 
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
     ```

2. **Training genotypes in QCTOOL dosage format.** Genotype data have to be prepared in the dosage format. NPS requires that all markers in this file overlap with GWAS summary statitics. We recommend to remove InDels, tri-allelic SNPs, rare variants with MAF < 5%, markers with any QC issue, and markers that are not found in the GWAS summary statistics. Markers with very different allele frequencies between GWAS and training cohort should be also discarded. NPS does not allow duplicated SNPs in the dosage file. Genotype data needs to be split by chromosomes for parallelization, and each file should be named as "chrom*N*.*CohortName*.dosage.gz." If you use UK Biobank data, all these QC processing steps can be done automatically using `ukbb_support` scripts [**LINK HERE**]. 

NPS matches SNPs using the combination of chromosome, position and alleles on the assumed forward strand (+) and does not rely on `SNPID` or `rsid`. NPS expects `alleleA` to match `ref` allele and `alleleB` to match `alt` allele in the GWAS summary statistics. The allelic dosage counts the genetic dosage of `alleleB`. Markers has to be sorted by `position`. [QCTOOL](https://www.well.ox.ac.uk/~gav/qctool/) can be used to generate the dosage files.  
```
chromosome SNPID rsid position alleleA alleleB trainI2 trainI3 trainI39 trainI41 trainI58
01 1_676118:676118:G:A 1_676118:676118:G:A 676118 G A 0 0 0 0 0
01 1_734349:734349:G:A 1_734349:734349:G:A 734349 G A 1 0 0 0 0
01 1_770886:770886:G:A 1_770886:770886:G:A 770886 G A 0 0 0 1 0
01 1_785050:785050:G:A 1_785050:785050:G:A 785050 G A 1 2 2 2 2
01 1_798400:798400:G:A 1_798400:798400:G:A 798400 G A 1 0 1 1 0
01 1_804759:804759:G:A 1_804759:804759:G:A 804759 G A 1 0 1 0 0
01 1_831489:831489:G:A 1_831489:831489:G:A 831489 G A 1 2 2 1 2
01 1_832318:832318:G:A 1_832318:832318:G:A 832318 G A 1 2 2 1 2
01 1_836924:836924:G:A 1_836924:836924:G:A 836924 G A 0 0 0 0 0
...
```
3. **Training sample IDs in PLINK .fam format.** The samples in the .fam file should appear in the exactly same order as the samples in the training genotype files. This is *space-separated* six-column text file *without a header*. Only first two columns (family ID and individual ID) will be used. See [here](https://www.cog-genomics.org/plink2/formats#fam) for additional details on the format. 
```
trainF2 trainI2 0 0 0 -9
trainF3 trainI3 0 0 0 -9
trainF39 trainI39 0 0 0 -9
trainF41 trainI41 0 0 0 -9
trainF58 trainI58 0 0 0 -9
...
```

4. **Training phenotypes in PLINK phenotype format.** NPS looks up phenotypes in the separate phenotype file and ignores the phenotype data provided in the .fam file above. The phenotype name has to be `Outcome` with cases and controls encoded with `1` and `0`, respectively. When the optional `TotalLiability` column is provided, NPS also reports *R2* on the liability scale, namely *R2* between polygenic scores and `TotalLiablity`. The combination of `FID` and `IID` are used to match samples to .fam file. Samples can appear in any order in this file. Missing phenotypes (e.g. encoded with `-9` or missing entry of samples described in .fam file) are not allowed. 
```
FID	IID	TotalLiability	Outcome
trainF68266	trainI68266	2.00037013482937	1
trainF77481	trainI77481	2.19146817328717	1
trainF39184	trainI39184	0.833227631267031	0
trainF76746	trainI76746	0.168893429255002	0
trainF65453	trainI65453	0.817881328612934	0
...
```
5. **Validation genotypes in QCTOOL dosage format.** This is for validation cohort. See the above training genotypes for the details. 
6. **Validation sample IDs in PLINK .fam format.** This is for validation cohort. See the above training sample IDs for the details. 
7. **Validation phenotypes in PLINK phenotype format.** This is for validation cohort. See the above training phenotypes for the details. Missing phenotypes are allowed in this file (encoded with `-9`). Such samples will be simply excluded when calculating the accuracy statistics. 


## Test cases
We provide two sets of simulated test cases. Due to their large file sizes, they are provided separately from the software distribution. Please download them from Sunyaev Lab (ftp://genetics.bwh.harvard.edu/download/schun/). Test set #1 is a relatively small dataset (225MB), and NPS can complete in less than 1 hour in total on a modest desktop PC without linear-algebra acceleration. In contrast, test set #2 is more realistic simulation (11GB) and will require serious computational resource: NPS will generate up to 1 TB of intermediary data and take 1/2 day to a day to run on computer clusters (the exact overall running time will depend on the degree of parallelization).  

Both simulation datasets were generated using our multivariate-normal simulator. See our NPS manuscript for the details.   
- **Test set #1.** The number of markers across the genome is limited to 100,449. We assume that all causal SNPs are included in the 100,449 SNPs. Note that this is unrealistic assumption; causal SNPs cannot be accurately tagged with such a sparse SNP set. The fraction of causal SNP is 0.005 (a total of 522 SNPs). The GWAS cohort size is 100,000. The training cohort has 2,500 cases and 2,500 controls. The valiation cohort is 5,000 samples without case over-sampling. The heritability is 0.5. The phenotype prevalence is 5%.  

- **Test set #2.** The number of markers across the genome is 5,012,500. The fraction of causal SNP is 0.001 (a total of 5,008 SNPs). The GWAS cohort size is 100,000. The training cohort has 2,500 cases and 2,500 controls. The valiation cohort is 5,000 samples without case over-sampling. The heritability is 0.5. The phenotype prevalence is 5%. This is one of the benchmark simulation datasets used in our manuscript, with 1/10-fold reduced validation cohort size. 

We assume that the test datasets will be downloaded and unpacked in the following directories: 
```bash
$ cd nps-1.0.0/testdata/
$ tar -zxvf NPS.Test1.tar.gz 
# This will create the following test data files in nps-1.0.0/testdata/Test1
# Test1/Test1.summstats.txt (GWAS summary statistics)
# Test1/chrom1.Test1.train.dosage.gz (training cohort genotypes)
# Test1/chrom2.Test1.train.dosage.gz (training cohort genotypes)
# ... 
# Test1/Test1.train.2.5K_2.5K.fam (training cohort sample IDs)
# Test1/Test1.train.2.5K_2.5K.phen (training cohort phenotypes)
# Test1/chrom1.Test1.val.dosage.gz (validation cohort genotypes)
# Test1/chrom2.Test1.val.dosage.gz (validation cohort genotypes)
# ... 
# Test1/Test1.val.5K.fam (validation cohort sample IDs)
# Test1/Test1.val.5K.phen (validation cohort phenotypes)

$ tar -zxvf NPS.Test2.tar.gz 
# This will create the following test data in nps-1.0.0/testdata/Test2
# Test2/Test2.summstats.txt (GWAS summary statistics)
# Test2/chrom1.Test2.train.dosage.gz (training cohort genotypes)
# Test2/chrom2.Test2.train.dosage.gz (training cohort genotypes)
# ... 
# Test2/Test2.train.2.5K_2.5K.fam (training cohort sample IDs)
# Test2/Test2.train.2.5K_2.5K.phen (training cohort phenotypes)
# Test2/chrom1.Test2.val.dosage.gz (validation cohort genotypes)
# Test2/chrom2.Test2.val.dosage.gz (validation cohort genotypes)
# ... 
# Test2/Test2.val.5K.fam (validation cohort sample IDs)
# Test2/Test2.val.5K.phen (validation cohort phenotypes)
```

### Running NPS on Test set #1

1. **Standardize genotypes.** The first step is to standardize the training genotypes to the mean of 0 and variance of 1. `sge/snps_stdgt.job` and `lsf/snps_stdgt.job` scripts run the computation in parallel across all 22 chromosomes. On a desktop computer, `batch_all_chroms.sh` can be used to run an SGE job script for all chromosomes one by one in batch. The first parameter (`testdata/Test1`) is the location of training cohort (.dosage.gz files), the second parameter (`Test1.train`) is the *CohortName* of training cohort, and the last parameter (`5000`) is the number of samples in the training cohort. 
```bash
$ cd nps-1.0.0/

# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_stdgt.job testdata/Test1 Test1.train 5000

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test1 Test1.train 5000
```
After all jobs are completed, `nps_check.sh` script can be used to make sure that all jobs are successful. If failure is detected, `FAIL` message will be printed. Note `nps_check.sh` script can be used on clusters (both SGE and LSF) and in batch mode as follows in the same way: 
```bash
$ ./nps_check.sh stdgt testdata/Test1 Test1.train 
Verifying nps_stdgt:
Checking testdata/Test1/chrom1.Test2.train ...OK
Checking testdata/Test1/chrom2.Test2.train ...OK
Checking testdata/Test1/chrom3.Test2.train ...OK
...
```

2. **Configure an NPS run.** NPS will create the directory to store intermediate data (`testdata/Test1/npsdat`) and save the NPS configuration file in the directory. For parameters, NPS needs the path to GWAS summary statistics file (`testdata/Test1/Test1.summstats.txt`), directory containing training genotypes (`testdata/Test1`), IDs of training samples (`testdata/Test1/Test1.train.2.5K_2.5K.fam`), phenotypes of training samples (`testdata/Test1/Test1.train.2.5K_2.5K.phen`), name of training cohort (`Test1.train`), analysis window size (`80` SNPs here), and directory for intermediate data (`testdata/Test1/npsdat`). The window size of 80 SNPs for ~100,000 genome-wide SNPs is comparable to 4,000 SNPs for ~5,000,000 genome-wide SNPs. 

```bash
# Same on clusters and for batch processing (no parallelization)
$ Rscript npsR/nps_init.R testdata/Test1/Test1.summstats.txt testdata/Test1 testdata/Test1/Test1.train.2.5K_2.5K.fam testdata/Test1/Test1.train.2.5K_2.5K.phen Test1.train 80 testdata/Test1/npsdat

# Check the results
$ ./nps_check.sh init testdata/Test1/npsdat/
Checking testdata/Test1/npsdat//args.RDS ...OK (version 1.0.0)
Checking testdata/Test1/npsdat//log ...OK
```

3. **Transform data to the decorrelated "eigenlocus" space.** In general, this is one of the most time-consuming steps in NPS. We recommend to run NPS four times on shifted windows and merge the results in the later steps. For the window shift, we recommend the shifts of 0, 1/WINSZ, 2/WINSZ and 3/WINSZ SNPs, where WINSZ is the size of analysis window. In the case of test set #1, they correspond to 0, 20, 40 and 60. The first parameter is the location of intermediary data (`testdata/Test1/npsdat/`), and the second parameter is the window shift (`0`, `20`, `40` or `60`). 
```bash
# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 0
$ ./batch_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 20
$ ./batch_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 40
$ ./batch_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 0 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 20 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 40 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 60 

# Check the results
$ ./nps_check.sh decor testdata/Test1/npsdat/ 0 
$ ./nps_check.sh decor testdata/Test1/npsdat/ 20 
$ ./nps_check.sh decor testdata/Test1/npsdat/ 40 
$ ./nps_check.sh decor testdata/Test1/npsdat/ 60 
```

4. **Prune correlations across windows.** This step prunes the correlation between genotypes across adjacent windows in the eigenlocus space.  
```bash
# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 0
$ ./batch_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 20
$ ./batch_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 40
$ ./batch_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh prune testdata/Test1/npsdat/ 0 
$ ./nps_check.sh prune testdata/Test1/npsdat/ 20 
$ ./nps_check.sh prune testdata/Test1/npsdat/ 40 
$ ./nps_check.sh prune testdata/Test1/npsdat/ 60 
```

5. **Separate GWAS-significant partition.** The partition of GWAS-significant associations will be separated out from the rest of association signals. NPS takes longer time to complete this step when there are more GWAS-significant signals.
```bash
# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 0
$ ./batch_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 20
$ ./batch_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 40
$ ./batch_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 0 
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 20 
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 40 
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 60 
```

6. **Define a partitioning scheme.** The partition scheme will be defined with `npsR/nps_prep_part.R` and then, partitioned genetic risk scores will be calculated for all training samples using `nps_part.job`. Specifically, for `npsR/nps_prep_part.R`, the first parameter is the location of intermediary data (`testdata/Test1/npsdat/`), the second is the window shift (`0`, `20`, `40` or `60`), the third is the number of partitions on intervals of eigenvalues of eigenlocus projection (`10`), and the last is the number of partitions on intervals of observed effect sizes in the eigenlocus space (`10`). For `nps_part.job`, the first parameter is the location of intermediary data (`testdata/Test1/npsdat/`), and the second is the window shift (`0`, `20`, `40` or `60`)
```
# Same on clusters and for batch processing (no parallelization)
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 0 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 20 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 40 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 60 10 10 

# Check the results
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 0 
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 20 
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 40 
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 60 

# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 0
$ ./batch_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 20
$ ./batch_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 40
$ ./batch_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh part testdata/Test1/npsdat/ 0 
$ ./nps_check.sh part testdata/Test1/npsdat/ 20 
$ ./nps_check.sh part testdata/Test1/npsdat/ 40 
$ ./nps_check.sh part testdata/Test1/npsdat/ 60 
```

7. **Estimate per-partition shrinkage weights.** Then, we estimate the per-partition shrinkage weights using `npsR/nps_weight.R`. We also provide two optional utilities: `npsR/nps_train_AUC.R`, which reports the AUC statistics of prediction in training cohort, and `npsR/nps_plot_shrinkage.R`, which plots the overall curve of GWAS effect sizes re-weighted by per-partition shrinkage. Both tools take the average of NPS run on shifted windows. `npsR/nps_plot_shrinkage.R` will save the shrinkage curve plot in the pdf file path given as second argument (`Test1.nps.pdf`). 

```bash
# Same on clusters and for batch processing (no parallelization)
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 0 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 20 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 40 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 60 

# Check the results
$ ./nps_check.sh weight testdata/Test1/npsdat/ 0 
$ ./nps_check.sh weight testdata/Test1/npsdat/ 20 
$ ./nps_check.sh weight testdata/Test1/npsdat/ 40 
$ ./nps_check.sh weight testdata/Test1/npsdat/ 60 
```

```bash
# Optional 
$ Rscript npsR/nps_train_AUC.R testdata/Test1/npsdat/ 0 20 40 60
```

With test set #1, `npsR/nps_train_AUC.R` will report the following AUC in the training cohort: 
```
Data: 2500 controls < 2500 cases.
Area under the curve: 0.8799
95% CI: 0.8707-0.8891 (DeLong)
```

`npsR/nps_train_AUC.R` will store the following plot of shrinkage curve. See our example [(Test1.nps.pdf)](https://github.com/sgchun/nps/blob/master/testdata/Test1.nps.pdf).
```
# Same on clusters and for batch processing (no parallelization)
$ Rscript npsR/nps_plot_shrinkage.R testdata/Test1/npsdat/ Test1.nps.pdf 0 20 40 60
```

8. **Convert back to per-SNP effect sizes.** Then, the re-weighted effect sizes should be converted back to the original per-SNP space from the eigenlocus space. This will generate Test1.train.adjbetahat.chrom*N*.txt for the shift of 0, and Test1.train.win_*shift*.adjbetahat.chrom*N*.txt for the rest of shifts. 
```bash
# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 0
$ ./batch_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 20
$ ./batch_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 40
$ ./batch_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 0 
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 20 
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 40 
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 60 
```

9. **Validate.**

```bash
# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ Test1.train testdata/ Test1.val
$ ./batch_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_20 testdata/ Test1.val
$ ./batch_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_40 testdata/ Test1.val
$ ./batch_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_60 testdata/ Test1.val

# SGE cluster
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train testdata/ Test1.val
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_20 testdata/ Test1.val
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_40 testdata/ Test1.val
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_60 testdata/ Test1.val

# Check the results
$ ./nps_check.sh score testdata/Test1/npsdat/ Test1.train testdata/Test1/ Test1.val
$ ./nps_check.sh score testdata/Test1/npsdat/ Test1.train.win_20 testdata/Test1/ Test1.val
$ ./nps_check.sh score testdata/Test1/npsdat/ Test1.train.win_40 testdata/Test1/ Test1.val
$ ./nps_check.sh score testdata/Test1/npsdat/ Test1.train.win_60 testdata/Test1/ Test1.val

# Calculate the overall prediction accuray in the validation cohort 
# Same on clusters and for batch processing (no parallelization)
$ Rscript npsR/nps_val.R testdata/Test1/npsdat/ testdata/Test1/ testdata/Test1/Test1.val.5K.fam testdata/Test1/Test1.val.5K.phen 0 20 40 60 
```

```
Non-Parametric Shrinkage 1.0.0 
Validation cohort:
Total  5000 samples
271  case samples
4729  control samples
0  samples with missing phenotype (-9)
Includes TotalLiability
Checking a prediciton model (winshift = 0 )...
Observed-scale R2 = 0.08795757 
Liability-scale R2 = 0.4207207 
Checking a prediciton model (winshift = 20 )...
Observed-scale R2 = 0.08940181 
Liability-scale R2 = 0.4150352 
Checking a prediciton model (winshift = 40 )...
Observed-scale R2 = 0.08912039 
Liability-scale R2 = 0.4187129 
Checking a prediciton model (winshift = 60 )...
Observed-scale R2 = 0.09010319 
Liability-scale R2 = 0.4182834 



Producing a combined prediction model...OK (saved in testdata/Test1.val.5K.phen.nps_score )
Observed-scale R2 = 0.09048684 
Liability-scale R2 = 0.4244738 
Loading required package: pROC
Type 'citation("pROC")' for a citation.

Attaching package: ‘pROC’

The following objects are masked from ‘package:stats’:

    cov, smooth, var

AUC:

Call:
roc.default(controls = prisk[vlY == 0], cases = prisk[vlY ==     1], ci = TRUE)

Data: 4729 controls < 271 cases.
Area under the curve: 0.8531
95% CI: 0.8321-0.8741 (DeLong)
Loading required package: DescTools
Nagelkerke's R2 = 0.2693255 

Call:  glm(formula = vlY ~ prisk, family = binomial(link = "logit"))

Coefficients:
(Intercept)        prisk  
    -7.2563       0.1908  

Degrees of Freedom: 4999 Total (i.e. Null);  4998 Residual
Null Deviance:	    2107 
Residual Deviance: 1621 	AIC: 1625
Done

```

### Running NPS on Test set #2

```bash
bsub -J stdgt[1-22] lsf/nps_stdgt.job testdata/Test2/ Test2.train 5000

./nps_check.sh stdgt testdata/Test2/ Test2.train 

Rscript npsR/nps_init.R testdata/Test2/Test2.summstats.txt testdata/Test2 testda
ta/Test2/Test2.train.2.5K_2.5K.fam testdata/Test2/Test2.train.2.5K_2.5K.phen Tes
t2.train 4000 testdata/Test2/npsdat 

./nps_check.sh init testdata/Test2/npsdat/

bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/T
est2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/T
est2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/T
est2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/T
est2/npsdat/ 3000 

./nps_check.sh decor testdata/Test2/npsdat/ 0 
./nps_check.sh decor testdata/Test2/npsdat/ 1000 
./nps_check.sh decor testdata/Test2/npsdat/ 2000 
./nps_check.sh decor testdata/Test2/npsdat/ 3000 

bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/T
est2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/T
est2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/T
est2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/T
est2/npsdat/ 3000 

./nps_check.sh prune testdata/Test2/npsdat/ 0 
./nps_check.sh prune testdata/Test2/npsdat/ 1000 
./nps_check.sh prune testdata/Test2/npsdat/ 2000 
./nps_check.sh prune testdata/Test2/npsdat/ 3000 

bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 3000 

./nps_check.sh gwassig testdata/Test2/npsdat/ 0
./nps_check.sh gwassig testdata/Test2/npsdat/ 1000
./nps_check.sh gwassig testdata/Test2/npsdat/ 2000
./nps_check.sh gwassig testdata/Test2/npsdat/ 3000

Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 0 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 1000 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 2000 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 3000 10 10

./nps_check.sh prep_part testdata/Test2/npsdat/ 0
./nps_check.sh prep_part testdata/Test2/npsdat/ 1000
./nps_check.sh prep_part testdata/Test2/npsdat/ 2000
./nps_check.sh prep_part testdata/Test2/npsdat/ 3000

bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 3000 

./nps_check.sh part testdata/Test2/npsdat/ 0
./nps_check.sh part testdata/Test2/npsdat/ 1000
./nps_check.sh part testdata/Test2/npsdat/ 2000
./nps_check.sh part testdata/Test2/npsdat/ 3000

Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 0 
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 1000 
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 2000 
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 3000 

Rscript npsR/nps_plot_shrinkage.R testdata/Test2/npsdat/ Test2.nps.pdf 0 1000 2000 3000 

Rscript npsR/nps_train_AUC.R testdata/Test2/npsdat/ 0 1000 2000 3000 
```

```
Data: 2500 controls < 2500 cases.
Area under the curve: 0.7843
95% CI: 0.7718-0.7968 (DeLong)
```

```bash
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 3000 

./nps_check.sh back2snpeff testdata/Test2/npsdat/ 0 
./nps_check.sh back2snpeff testdata/Test2/npsdat/ 1000 
./nps_check.sh back2snpeff testdata/Test2/npsdat/ 2000 
./nps_check.sh back2snpeff testdata/Test2/npsdat/ 3000 

bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ Test2.train testdata/Test2/ Test2.val
bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ Test2.train.win_1000 testdata/Test2/ Test2.val
bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ Test2.train.win_2000 testdata/Test2/ Test2.val
bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ Test2.train.win_3000 testdata/Test2/ Test2.val

./nps_check.sh score testdata/Test2/npsdat/ Test2.train testdata/Test2/ Test2.val
./nps_check.sh score testdata/Test2/npsdat/ Test2.train.win_1000 testdata/Test2/ Test2.val
./nps_check.sh score testdata/Test2/npsdat/ Test2.train.win_2000 testdata/Test2/ Test2.val
./nps_check.sh score testdata/Test2/npsdat/ Test2.train.win_3000 testdata/Test2/ Test2.val

Rscript npsR/nps_val.R testdata/Test2/npsdat/ testdata/Test2/ testdata/Test2/Test2.val.5K.fam testdata/Test2/Test2.val.5K.phen 0 1000 2000 3000 
```

```
Non-Parametric Shrinkage 1.0.0 
Validation cohort:
Total  5000 samples
240  case samples
4760  control samples
0  samples with missing phenotype (-9)
Includes TotalLiability
Checking a prediciton model (winshift = 0 )...
Observed-scale R2 = 0.04862955 
Liability-scale R2 = 0.2303062 
Checking a prediciton model (winshift = 1000 )...
Observed-scale R2 = 0.04994584 
Liability-scale R2 = 0.2298484 
Checking a prediciton model (winshift = 2000 )...
Observed-scale R2 = 0.05150205 
Liability-scale R2 = 0.2268046 
Checking a prediciton model (winshift = 3000 )...
Observed-scale R2 = 0.05258402 
Liability-scale R2 = 0.2298871 



Producing a combined prediction model...OK (saved in testdata/Test2/Test2.val.5K.phen.nps_score )
Observed-scale R2 = 0.05253146 
Liability-scale R2 = 0.2376991 
Loading required package: pROC
Type 'citation("pROC")' for a citation.

Attaching package: ‘pROC’

The following objects are masked from ‘package:stats’:

    cov, smooth, var

AUC:

Call:
roc.default(controls = prisk[vlY == 0], cases = prisk[vlY ==     1], ci = TRUE)

Data: 4760 controls < 240 cases.
Area under the curve: 0.7886
95% CI: 0.7617-0.8154 (DeLong)
Loading required package: DescTools
Nagelkerke's R2 = 0.1668188 

Call:  glm(formula = vlY ~ prisk, family = binomial(link = "logit"))

Coefficients:
(Intercept)        prisk  
    -5.2888       0.2387  

Degrees of Freedom: 4999 Total (i.e. Null);  4998 Residual
Null Deviance:      1926 
Residual Deviance: 1652         AIC: 1656
Done
```
