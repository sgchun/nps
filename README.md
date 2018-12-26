
# Non-Parametric Shrinkage (NPS)
NPS is a non-parametric polygenic risk prediction algorithm described in Chun et al. (2018) BioRxiv [(preprint)](https://www.biorxiv.org/content/early/2018/07/16/370064). 

## How to install

1. The core NPS module is implemented in R. R can be downloaded from [here](https://www.r-project.org/). R-3.0 or later is required to run NPS. Although NPS can run on a plain-vanilla version of R, we strongly recommend to use R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [R open](https://mran.microsoft.com/open). Since the most time-consuming steps in NPS are matrix manipulations, this can significantly reduce the running time of NPS. 

2. (Optional) NPS relies on R modules, `pROC` and `DescTools`, to calculate the AUC and Nagelkerke's R2 statistics. These modules are optional; if they are not installed, AUC and Nagelkerke's R2 will not be reported. To enable this feature, please install these packages by running the following on command line: 

```
$ Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
$ Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
```

In case that it is preferred to install these R extensions in your home directory (e.g. ~/R) instead of the default system path, please do the following instead:
```
$ Rscript -e 'install.packages("pROC", "~/R", repos="http://cran.r-project.org")' 
$ Rscript -e 'install.packages("DescTools", "~/R", repos="http://cran.r-project.org")' 

# Add "~/R" to your local R library path in the login shell start up file 
# For examaple in case of bash, add this: 
$ export R_LIBS=~/R:$R_LIBS
```

3. Download and unpack NPS package as below. Some of NPS codes are optimized in C++ and need to be compiled. For this, you need GNU C++ compiler with C++0x support (GCC 4.4 or later). This step will create two executable binaries, `stdgt` and `grs`, in the top-level NPS directory. `stdgt` converts an imputed dosage file to standardized genotypes with the mean of 0 and variance of 1. `grs` calculates genetic risk scores for all indiviudals in an imputed dosage file using per-SNP genetic effects computed by NPS.

```
$ tar -zxvf nps-1.0.0.tar.gz
$ cd nps-1.0.0/
$ make
```

4. We recommend to run NPS on computer clusters, processing all chromosomes in parallel. To make this easier, we provide job scripts for SGE and LSF clusters. Please see `sge` and `lsf` directories along with provided examples below. You may still need to modify the provided job scripts, for example, to load necessary modules. 

5. (Optional) We provide sample SGE job scripts to prepare UK Biobank data for NPS training and validation cohorts. To gain access to UK Biobank data, please see [UK Biobank data access application procedure](https://www.ukbiobank.ac.uk/). Our scripts rely on [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) and [QCTOOL v2](https://www.well.ox.ac.uk/~gav/qctool/). If you use different cohorts, please convert them to .bgen file and follow our instructions. 

## Input files for NPS
To run NPS, you need the following set of input files: 

1. **GWAS summary statistics.** This is a tab-delimited text file. NPS requires the following seven essential columns: 
- `chr`: chromosome name starting with "chr." NPS expects only chromosomes chr1-chr22.
- `pos`: base positions of SNP.
- `ref` and `alt`: reference and alternative alleles of a SNP. NPS does not allow InDels or tri-allelic SNPs. There should not be duplicated SNPs in the file. 
- `reffreq`: allele frequency of reference allele in the discovery GWAS cohort. 
- `pval`: p-value of association. 
- `effalt`: estimated *per-allele* effect size of alternative allele. For case/control GWAS, log(OR) should be used for this. NPS will convert this value to an effect size relative to standardized genotype using `reffreq`.  

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
2. **Training genotypes in QCTOOL dosage format.** Genotype data have to be prepared in the dosage format. NPS requires that all markers in this file to overlap with GWAS summary statitics. We recommend to remove InDels, tri-allelic SNPs, rare variants with MAF < 5%, markers with any QC issue, and markers that are not found in the GWAS summary statistics. Markers with very different allele frequencies between GWAS and training cohort should be also discarded. NPS does not allow duplicated SNPs in the dosage file. Genotype data needs to be split by chromosomes for parallelization, and each file should be named as "chrom*N*.*CohortName*.dosage.gz." NPS matches SNPs using the combination of chromosome, position and alleles on the assumed forward strand (+) and does not rely on `SNPID` or `rsid`. NPS expects `alleleA` to match `ref` allele and `alleleB` to match `alt` allele in the GWAS summary statistics. The allelic dosage counts the genetic dosage of `alleleB`. [QCTOOL](https://www.well.ox.ac.uk/~gav/qctool/) can be used to generate the dosage files.  
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

4. **Training phenotypes in PLINK phenotype format.** NPS looks up phenotypes in the separate phenotype file and ignores the phenotype data provided in the .fam file above. The phenotype name has to be `Outcome` with cases and controls encoded with `1` and `0`, respectively. When the optional `TotalLiability` column is provided, NPS also reports *R2* on the liability scale, namely *R2* between polygenic scores and `TotalLiablity`. The combination of `FID` and `IID` are used to match samples to .fam file. Samples can appear in any order in this file. Missing phenotypes (e.g. encoded with `-9` or missing samples listed in .fam file) are not allowed. See [here](http://zzz.bwh.harvard.edu/plink/data.shtml#pheno) for additional details on the format.
```
FID	IID	TotalLiability	Outcome
trainF68266	trainI68266	2.00037013482937	1
trainF77481	trainI77481	2.19146817328717	1
trainF39184	trainI39184	0.833227631267031	0
trainF76746	trainI76746	0.168893429255002	0
trainF65453	trainI65453	0.817881328612934	0
...
```
5. **Validation genotypes in QCTOOL dosage format.** This is for validation cohort. See the above for the explanation. 
6. **Validation sample IDs in PLINK .fam format.** This is for validation cohort. See the above for the explanation. 
7. **Validation phenotypes in PLINK phenotype format.** This is for validation cohort. See the above for the explanation. Missing phenotypes are allowed in this file (encoded with `-9`). Such samples will be simply excluded when calculating the accuracy statistics. 


## Test cases
We provide two sets of simulated test cases. Due to their file large sizes, they are provided separately from the software distribution. Please download them from [Sunyaev Lab website (LINKS WILL BE PROVIDED SOON)]. Test set #1 is a relatively small dataset (225MB), and NPS can complete in less than 1 hour in total on a modest desktop PC without linear-algebra acceleration. In contrast, test set #2 is more realistic simulation (11GB) and will require serious computational resource: NPS will generate up to 1 TB of intermediary data and take 1/2 day to a day to run on computer clusters (the exact overall running time will depend on the degree of parallelization).  

Both simulation datasets were generated using our multivariate-normal simulator. See our NPS manuscript for the details.   
- **Test set #1.** The number of markers across the genome is limited to 100,449. We assume that all causal SNPs are included in the 100,449 SNPs. Note that this is unrealistic assumption; causal SNPs cannot be accurately tagged with such a sparse SNP set. The fraction of causal SNP is 0.005 (a total of 522 SNPs). The GWAS cohort size is 100,000. The training cohort has 2,500 cases and 2,500 controls. The valiation cohort is 5,000 samples without case over-sampling. The heritability is 0.5. The phenotype prevalence is 5%.  

- **Test set #2.** The number of markers across the genome is 5,012,500. The fraction of causal SNP is 0.001 (a total of 5,008 SNPs). The GWAS cohort size is 100,000. The training cohort has 2,500 cases and 2,500 controls. The valiation cohort is 5,000 samples without case over-sampling. The heritability is 0.5. The phenotype prevalence is 5%. This is one of the benchmark simulation datasets used in our manuscript, with 1/10-fold reduced validation cohort size. 

We assume that you download and unpack test datasets in the following directories.
```
$ cd nps-1.0.0/testdata/
$ tar -zxvf NPS.Test1.tar.gz 
# This will create test data in nps-1.0.0/testdata/Test1
$ tar -zxvf NPS.Test2.tar.gz 
# This will create test data in nps-1.0.0/testdata/Test2
```

### Running NPS on Test set #1

1. **Standardize genotypes.** The first step is to standardize the training genotypes to the mean of 0 and variance of 1. `sge/snps_stdgt.job` and `lsf/snps_stdgt.job` scripts run the computation in parallel across all 22 chromosomes. On a desktop computer, `batch_all_chroms.sh` can be used to run an SGE job script for all chromosomes one by one in batch. The first parameter (`testdata/Test1`) is the location of training cohort (.dosage.gz files), the second parameter (`Test1.train`) is the *CohortName* of training cohort, and the last parameter (`5000`) is the number of samples in the training cohort. 
```
$ cd nps-1.0.0/

# SGE cluster
$ qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test1 Test1.train 5000

# Batch processing (on desktop)
$ ./batch_all_chroms.sh sge/nps_stdgt.job testdata/Test1 Test1.train 5000
```
After all jobs are completed, `nps_check.sh` script can be used to make sure that all computations were finished successfully as follows (`FAIL` message will be displayed if any failure is detected): 
```
$ ./nps_check.sh stdgt testdata/Test1 Test1.train 
Verifying nps_stdgt:
Checking testdata/Test1/chrom1.Test2.train ...OK
Checking testdata/Test1/chrom2.Test2.train ...OK
Checking testdata/Test1/chrom3.Test2.train ...OK
...
```

2. **Configure an NPS run.**

```
$ Rscript npsR/nps_init.R testdata/Test1/Test1.summstats.txt testdata/Test1 testdata/Test1/Test1.train.2.5K_2.5K.fam testdata/Test1/Test1.train.2.5K_2.5K.phen Test1.train 80 testdata/Test1/npsdat

$ ./nps_check.sh init testdata/Test1/npsdat/
```

3. **Transform data to the decorrelated "eigenlocus" space.**
```
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 0 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 20 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 40 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 60 

$ ./nps_check.sh decor testdata/Test1/npsdat/ 0 
$ ./nps_check.sh decor testdata/Test1/npsdat/ 20 
$ ./nps_check.sh decor testdata/Test1/npsdat/ 40 
$ ./nps_check.sh decor testdata/Test1/npsdat/ 60 
```

4. **Prune.**
```
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 60

$ ./nps_check.sh prune testdata/Test1/npsdat/ 0 
$ ./nps_check.sh prune testdata/Test1/npsdat/ 20 
$ ./nps_check.sh prune testdata/Test1/npsdat/ 40 
$ ./nps_check.sh prune testdata/Test1/npsdat/ 60 
```

5. **Separate GWAS-significant partition.**
```
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 60

$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 0 
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 20 
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 40 
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 60 
```

6. **Define a partitioning scheme.**
```
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 0 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 20 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 40 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 60 10 10 

$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 0 
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 20 
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 40 
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 60 

$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 60

$ ./nps_check.sh part testdata/Test1/npsdat/ 0 
$ ./nps_check.sh part testdata/Test1/npsdat/ 20 
$ ./nps_check.sh part testdata/Test1/npsdat/ 40 
$ ./nps_check.sh part testdata/Test1/npsdat/ 60 
```

7. **Estimate per-partition shrinkage weights.**
```
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 0 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 20 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 40 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 60 
```

```
$ Rscript npsR/nps_train_AUC.R testdata/Test1/npsdat/ 0 20 40 60
Data: 2500 controls < 2500 cases.
Area under the curve: 0.8799
95% CI: 0.8707-0.8891 (DeLong)
```

```
$ Rscript npsR/nps_plot_shrinkage.R testdata/Test1/npsdat/ Test1.nps.pdf 0 20 40 60
```

8. **Convert back to per-SNP effect sizes.**
```
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 60

$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 0 
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 20 
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 40 
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 60 
```

9. **Validate.**
```
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train testdata/ Test1.val
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_20 testdata/ Test1.val
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_40 testdata/ Test1.val
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ Test1.train.win_60 testdata/ Test1.val

$ Rscript npsR/nps_val.R testdata/Test1/npsdat/ testdata/ testdata/Test1.val.5K.fam testdata/Test1.val.5K.phen 0 20 40 60 

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
