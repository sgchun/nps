
# Non-Parametric Shrinkage (NPS)
NPS is a non-parametric polygenic risk prediction model described in Chun et al. (2018) [(preprint)](https://www.biorxiv.org/content/early/2018/07/16/370064). NPS starts with a set of summary statistics in the form of SNP effect sizes from a large GWAS cohort. It then removes the correlation structure across summary statistics arising due to linkage disequilibrium and applies a piecewise linear interpolation on conditional mean effects. The conditional mean effects are estimated by partitioning-based non-parametric shrinkage algorithm using a training cohort with individual-level genotype data. 

For inquiries on using the software, please contact: Sung Chun (sgchun@bwh.harvard.edu), Nathan Stitziel (nstitziel@wustl.edu) or Shamil Sunyaev (ssunyaev@rics.bwh.harvard.edu). 

For citation: 
> Chun et al. Non-parametric polygenic risk prediction using partitioned GWAS summary statistics. 
> BioRxiv 370064, doi: https://doi.org/10.1101/370064 (preprint).

## How to Install
1. Download and unpack NPS package as below. Some of NPS codes are optimized in C++ and need to be compiled with GNU C++ compiler (GCC-4.4 or later). This will create two executable binaries, **stdgt** and **grs**, in the top-level NPS directory. **stdgt** is used to convert allelic dosages to standardized genotypes with the mean of 0 and variance of 1. **grs** calculates genetic risk scores using per-SNP genetic effects computed by NPS.

   ```bash
   tar -zxvf nps-1.0.0.tar.gz
   cd nps-1.0.0/
   make
   ```

   **Note on computer clusters: If you loaded a GCC module to compile NPS, you need to load the module also in job scripts - `nps_stdgt.job` and `nps_score.job` - as stdgt and grs will depend on GCC shared libraries in the run time.**

2. The core NPS module was implemented in R and requires R-3.0 or later (available for download from [here](https://www.r-project.org/)). Although NPS can run on a standard version of R, we strongly recommend to use R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [R open](https://mran.microsoft.com/open). These libraries can speed up NPS substantially.  

3. (*Optional*) NPS relies on R modules, [pROC](https://cran.r-project.org/web/packages/pROC/index.html) and [DescTools](https://cran.r-project.org/web/packages/DescTools/index.html), to calculate the AUC and Nagelkerke's *R2* statistics. These modules are optional; if they are not installed, AUC and Nagelkerke's *R2* will simply not be reported. To enable this feature, please install these packages by running the following on command line: 

   ```bash
   Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
   Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
   ```

   In case that you prefer to install the R extensions in your home directory (e.g. ~/R), please do the following instead:

   ```bash
   Rscript -e 'install.packages("pROC", "~/R", repos="http://cran.r-project.org")' 
   Rscript -e 'install.packages("DescTools", "~/R", repos="http://cran.r-project.org")' 
   
   # Add "~/R" to the local R library path in your login shell's start-up file.
   # For example, in case of bash, add the following to .bash_profile or .bashrc: 
   export R_LIBS="~/R:$R_LIBS"
   ```

4. Although we provide a command line tool to run NPS on desktop computers without parallelization (see `run_all_chroms.sh`), we strongly recommend to run it on computer clusters processing all chromosomes in parallel. To make this easier, we provide job scripts for SGE and LSF clusters (see `sge/` and `lsf/` directories). You may still need to modify the provided job scripts to load necessary modules similarly as the following example of [sge/nps_score.job](https://github.com/sgchun/nps/blob/master/sge/nps_score.job):

   ```bash
   ###
   # ADD CODES TO LOAD MODULES HERE
   #
   # Load R module if necessary. 
   # 
   # If you loaded a GCC module to compile grs, you also need to load 
   # the GCC module here. 
   #
   # ---------------------- EXAMPLE ----------------------------
   # On clusters running environment modules and providing R-mkl
   module add gcc/5.3.0 
   module add R-mkl/3.3.2
   
   # On clusters running DotKit instead and supporting OpenblasR
   use GCC-5.3.0 
   use OpenblasR
   # -----------------------------------------------------------
   ...
   ```

   **Note: Do not blindly add the above lines. The details will depend on individual system configurations.** 

5. We provide job scripts to prepare training and validation cohorts for NPS. These scripts require [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) and [QCTOOL v2](https://www.well.ox.ac.uk/~gav/qctool/). Please modify the job scripts (`ukbb_support/*.job`) to load **bgen** and **qctool** modules as necessary.

## Input files for NPS
To run NPS, you need the following set of input files: 

1. **GWAS summary statistics.** NPS supports two summary statistics formats: *minimal* and *preformatted*. 
   - The summary statistics in the *minimal* format can be automatically converted into the *preformatted* format and harmonized with training genotype data using provided NPS scripts (See [here](https://github.com/sgchun/nps#how-to-prepare-training-and-validation-cohorts-for-nps) for the step-by-step instruction). The minimal format is a tab-delimited text file with the following seven or eight columns: 
     - **chr**: chromosome number. NPS expects only chromosomes 1-22.
     - **pos**: base position of SNP.
     - **a1** and **a2**: Alleles at each SNP in any order.
     - **effal**: effect allele. It should be either a1 or a2 allele. 
     - **pval**: p-value of association. 
     - **effbeta**: estimated *per-allele* effect size of *the effect allele*. For case/control GWAS, log(OR) should be used. 
     - **effaf**: (*Optional*) allele frequency of *effect allele* in the discovery GWAS cohort. If this column is provided, the effect allele frequency will be compared between GWAS data and training cohort, and markers with too discrepant allele frequencies will be filtered out. If this information is not available, allele frequencies of training cohort will be copied to the summary statistics file. 
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
     
   - The *preformatted* format is the native format for NPS. We provide all summary statistics of our [test cases](https://github.com/sgchun/nps#test-cases) in this format. This is a tab-delimited text file format, and rows are sorted by chromsome numbers and positions. The following seven columns are required: 
     - **chr**: chromosome name starting with "chr." NPS expects only chromosomes chr1-chr22.
     - **pos**: base position of SNP.
     - **ref** and **alt**: reference and alternative alleles of SNP. NPS does not allow InDels, tri-allelic SNPs, or duplicated markers. 
     - **reffreq**: allele frequency of reference allele in the discovery GWAS cohort. 
     - **pval**: p-value of association. 
     - **effalt**: estimated *per-allele* effect size of *the alternative allele*. For case/control GWAS, log(OR) should be used. NPS will convert **effalt** to effect sizes relative to the standardized genotype using **reffreq**.  
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

2. **Training genotypes in the QCTOOL dosage format.** Genotype data of the training cohort are expected to follow the "dosage" format. We use [QCTOOL](https://www.well.ox.ac.uk/~gav/qctool/) to generate these files (See [instructions](https://github.com/sgchun/nps#how-to-prepare-training-and-validation-cohorts-for-nps)). The genotypes are split by chromosomes, and for each chromosome, the file is named as "chrom*N*.*TrainSetTag*.dosage.gz." These files are space-delimited compressed text files with the first six columns specifying the marker and rest of columns reporting its allelic dosage in each individual as follows:

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
   
   In order to match SNPs between GWAS summary statistics and training and validation cohorts, NPS relies on the combination of chromosome, base position and alleles. All alleles are designated on the forward strand (+). **SNPID** or **rsid** will be ignored. The **alleleA** has to match **ref** allele, and the **alleleB** has to match **alt** allele in the GWAS summary statistics. The allelic dosage counts the genetic dosage of **alleleB** in each individual. Markers has to be ordered by **position**. All SNPs in the training genotype file are expected to have non-missing GWAS summary statistics.  

3. **Training sample IDs in PLINK .fam format.** The samples in the .fam file should appear in the exactly same order as the samples in the training genotype files. This is *space-separated* six-column text file without a header. The phenotype information in this file is ignored. See [here](https://www.cog-genomics.org/plink2/formats#fam) for additional details on the format. 
   ```
   trainF2 trainI2 0 0 0 -9
   trainF3 trainI3 0 0 0 -9
   trainF39 trainI39 0 0 0 -9
   trainF41 trainI41 0 0 0 -9
   trainF58 trainI58 0 0 0 -9
   ```
4. **Training phenotypes in PLINK phenotype format.** NPS looks up phenotypes in a separately prepared phenotype file. The phenotype name has to be **Outcome** with cases and controls encoded by **1** and **0**, respectively. The combination of **FID** and **IID** are used to match samples with .fam file. This file is *tab-delimited*, and samples can appear in any order. Missing phenotypes (e.g. encoded with **-9** or missing entry of samples described in .fam file) are not allowed.
   ```
   FID	IID	Outcome
   trainF2	trainI2  1
   trainF39 trainI39 1
   trainF3	trainI3  0
   trainF41 trainI41 0
   trainF58 trainI58 0
   ```
5. **Validation genotypes in QCTOOL dosage format.** Same as the training genotype dosage format. 
6. **Validation sample IDs in PLINK .fam format.** Same as the training sample ID file format.
7. **Validation phenotypes in PLINK phenotype format.** Similar to the training phenotype file format. Unlike the training cohort phenotype file, missing phenotypes (encoded by **-9**) are allowed in a validation cohort phenotype file. Samples with missing phenotypes will be simply excluded when evaluating the accuracy of prediction model.

## How to prepare training and validation cohorts for NPS
We take an example of UK Biobank to show how to prepare training and validation cohorts for NPS. In principle, however, NPS can work with other cohorts as far as the genotype data are prepared in .bgen file format. To gain access to UK Biobank data, please see [UK Biobank data access application procedure](https://www.ukbiobank.ac.uk/). 

### Using UK Biobank as a training cohort
UK Biobank data consist of the following files: 
- **ukb_imp_chrN_v3.bgen**: imputed allelic for each chromosome chrN
- **ukb_mfi_chrN_v3.txt**: information of markers in the bgen file
- **ukb31063.sample**: bgen sample information file 

Assuming that UK Biobank dataset is located in `<path_to_ukbb>/` directory, we first exclude SNPs with minor allele frequency < 5% or imputation quality (INFO) score < 0.4 by running the following:   
```bash
qsub -l h_vmem=4G -t 1-22 ukbb_support/common_snps.job <path_to_ukbb>/ukb_imp_chr#_v3.bgen <path_to_ukbb>/ukb_mfi_chr#_v3.txt <work_dir>
```
The output files will be stored in `<work_dir>/`.

Next, we filter bgen files to include only samples with identifier listed in `<sample_id_file>` and convert bgen files to dosage file format. `<sample_id_file>` is simply a list of sample IDs, with one sample in each line, and looks like the following: 
```
2959669
1774228
3227484
...
```

With the `<sample_id_file>`, the following step will generate filtered dosage files using qctool in `<work_dir>/`: 
```bash
qsub -l h_vmem=4G -t 1-22 ukbb_support/filter_samples.job <path_to_ukbb>/ukb31063.sample <work_dir> <sample_id_file> <training_cohort_name>
```  

Then, we harmonize GWAS summary statitics with training cohort data. When `<summary_statistics_file>` is a raw summary statistics file in the *minimal* format, running the following step will harmonize it with training cohort data and generate the harmonized GWAS summary statistics in the *preformatted* format: 
```bash
Rscript ukbb_support/harmonize_summstats.R <summary_statistics_file> <work_dir> <training_cohort_name>
```

We need to filter out SNPs that were flagged for removal during the harmonization step by running the following line: 
```bash
qsub -l h_vmem=4G -t 1-22 ukbb_support/filter_variants.job <work_dir> <training_cohort_name>
```

Finally, the following step will create `<work_dir>/<training_cohort_name>.fam`, which keeps tracks of the IDs of all training cohort samples: 
```bash
ukbb_support/make_fam.sh <work_dir> <training_cohort_name>
```

After the above steps, ukbb_support job scripts will automatially generate the following files ready for NPS: 
- `<work_dir>/<training_cohort_name>.preformatted_summstats.txt`
- `<work_dir>/chromN.<training_cohort_name>.QC2.dosage.gz`
- `<work_dir>/<training_cohort_name>.fam`

* **Note: common_snps.job and filter_samples.job will use bgenix and qctool, respectively. The job scripts may need to be moditifed to load these modules.** 
* **Note: The job scripts in `ukbb_support/` directory is for SGE clusters but can be easily modified for LSF or other cluster systems.**
* **Note: Some steps take long running time and demand memory space up to 4GB (`qsub -l h_vmem=4G`). `ukbb_support/harmonize_summstats.R` may take memory up to ~8GB. `ukbb_support/harmonize_summstats.R` will terminate abruptly without explanation if it runs out of memory.**

### Using other cohort as a training cohort 

```bash

```

### Using UK Biobank for validation as well as training cohorts

```
cp <training_cohort_name>.UKBB_rejected_SNPIDs <validation_cohort_name>.UKBB_rejected_SNPIDs
qsub -t 1-22 ukbb_support/filter_samples.job <path_to_ukbb>/ukb31063.sample <work_dir> <sample_id_file> <validation_cohort_name>
qsub -t 1-22 ukbb_support/filter_variants.job <work_dir> <validation_cohort_name>
ukbb_support/make_fam.sh <work_dir> <validation_cohort_name>
```

### Using a different cohort for validation


## Running NPS 

### Test cases
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

For Test set #1, we provide an instruction on running it on desktop without parallelization or on SGE clusters. For LSF clusters, see the example of [Test set #2](https://github.com/sgchun/nps#running-nps-on-test-set-2). 

1. **Standardize genotypes.** The first step is to standardize the training genotypes to the mean of 0 and variance of 1. `sge/snps_stdgt.job` and `lsf/snps_stdgt.job` scripts run the computation in parallel across all 22 chromosomes. On a desktop computer, `run_all_chroms.sh` can be used to run an SGE job script for all chromosomes one by one sequentially. The first parameter (`testdata/Test1`) is the location of training cohort (.dosage.gz files), the second parameter (`Test1.train`) is the *GenotypeSetID* of training cohort, and the last parameter (`5000`) is the number of samples in the training cohort. 
```bash
$ cd nps-1.0.0/

# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_stdgt.job testdata/Test1 Test1.train

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test1 Test1.train
```
After all jobs are completed, `nps_check.sh` script can be used to make sure that all jobs are successful. If failure is detected, `FAIL` message will be printed. Note `nps_check.sh` script can be used on clusters (both SGE and LSF) and with serial processing as follows in the same way: 
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
# Same on clusters and for serial processing (no parallelization)
$ Rscript npsR/nps_init.R testdata/Test1/Test1.summstats.txt testdata/Test1 testdata/Test1/Test1.train.2.5K_2.5K.fam testdata/Test1/Test1.train.2.5K_2.5K.phen Test1.train 80 testdata/Test1/npsdat

# Check the results
$ ./nps_check.sh init testdata/Test1/npsdat/
Checking testdata/Test1/npsdat//args.RDS ...OK (version 1.0.0)
Checking testdata/Test1/npsdat//log ...OK
```

3. **Transform data to the decorrelated "eigenlocus" space.** In general, this is one of the most time-consuming steps in NPS. We recommend to run NPS four times on shifted windows and merge the results in the later steps. For the window shift, we recommend the shifts of 0, 1/WINSZ, 2/WINSZ and 3/WINSZ SNPs, where WINSZ is the size of analysis window. In the case of test set #1, they correspond to 0, 20, 40 and 60. The first parameter is the location of intermediary data (`testdata/Test1/npsdat/`), and the second parameter is the window shift (`0`, `20`, `40` or `60`). 
```bash
# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 0
$ ./run_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 20
$ ./run_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 40
$ ./run_all_chroms.sh sge/nps_decor.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 0 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 20 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 40 
$ qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 60 

# Check the results
$ ./nps_check.sh decor testdata/Test1/npsdat/ 0 20 40 60 
```

4. **Prune correlations across windows.** This step prunes the correlation between genotypes across adjacent windows in the eigenlocus space.  
```bash
# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 0
$ ./run_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 20
$ ./run_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 40
$ ./run_all_chroms.sh sge/nps_prune.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_prune.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh prune testdata/Test1/npsdat/ 0 20 40 60
```

5. **Separate GWAS-significant partition.** The partition of GWAS-significant associations will be separated out from the rest of association signals. NPS takes longer time to complete this step when there are more GWAS-significant signals.
```bash
# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 0
$ ./run_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 20
$ ./run_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 40
$ ./run_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh gwassig testdata/Test1/npsdat/ 0 20 40 60
```

6. **Define a partitioning scheme.** The partition scheme will be defined with `npsR/nps_prep_part.R` and then, partitioned genetic risk scores will be calculated for all training samples using `nps_part.job`. Specifically, for `npsR/nps_prep_part.R`, the first parameter is the location of intermediary data (`testdata/Test1/npsdat/`), the second is the window shift (`0`, `20`, `40` or `60`), the third is the number of partitions on intervals of eigenvalues of eigenlocus projection (`10`), and the last is the number of partitions on intervals of observed effect sizes in the eigenlocus space (`10`). For `nps_part.job`, the first parameter is the location of intermediary data (`testdata/Test1/npsdat/`), and the second is the window shift (`0`, `20`, `40` or `60`)
```
# Same on clusters and for serial processing (no parallelization)
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 0 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 20 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 40 10 10 
$ Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 60 10 10 

# Check the results
$ ./nps_check.sh prep_part testdata/Test1/npsdat/ 0 20 40 60 

# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 0
$ ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 20
$ ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 40
$ ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh part testdata/Test1/npsdat/ 0 20 40 60
```

7. **Estimate per-partition shrinkage weights.** Then, we estimate the per-partition shrinkage weights using `npsR/nps_weight.R`. We also provide two optional utilities: `npsR/nps_train_AUC.R`, which reports the AUC statistics of prediction in training cohort, and `npsR/nps_plot_shrinkage.R`, which plots the overall curve of GWAS effect sizes re-weighted by per-partition shrinkage. Both tools take the average of NPS run on shifted windows. `npsR/nps_plot_shrinkage.R` will save the shrinkage curve plot in the pdf file path given as second argument (`Test1.nps.pdf`). 

```bash
# Same on clusters and for serial processing (no parallelization)
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 0 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 20 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 40 
$ Rscript npsR/nps_weight.R testdata/Test1/npsdat/ 60 

# Check the results
$ ./nps_check.sh weight testdata/Test1/npsdat/ 0 20 40 60
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
# Same on clusters and for serial processing (no parallelization)
$ Rscript npsR/nps_plot_shrinkage.R testdata/Test1/npsdat/ Test1.nps.pdf 0 20 40 60
```

8. **Convert back to per-SNP effect sizes.** Then, the re-weighted effect sizes should be converted back to the original per-SNP space from the eigenlocus space. This will generate Test1.train.adjbetahat.chrom*N*.txt for the shift of 0, and Test1.train.win_*shift*.adjbetahat.chrom*N*.txt for the rest of shifts in `testdata/Test1/npsdat/` directory. 
```bash
# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 0
$ ./run_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 20
$ ./run_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 40
$ ./run_all_chroms.sh sge/nps_back2snpeff.job testdata/Test1/npsdat/ 60

# Or on SGE cluster
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 0
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 20
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 40
$ qsub -cwd -t 1-22 sge/nps_back2snpeff.job testdata/Test1/npsdat/ 60

# Check the results
$ ./nps_check.sh back2snpeff testdata/Test1/npsdat/ 0 20 40 60
```

9. **Validate the accuracy of prediction model in a validation cohort.** Polygenic risk scores will be calculated for each chromosome for all individuals in the validation cohort, named `Test1.val`, found in `testdata/` directory. The per-SNP effect sizes calculated by NPS will be looked up in `testdata/Test1/npsdat/` folder, and they will be designed by `Test1.train` for no window shift, and `Test1.train.win_#` for shifted windows. Then, `npsR/nps_val.R` script will be used to merge all information and compute the accuracy statistics. `npsR/nps_val.R`will require the NPS work directory (`testdata/Test1/npsdat/`), directory containing validation cohort data (`testdata/Test1/`), file containing sample IDs in the validation cohort (`testdata/Test1/Test1.val.5K.fam`), phenotypes for validation samples (`testdata/Test1/Test1.val.5K.phen`), and window shifts used in the prediction model (`0`, `20`, `40` or `60`). 

```bash
# Serial processing (on desktop)
$ ./run_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0 
$ ./run_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 20
$ ./run_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 40
$ ./run_all_chroms.sh sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 60

# SGE cluster
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 20
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 40
$ qsub -cwd -t 1-22 sge/nps_score.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 60

# Check the results
$ ./nps_check.sh score testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0 20 40 60

# Calculate the overall prediction accuray in the validation cohort 
# Same on clusters and for serial processing (no parallelization)
$ Rscript npsR/nps_val.R testdata/Test1/npsdat/ testdata/Test1/ testdata/Test1/Test1.val.5K.fam testdata/Test1/Test1.val.5K.phen 0 20 40 60 
```
`npsR/nps_val.R` will print out the following output on Test case #1. It reports here the liability-scale R2 of 0.4244738 and AUC of 0.8531, Nagelkerke's R2 of 0.2693255. The individual NPS polygenic score is stored in the file `testdata/Test1.val.5K.phen.nps_score`. 

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

We provide the instruction to run on Test set #2 with LSF. Running it on SGE clusters will be similar. For ~5,000,000 genome-wide SNPs, we recommend to use the window size of 4,000 SNPs, and window shifts of 0, 1,000, 2,000, and 3,000 SNPs. In addition, some of NPS steps now require large memory space. For the most of data we tested, 4GB memory limit is sufficient to run individual NPS tasks. The memory limit can be specified by `-R 'rusage[mem=4000]'` in LSF and `-l h_vmem=4G` in SGE clusters. 

```bash
bsub -J stdgt[1-22] lsf/nps_stdgt.job testdata/Test2/ Test2.train

./nps_check.sh stdgt testdata/Test2/ Test2.train 

Rscript npsR/nps_init.R testdata/Test2/Test2.summstats.txt testdata/Test2 testdata/Test2/Test2.train.2.5K_2.5K.fam testdata/Test2/Test2.train.2.5K_2.5K.phen Test2.train 4000 testdata/Test2/npsdat 

./nps_check.sh init testdata/Test2/npsdat/

bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J decor[1-22] lsf/nps_decor.job testdata/Test2/npsdat/ 3000 

./nps_check.sh decor testdata/Test2/npsdat/ 0 1000 2000 3000

bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J prune[1-22] lsf/nps_prune.job testdata/Test2/npsdat/ 3000 

./nps_check.sh prune testdata/Test2/npsdat/ 0 1000 2000 3000

bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] lsf/nps_gwassig.job testdata/Test2/npsdat/ 3000 

./nps_check.sh gwassig testdata/Test2/npsdat/ 0 1000 2000 3000

Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 0 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 1000 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 2000 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 3000 10 10

./nps_check.sh prep_part testdata/Test2/npsdat/ 0 1000 2000 3000

bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J part[1-22] lsf/nps_part.job testdata/Test2/npsdat/ 3000 

./nps_check.sh part testdata/Test2/npsdat/ 0 1000 2000 3000

Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 0 
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 1000 
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 2000 
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 3000 

Rscript npsR/nps_plot_shrinkage.R testdata/Test2/npsdat/ Test2.nps.pdf 0 1000 2000 3000 

Rscript npsR/nps_train_AUC.R testdata/Test2/npsdat/ 0 1000 2000 3000 
```

The following is the training AUC reported by `npsR/nps_train_AUC.R` with test case #2. 
```
Data: 2500 controls < 2500 cases.
Area under the curve: 0.7843
95% CI: 0.7718-0.7968 (DeLong)
```

The validation can be proceeded as follows. 
```bash
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 0 
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 1000 
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 2000 
bsub -R 'rusage[mem=4000]' -J back2snpeff[1-22] lsf/nps_back2snpeff.job testdata/Test2/npsdat/ 3000 

./nps_check.sh back2snpeff testdata/Test2/npsdat/ 0 1000 2000 3000

bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0
bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 1000
bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 2000
bsub -J score[1-22] lsf/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 3000

./nps_check.sh score testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0 1000 2000 3000

Rscript npsR/nps_val.R testdata/Test2/npsdat/ testdata/Test2/ testdata/Test2/Test2.val.5K.fam testdata/Test2/Test2.val.5K.phen 0 1000 2000 3000 
```

`npsR/nps_val.R` reports the following performance statistics in the validation cohort. 
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


