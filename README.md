
# Non-Parametric Shrinkage (NPS)
NPS implements a non-parametric polygenic risk prediction algorithm described in [Chun et al. 2020 (preprint)](https://doi.org/10.1101/370064). NPS projects genetic data into an orthogonal domain called "eigenlocus space". Then, it re-weights GWAS effect sizes by partitioning genetic variations into trenches and measuring the predictive power of each trench in an independent training cohort. To run NPS, two sets of data are required: GWAS summary statistics and small individual-level training cohort with both genotype and phenotype data. 

We recommend running NPS in parallelized computer clusters. We support the following platforms: 
* SGE / UGER
* LSF 
* Slurm
* MacOS (not recommended for genome-wide datasets)

For citation: 
> Chun et al. Non-parametric polygenic risk prediction using partitioned GWAS summary statistics.  
> BioRxiv 2020. doi: 10.1101/370064 (preprint).

For inquiries on software, please contact: 
* Sung Chun (SungGook.Chun@childrens.harvard.edu)
* Nathan Stitziel (nstitziel@wustl.edu) 
* Shamil Sunyaev (ssunyaev@rics.bwh.harvard.edu). 

## How to Install
1. Download and unpack NPS package ([version 1.1.1](https://github.com/sgchun/nps/archive/1.1.1.tar.gz)). Part of NPS codes are optimized in C++ and have to be compiled using GNU C++ compiler (GCC-4.4 or later). This will create two executable binaries, **stdgt** and **grs**, in the top-level directory. 

   ```bash
   tar -zxvf nps-1.1.1.tar.gz
   cd nps-1.1.1/
   make
   ```

2. The core NPS module was implemented in R (version 3.3 or higher required). Although NPS can run on a standard version of R, we strongly recommend using R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [Microsoft R open](https://mran.microsoft.com/open). These libraries can substantially speed up NPS operations.  

3. NPS relies on R libraries, [pROC](https://cran.r-project.org/web/packages/pROC/index.html) and [DescTools](https://cran.r-project.org/web/packages/DescTools/index.html), to report the accuracy of polygenic scores in AUC and Nagelkerke's R^2. These modules are optional; if they are not installed, AUC and Nagelkerke's R^2 calculation will be skipped. To install these packages, run the following on command line: 

   ```bash
   Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
   Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
   ```

   To install the R extensions in the home directory (e.g. ~/R) rather than in the default system path, use the following commands instead:

   ```bash
   Rscript -e 'install.packages("pROC", "~/R", repos="http://cran.r-project.org")' 
   Rscript -e 'install.packages("DescTools", "~/R", repos="http://cran.r-project.org")' 
   
   # Add "~/R" to the local R library path in your login shell's start-up file.
   # For example, in case of bash, add the following to .bash_profile or .bashrc: 
   export R_LIBS="~/R:$R_LIBS"
   ```

4. Although we provide a command line tool to run NPS on desktop computers [without parallelization], we strongly recommend running it on computer clusters, processing all chromosomes in parallel. To make this easier, we provide job script templates in the [nps-1.1.1/sge/](https://github.com/sgchun/nps/tree/master/sge) directory. These scripts run not only with SGE but also with UGER, LSF and Slurm schedulers. You may still need to modify the provided job scripts to configure and load necessary modules, for example:

   ```bash
   ###
   # ADD CODES TO LOAD MODULES HERE
   # ---------------------- EXAMPLE ----------------------------
   # On clusters running environment modules and providing R-mkl
   module add gcc/5.3.0 
   module add R-mkl/3.3.2
   # -----------------------------------------------------------
   ...
   ```

   **The details will depend on specific system configurations.** 

## NPS Test datsets
We provide two sets of simulated test cases. They are provided separately from the software distribution and can be downloaded from Sunyaev Lab FTP server (ftp://genetics.bwh.harvard.edu/download/schun/). Test set #1 is relatively small and can be easily run on desktop computer whereas test set #2 is a more realistic dataset but requires serious computational resource to run NPS.

Test set | Total # of simulated SNPs | # of simulated causal SNPs | NPS disk space requirement | NPS run time
--- | --- | --- | --- | --- 
#1 | 100,449 | 522 (0.5%) | XX GB | < 30 mins* 
#2 | 5,012,500 | 5,008 (0.1%) | 1 TB | 3-6 hours**

[*] On desktop computer, without parallelization  
[**] On computer clusters, parallelizing up to 88 CPUs, with linear algebra acceleration.  

Download the test dataset and unpacked it as below. See [File Formats](https://github.com/sgchun/nps/blob/master/FileFormats.md) for the description on the input file formats.
```bash
cd nps-1.1.1/testdata/

tar -zxvf NPS.Test1.tar.gz 
# This will create the following test data files in nps-1.1.1/testdata/Test1
# Test1/Test1.summstats.txt (PREFORMATTED GWAS summary statistics)
# Test1/Test1.train.2.5K_2.5K.fam (training cohort sample IDs)
# Test1/Test1.train.2.5K_2.5K.phen (training cohort phenotypes)
# Test1/chrom1.Test1.train.dosage.gz (training cohort genotypes)
# Test1/chrom2.Test1.train.dosage.gz (training cohort genotypes)
# ... 
# Test1/Test1.val.5K.fam (validation cohort sample IDs)
# Test1/Test1.val.5K.phen (validation cohort phenotypes)
# Test1/chrom1.Test1.val.dosage.gz (validation cohort genotypes)
# Test1/chrom2.Test1.val.dosage.gz (validation cohort genotypes)
# ... 
```

## Running NPS
Test set #1 is small enough to run on desktop computers (MacOS and Linux are supported) without parallel processing. We provide a wrapper script (`run_all_chroms.sh`) to drive cluster jobs sequentially, by processing one chromosome at a time. To run test set #1 on computer clusters, see the instructions page for [SGE] and [LSF] schedulers. 

1. **Standardize genotypes.** This step standardizes the genotype data of training cohort to the mean of 0 and variance of 1 . The training genotype files should be in [the dosage format](https://github.com/sgchun/nps/blob/master/FileFormats.md) and named as chrom*N*.*DatasetID*.dosage.gz. The command arguments are:
    * (1) directory where training genotype files are: `testdata/Test1`
    * (2) *DatasetID* of training genotype files: `Test1.train`. 
   ```bash
   cd nps-1.1.1/
   
   ./run_all_chroms.sh sge/nps_stdgt.job testdata/Test1 Test1.train
   ```

2. **Configure an NPS run.** For test set #1, which has ~100,000 genomewide SNPs, we recommend a window size of 80 SNPs. In general, for ~5,000,000 genome-wide SNPs we recommend to use 4,000-SNP windows. The command arguments are:
    * (1) GWAS summary statistics file: `testdata/Test1/Test1.summstats.txt`
    * (2) directory where training genotype files are: `testdata/Test1`
    * (3) sample information of training cohort: `testdata/Test1/Test1.train.2.5K_2.5K.fam`
    * (4) phenotypes information of training samples: `testdata/Test1/Test1.train.2.5K_2.5K.phen`
    * (5) *DatasetID* of training genotype files: `Test1.train`
    * (6) analysis window size: `80`.
    * (7) directory to store NPS data: `testdata/Test1/npsdat` (All NPS output files will be stored in this directory.)
   ```bash
   Rscript npsR/nps_init.R testdata/Test1/Test1.summstats.txt testdata/Test1 testdata/Test1/Test1.train.2.5K_2.5K.fam testdata/Test1/Test1.train.2.5K_2.5K.phen Test1.train 80 testdata/Test1/npsdat
   ```

3. **Set up a special partition for GWAS-significant SNPs.** The command argument is: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
```bash
   ./run_all_chroms.sh sge/nps_gwassig.job testdata/Test1/npsdat/
   ```

4. **Set up the decorrelated "eigenlocus" space.** This step transform the genetic data into an orthogonalized domain and one of the most time-consuming steps to run NPS. We recommend running NPS four times on shifted overlapping windows and merging the results in the last step. The recommended window shifts are 0, *WindowSize* \* 1/4, *WindowSize* \* 2/4 and WindowSize \* 3/4 SNPs, where *WindowSize* is the size of analysis window. For test set #1, set the *WindowSize* to 80, thus the window shifts should be `0`, `20`, `40` and `60`. For the default window size of `4000`, we recommend the window shifts of `0`, `1000`, `2000` and `3000`. The command arguments are: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) window shift: `0`, `20`, `40` or `60`  
   ```bash
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 0
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 20
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 40
   ./run_all_chroms.sh sge/nps_decor_prune.job testdata/Test1/npsdat/ 60
    ```
   
5. **Partition the rest of genome.** First, we define the partition cut-offs by running `npsR/nps_prep_part.R`. We recommend 10-by-10 double-partitioning on the intervals of eigenvalues of projection and estimated effect sizes in the eigenlocus space. The command arguments are:
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) Number of partitions on eigenvalues: `10`
    * (3) Number of partitions on estimated effects: `10` 
   ```
   Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 10 10 
   ```
   
   Then, we calculate partitioned polygenic risk scores in the training cohort by running `nps_part.job`. The command arguments are: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
    * (2) window shift: `0`, `20`, `40` or `60`  
   ```
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 0
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 20
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 40
   ./run_all_chroms.sh sge/nps_part.job testdata/Test1/npsdat/ 60
   ```

6. **Estimate shrinkage weights for each partition.** The command argument is: 
    * (1) NPS data directory: `testdata/Test1/npsdat`
   ```bash
   Rscript npsR/nps_reweight.R testdata/Test1/npsdat/ 
   ```

7. **Evaluate the accuracy of trained prediction model in a validation cohort.** 


Last, polygenic risk scores will be calculated for each chromosome and for each individual in the validation cohort using `sge/nps_score.dosage.job` as follows: 
   ```bash
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0 
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 20
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 40
   ./run_all_chroms.sh sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 60
   ```
   Here, the first argument for `sge/nps_score.dosage.job` is the NPS data directory (`testdata/Test1/npsdat/`), the second argument is the directory containing validation cohort data (`testdata/Test1/`), and the third argument is the DatasetID for validation genotypes. Since the genotype files for validation cohorts are named as chrom*N*.*Test1.val*.dosage.gz, DatasetID has to be `Test1.val`. The last argument is the window shift (`0`, `20`, `40` or `60`). 
   
   Finally, `npsR/nps_val.R` will combine polygenic risk scores across all shifted windows and report per-individual scores along with overall accuracy statistics: 
   ```
   Rscript npsR/nps_val.R testdata/Test1/npsdat/ Test1.val testdata/Test1/Test1.val.5K.fam testdata/Test1/Test1.val.5K.phen
   ```
   The command arguments are: 
   - NPS data directory: `testdata/Test1/npsdat/`
   - DatasetID for validation genotypes: `Test1.val`
   - sample IDs of validation cohort: `testdata/Test1/Test1.val.5K.fam`
   - phenotypes of validation samples: `testdata/Test1/Test1.val.5K.phen`
   
   `npsR/nps_val.R` will print out the following. Here, it reports the AUC of 0.8776 and Nagelkerke's R2 of 0.3172322 in the validation cohort. The polygenic risk score for each individuals in the cohort are stored in the file `testdata/Test1/Test1.val.5K.phen.nps_score`. 

   > Producing a combined prediction model...OK ( saved in testdata/Test1/Test1.val.5K.phen.nps_score )  
   > Observed-scale R2 = 0.1061336  
   > Liability-scale R2 = 0.4693835  
   > ...   
   > Data: 4729 controls < 271 cases.  
   > Area under the curve: 0.8776  
   > 95% CI: 0.8589-0.8963 (DeLong)  
   > Nagelkerke's R2 = 0.3172322  

## Running NPS on test set #2
To run test set #2 on computer clusters, see the instructions for [SGE] and [LSF]. 

## How to prepare training and validation cohorts for NPS
We take an example of UK Biobank to show how to prepare training and validation cohorts for NPS. In principle, however, NPS can work with other cohorts as far as the genotype data are prepared in the bgen file format. To gain access to UK Biobank, please check [UK Biobank data access application procedure](https://www.ukbiobank.ac.uk/). 

### Using UK Biobank as a training cohort
UK Biobank data consist of the following files: 
- **ukb_imp_chrN_v3.bgen**: imputed allelic dosages by chromosomes
- **ukb_mfi_chrN_v3.txt**: marker information by chromosomes
- **ukb31063.sample**: bgen sample information file for the entire cohort 

Assuming that UK Biobank dataset is located in `<path_to_ukbb>/`, we first exclude all Indels and SNPs with minor allele frequency < 5% or imputation quality (INFO) score < 0.4 by running the following:   
```bash
qsub -l h_vmem=4G -t 1-22 support/common_snps.job <path_to_ukbb>/ukb_imp_chr#_v3.bgen <path_to_ukbb>/ukb_mfi_chr#_v3.txt <work_dir>
```
The command arguments are: 
* file path to bgen files: `<path_to_ukbb>/ukb_imp_chr#_v3.bgen` with chromosome numbers replacing `#`
* file path to marker information files: `<path_to_ukbb>/ukb_mfi_chr#_v3.txt` with chromosome numbers replacing `#`
* work directory: `<work_dir>`, where output files will be saved.    

Next, we filter bgen files to include only training cohort samples (as specified in `<sample_id_file>`) and then export the filtered genotypes into dosage files. The `<sample_id_file>` is simply a list of sample IDs, with one sample in each line. This can be done as follows: 

```bash
qsub -l h_vmem=4G -t 1-22 support/filter_samples.job <path_to_ukbb>/ukb31063.sample <work_dir> <sample_id_file> <training_cohort_name>
```  
The command arguments are: 
* bgen sample file for the entire cohort: `<path_to_ukbb>/ukb31063.sample`
* work directory: `<work_dir>`, where output files will be saved
* file listing the sample IDs for training cohort: `<sample_id_file>`
* name of the training cohort: `<training_cohort_name>`. Name should not contain a whitespace character. Output files will be named after `<training_cohort_name>`.  

Then, we harmonize GWAS summary statitics with training cohort data: 
```bash
# CAUTION: This script uses large memory space.
Rscript support/harmonize_summstats.R <summary_statistics_file> <work_dir> <training_cohort_name>
```
The command arguments are: 
* GWAS summary statistics: `<summary_statistics_file>` in the *MINIMAL* format 
* work directory: `<work_dir>`, where output files will be saved
* name of the training cohort: `<training_cohort_name>` used in the previous step with `support/filter_samples.job`

`support/harmonize_summstats.R` will run QC filters and generate the harmonized GWAS summary statistics in the *PREFORMATTED* format (`<work_dir>/<training_cohort_name>.preformatted_summstats.txt`), which can be now used for core NPS modules. Specifically, the following QCs measures will be taken: 
* check missing values or numerical underflows in summary statistics
* remove tri-allelic SNPs
* assign reference and alternative alleles 
* remove duplicated markers 
* restrict to the markers overlapping between GWAS and training data
* cross-check allele frequencies between datasets if **effaf** is provided in the summary statistics file. Variants with too discordant allele frequencies will be rejected to prevent potential allele flips.

After that, we need to filter out SNPs in the training genotype files that were flagged for removal during the above harmonization step as follows: 
```bash
qsub -l h_vmem=4G -t 1-22 support/filter_variants.job <work_dir> <training_cohort_name>
```

Finally, the following step will create `<work_dir>/<training_cohort_name>.QC2.fam`, which keeps tracks of the IDs of all training cohort samples: 
```bash
support/make_fam.sh <work_dir> <training_cohort_name>.QC2
```
Here, we extract the sample IDs from the column header of training genotype dosage file and fill **both FID and IID** of .fam file with the same sample IDs. *If this behavior is not desirable, the .fam file has to be manually created.*

Overall, the job scripts will automatially generate the following set of NPS input files: 
- `<work_dir>/<training_cohort_name>.preformatted_summstats.txt`
- `<work_dir>/chromN.<training_cohort_name>.QC2.dosage.gz`
- `<work_dir>/<training_cohort_name>.QC2.fam`

**Note:**
* `common_snps.job` and `filter_samples.job` use bgen and qctool, respectively. The job scripts may need to be moditifed to load these modules.
* The job scripts in `support/` directory is written for SGE clusters but can be easily ported to LSF or other cluster systems.
* The job scripts use memory space up to 4GB (run with `qsub -l h_vmem=4G`). Depending on data,`support/harmonize_summstats.R` can take memory up to ~8GB. It will terminate abruptly if it runs out of memory.

### Using a different cohort as a training cohort 

To use other cohort as a training cohort, you will need to generate marker information files similar to **ukb_mfi_chrN_v3.txt** of UK Biobank. We can generate these files from bgen files (**<dataset_dir>/chromN.bgen**) as follows:
```bash
qsub -t 1-22 support/make_snp_info.job <dataset_dir>/chrom#.bgen <work_dir>
```
The first parameter (`<dataset_dir>/chrom#.bgen`) is the path to bgen genotype files, with `#` replacing a chromosome number. Internally, `support/make_snp_info.job` relies on **qctool**, thus the job script may need to be modified to load the module if needed. The output marker information files will be saved in `<work_dir>` and named as "**chrom*N*.mfi.txt**". 

The rest of steps are straight-forward and similar to using UK Biobank data: 
```bash
# Filter out InDels and SNPs with MAF < 5% or INFO < 0.4
qsub -l h_vmem=4G -t 1-22 support/common_snps.job <dataset_dir>/chrom#.bgen <work_dir>/chrom#.mfi.txt <work_dir>

# Restrict samples to those specified in <sample_id_file> and export .dosage.gz files
qsub -l h_vmem=4G -t 1-22 support/filter_samples.job <bgen_sample_file_of_entire_cohort> <work_dir> <sample_id_file> <training_cohort_name>

# Harmonize GWAS summary statistics with training genotype data 
Rscript support/harmonize_summstats.R <summary_statistics_file> <work_dir> <training_cohort_name>

# Run extra variant filtering
qsub -l h_vmem=4G -t 1-22 support/filter_variants.job <work_dir> <training_cohort_name>

# Generate .fam file with identical IIDs and FIDs
support/make_fam.sh <work_dir> <training_cohort_name>.QC2
```

### Using UK Biobank for a validation as well as a training cohort
UK Biobank can be split into two and used as a validation as well as a training cohort. Assume that `<sample_id_file>` contains the IDs of samples to include in the validation cohort. Then, validation cohort can be prepared for NPS similarly: 
```bash
# import the list of SNPs rejected while harmonizing the training cohort into the validation cohort
cp <work_dir>/<training_cohort_name>.UKBB_rejected_SNPIDs <work_dir>/<validation_cohort_name>.UKBB_rejected_SNPIDs

# Restrict samples to those specified in <sample_id_file> and export .dosage.gz files
qsub -t 1-22 support/filter_samples.job <path_to_ukbb>/ukb31063.sample <work_dir> <sample_id_file> <validation_cohort_name>

# Run extra variant filtering with <work_dir>/<validation_cohort_name>.UKBB_rejected_SNPIDs
qsub -t 1-22 support/filter_variants.job <work_dir> <validation_cohort_name>

# Generate .fam file with identical IIDs and FIDs
support/make_fam.sh <work_dir> <validation_cohort_name>.QC2
```
The `<training_cohort_name>.UKBB_rejected_SNPIDs` file contains the list of SNPs that were rejected while harmonizing the GWAS summary statistics with training cohort data. This file has to be copied to the validation cohort so that training and validation cohorts will have the same set of markers after running `support/filter_variants.job`.  

### Using a validation cohort that is independent from a training cohort
To help deploying NPS polygenic scores to a cohort that is independent from a training cohort, we provide `sge/nps_harmonize_val.job` ( `lsf/nps_harmonize_val.job` for LSF). This will convert bgen files into the dosage file format and perform the following: 
* Remove SNPs that were not defined in the trained polygenic score model. 
* Fill in SNPs that are missing in the validation cohort with 0.
* Make sure that the validation cohort genotype files and polygenic model have the consistent ordering of markers. 

This will be done by running `nps_harmonize_val.job` as follows: 
```
# Generate <work_dir>/chromN.<cohort_name>.dosage.gz files
qsub -t 1-22 -l h_vmem=4G sge/nps_harmonize_val.job <nps_data_dir> <dataset_dir>/chrom#.bgen <bgen_sample_file> <work_dir> <cohort_name>

# Generate <work_dir>/<cohort_name>.fam file with identical IIDs and FIDs
support/make_fam.sh <work_dir> <cohort_name>
```

The command arguments are: 
* NPS data directory: `<nps_data_dir>` to locate the marker information of trained polygenic risk score
* file path to bgen files of validation cohort: `<dataset_dir>/chrom#.bgen` with chromosome number replacing `#`. 
* bgen sample file: `<bgen_sample_file>` contains the sample information of cohort
* work directory: `<work_dir>`, where output files will be saved. 

**Note:**
* `nps_harmonize_val.job` relies on **qctool** internally. The job script may need to be modified to load the module.
* This script will generate harmonized genotype files named as "chrom*N*.*<cohort_name>*.dosage.gz". DatasetID to designate this files in NPS will be just *<cohort_name>*.
* Currently, `support/make_fam.sh` produces .fam files with identical FID and IID, which are extracted from .dosage.gz files. 
