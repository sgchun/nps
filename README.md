
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

6. **Estimate shrinkage weights for each partition.** If the sex covariate is included in the training data (.fam file), ...  The command argument is: 
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

