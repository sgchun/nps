## Running NPS on test set #1 using SGE or UGER scheduler
 
All steps have to run in the top-level NPS directory (nps-1.1.1/), and the jobs should be launched with the `qsub -cwd` option. The option `-t 1-22` will run NPS jobs over all 22 chromosomes in parallel. The job scripts are located in the [nps-1.1.1/sge/ directory](https://github.com/sgchun/nps/tree/master/sge). These scripts run not only with SGE but also with UGER, LSF and Slurm schedulers. Depending on the system, it may be necessary to modify the provided job scripts to load required modules, for example:
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
   
We provide a command line tool (`nps_check.sh`) to verify the integrity of each job. We strongly recommend to check the completion of all jobs using this tool before moving onto the next step: If `nps_check.sh` detects an error, `FAIL` message will be printed. Otherwise, only `OK` messages will be reported. 

```
cd nps-1.1.1/

# Standardize genotypes
qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test1 Test1.train

# Configure
Rscript npsR/nps_init.R --gwas testdata/Test1/Test1.summstats.txt \
       --train-dir testdata/Test1 \
       --train-dataset Test1.train \
       --window-size 80 \
       --out testdata/Test1/npsdat

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Separate GWAS-significant SNPs
qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/Test1/npsdat/

# Set up the eigenlocus space 
qsub -cwd -t 1-22 sge/nps_decor_prune.job testdata/Test1/npsdat/ 0 
qsub -cwd -t 1-22 sge/nps_decor_prune.job testdata/Test1/npsdat/ 20 
qsub -cwd -t 1-22 sge/nps_decor_prune.job testdata/Test1/npsdat/ 40 
qsub -cwd -t 1-22 sge/nps_decor_prune.job testdata/Test1/npsdat/ 60 

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Partitioning the rest of genetic variations
Rscript npsR/nps_prep_part.R testdata/Test1/npsdat/ 10 10

# Calculate partitioned risk scores in the training cohort
qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 0
qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 20
qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 40
qsub -cwd -t 1-22 sge/nps_part.job testdata/Test1/npsdat/ 60

# Check the results
./nps_check.sh testdata/Test1/npsdat/

# Estimate per-partition shrinkage weights
Rscript npsR/nps_reweight.R testdata/Test1/npsdat/ 0 

# Calculate polygenic scores for each chromosome and for each individual in the validation cohort
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 0 
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 20 
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 40 
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test1/npsdat/ testdata/Test1/ Test1.val 60 

# Check the results 
./nps_check.sh testdata/Test1/npsdat/ 

# Calculate overall polygenic scores and report prediction accuracies
Rscript npsR/nps_val.R --out testdata/Test1/npsdat --val-dataset Test1.val 
```

