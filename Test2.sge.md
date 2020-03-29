## Running NPS on test set #2 using SGE or UGER scheduler
With most datasets, 4GB memory space per a task will be sufficient. On SGE and UGER, the memory requirement can be specified by `qsub -l h_vmem=4G`. For nps_decor_prune.job, the hard running time limit may need to be extended, for example, with `-l h_rt=06:00:00` argument if necessary. 

```bash
cd nps-1.1.1/

# Standardize genotypes
qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test2/ Test2.train

# Configure
# CAUTION: This step needs at least 4G memory
# You need to run this as a job or open an interative session with "-l h_vmem=4G"
Rscript npsR/nps_init.R --gwas testdata/Test2/Test2.summstats.txt \
    --train-dir testdata/Test2 \
    --train-dataset Test2.train \
    --out testdata/Test2/npsdat

# Check the results
./nps_check.sh testdata/Test2/npsdat/ 

# Separate GWAS-significant SNPs
qsub -cwd -t 1-22 -l h_vmem=4G sge/nps_gwassig.job testdata/Test2/npsdat/

# Check the results
./nps_check.sh testdata/Test2/npsdat/ 

# Set up the eigenlocus space
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 0
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 1000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 2000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Partition the rest of genetic variations
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 10 10

# Calculate partitioned risk scores in the training cohort
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 0
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 1000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 2000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Estimate per-partition shrinkage weights
Rscript npsR/nps_reweight.R testdata/Test2/npsdat/

# Calculate polygenic scores for each chromosome and for each individual in the validation cohort
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0 
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 1000   
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 2000   
qsub -cwd -t 1-22 sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 3000   

# Check the results 
./nps_check.sh testdata/Test2/npsdat/ 

# Calculate overall polygenic scores and report prediction accuracies
Rscript npsR/nps_val.R --out testdata/Test2/npsdat --val-dataset Test2.val 
```
