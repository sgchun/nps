## Running NPS on test set #2 using LSF schedulers

In most cases, 4GB memory space per a task will be sufficient for running NPS jobs. On LSF clusters, the memory requirement can be specified by `bsub -R 'rusage[mem=4000]'`. 

```bash
cd nps-1.1.1/

# Standardize genotypes
bsub -J stdgt[1-22] sge/nps_stdgt.job testdata/Test2/ Test2.train

# Configure
# CAUTION: This step needs at least 4G memory
# You need to run this as a job or open an interative session with "-R 'rusage[mem=4000]'"
Rscript npsR/nps_init.R --gwas testdata/Test2/Test2.summstats.txt \
    --train-dir testdata/Test2 \
    --train-dataset Test2.train \
    --out testdata/Test2/npsdat

# Check the results
./nps_check.sh testdata/Test2/npsdat/ 

# Separate GWAS-significant SNPs
bsub -R 'rusage[mem=4000]' -J gwassig[1-22] sge/nps_gwassig.job testdata/Test2/npsdat/

# Check the results
./nps_check.sh testdata/Test2/npsdat/ 

# Set up the eigenlocus space
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 0
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 1000
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 2000
bsub -R 'rusage[mem=4000]' -J decor[1-22] sge/nps_decor_prune.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Partition the rest of genetic variations
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 10 10

# Calculate partitioned risk scores in the training cohort
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 0
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 1000
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 2000
bsub -J part[1-22] sge/nps_part.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Estimate per-partition shrinkage weights
Rscript npsR/nps_reweight.R testdata/Test2/npsdat/

# Calculate polygenic scores for each chromosome and for each individual in the validation cohort
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0 
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 1000   
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 2000   
bsub -J score[1-22] sge/nps_score.dosage.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 3000   

# Check the results 
./nps_check.sh testdata/Test2/npsdat/ 

# Calculate overall polygenic scores and report prediction accuracies
Rscript npsR/nps_val.R --out testdata/Test2/npsdat --val-dataset Test2.val
```

