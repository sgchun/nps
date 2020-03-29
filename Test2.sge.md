
## Running NPS on test set #2 using SGE clusters

NPS can be run on test set #2 similarly as test set #1 except: 
* Since test set #2 has a total of ~5,000,000 genome-wide common SNPs, we recommend to use the window size of **4,000 SNPs**. And accordingly, window shifts should be set to **0, 1,000, 2,000, and 3,000 SNPs**. This is the setting we generally recommend for real data as well. 
* Some of NPS steps now require large memory space. With most datasets, 4GB memory space per a task is sufficient for NPS. On SGE, the memory requirement can be specified by `qsub -l h_vmem=4G`.
* For test set #2 and real data sets, we do not recommend running NPS without parallelization because of heavy computational requirements.

```bash
cd nps-1.1.1/

# Standardize genotypes
qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test2/ Test2.train

# Configure
# CAUTION: This step requires large memory space
# You may need to run this as a job or open an interative session for it by:
# qlogin -l h_vmem=4G
Rscript npsR/nps_init.R --gwas testdata/Test2/Test2.summstats.txt \
    --train-dir testdata/Test2 \
    --train-dataset Test2.train \
    --out testdata/Test2/npsdat


qsub -cwd -t 1-22 -l h_vmem=4G sge/nps_gwassig.job testdata/Test2/npsdat/

./nps_check.sh testdata/Test2/npsdat/ 

# Set up the decorrelate eigenlocus space
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 0
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 1000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 2000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_decor_prune.job testdata/Test2/npsdat/ 3000

# Check the results
./nps_check.sh testdata/Test2/npsdat/

# Define partitioning boundaries 
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 0 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 1000 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 2000 10 10
Rscript npsR/nps_prep_part.R testdata/Test2/npsdat/ 3000 10 10

# Calculate partitioned risk scores in the training cohort
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 0
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 1000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 2000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_part.job testdata/Test2/npsdat/ 3000

# Check the results of last step
./nps_check.sh last testdata/Test2/npsdat/ 0 1000 2000 3000

# Estimate per-partition shrinkage weights
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 0
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 1000
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 2000
Rscript npsR/nps_weight.R testdata/Test2/npsdat/ 3000

# (Optional) Report the overall AUC of prediction in the training cohort
Rscript npsR/nps_plot_shrinkage.R testdata/Test2/npsdat/ Test2.nps.pdf 0 1000 2000 3000

# (Optional) Generate a plot of overall shrinkage curves
Rscript npsR/nps_train_AUC.R testdata/Test2/npsdat/ 0 1000 2000 3000

qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_back2snpeff.job testdata/Test2/npsdat/ 0
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_back2snpeff.job testdata/Test2/npsdat/ 1000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_back2snpeff.job testdata/Test2/npsdat/ 2000
qsub -cwd -l h_vmem=4G -t 1-22 sge/nps_back2snpeff.job testdata/Test2/npsdat/ 3000

# Check the results of last step
./nps_check.sh last testdata/Test2/npsdat/ 0 1000 2000 3000

# Convert back to per-SNP effect sizes
qsub -cwd -t 1-22 sge/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0
qsub -cwd -t 1-22 sge/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 1000
qsub -cwd -t 1-22 sge/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 2000
qsub -cwd -t 1-22 sge/nps_score.job testdata/Test2/npsdat/ testdata/Test2/ Test2.val 3000

# Check the results of nps_score
./nps_check.sh score testdata/Test2/npsdat/ testdata/Test2/ Test2.val 0 1000 2000 3000

# Calculate polygenic scores for each chromosome and for each individual in the validation cohort
Rscript npsR/nps_val.R testdata/Test2/npsdat/ testdata/Test2/ testdata/Test2/Test2.val.5K.fam testdata/Test2/Test2.val.5K.phen 0 1000 2000 3000
```

`nps_train_AUC.R` will report the following AUC in the training cohort:
> Data: 2500 controls < 2500 cases.  
> Area under the curve: **0.7843**  
> 95% CI: 0.7718-0.7968 (DeLong)  

`nps_plot_shrinkage.R` will plot [the curves of estimated conditional mean effects](https://github.com/sgchun/nps/blob/master/testdata/Test2.nps.pdf) and save it to a pdf file (`Test2.nps.pdf`). 

`nps_val.R` will report the following overall prediction accuracy in the validation cohort: 
> Non-Parametric Shrinkage 1.1  
> Validation cohort:  
> Total  5000 samples  
> 240  case samples  
> 4760  control samples  
> 0  samples with missing phenotype (-9)  
> Includes TotalLiability  
> Checking a prediction model (winshift = 0 )...  
> Observed-scale R2 = 0.04862955   
> Liability-scale R2 = 0.2303062   
> Checking a prediction model (winshift = 1000 )...  
> Observed-scale R2 = 0.04994584  
> Liability-scale R2 = 0.2298484  
> Checking a prediction model (winshift = 2000 )...  
> Observed-scale R2 = 0.05150205  
> Liability-scale R2 = 0.2268046  
> Checking a prediction model (winshift = 3000 )...  
> Observed-scale R2 = 0.05258402  
> Liability-scale R2 = 0.2298871  
>  
>  
>  
> Producing a combined prediction model...OK (saved in **testdata/Test2/Test2.val.5K.phen.nps_score** )  
> Observed-scale R2 = 0.05253146  
> Liability-scale R2 = 0.2376991  
> Loading required package: pROC  
> Type 'citation("pROC")' for a citation.  
>  
> Attaching package: ‘pROC’  
>  
> The following objects are masked from ‘package:stats’:  
>  
> cov, smooth, var  
>  
> AUC:  
>  
> Call:  
> roc.default(controls = prisk[vlY == 0], cases = prisk[vlY ==     1], ci = TRUE)  
>   
> Data: 4760 controls < 240 cases.  
> **Area under the curve: 0.7886**  
> 95% CI: 0.7617-0.8154 (DeLong)  
> Loading required package: DescTools  
> **Nagelkerke's R2 = 0.1668188**   
>  
> Call: &nbsp;glm(formula = vlY ~ prisk, family = binomial(link = "logit"))  
>  
> Coefficients:  
> (Intercept) &nbsp; &nbsp; &nbsp; &nbsp; &nbsp;prisk  
> &nbsp; &nbsp; -5.2888 &nbsp; &nbsp; &nbsp; 0.2387  
>  
> Degrees of Freedom: 4999 Total (i.e. Null);  4998 Residual  
> Null Deviance: &nbsp; &nbsp; &nbsp;1926  
> Residual Deviance: 1652 &nbsp; &nbsp; &nbsp; &nbsp; AIC: 1656  
> Done  

