
# Non-Parametric Shrinkage (NPS)



## How to install

```
make
```

pROC

## How to run JLIM on provided example  

```

for chrom in `seq 1 22`; do env SGE_TASK_ID=$chrom sge/nps_stdgt.job testdata/ Test1.train 5000; done

# In SGE Cluster
qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/ Test1.train 5000

sge/nps_check.sh stdgt testdata/ Test1.train 


```

```
Rscript npsR/nps_init.R testdata/Test1.summstats.txt testdata/ testdata/Test1.train.2.5K_2.5K.fam testdata/Test1.train.2.5K_2.5K.phen Test1.train 80 testdata/npsdat

./sge/nps_check.sh init testdata/npsdat/

```
qsub -cwd -t 1-22 sge/nps_decor.job testdata/npsdat/ 0 
qsub -cwd -t 1-22 sge/nps_decor.job testdata/npsdat/ 20 
qsub -cwd -t 1-22 sge/nps_decor.job testdata/npsdat/ 40 
qsub -cwd -t 1-22 sge/nps_decor.job testdata/npsdat/ 60 

for chrom in `seq 1 22`; do env SGE_TASK_ID=$chrom sge/nps_decor.job testdata/npsdat/ 0; done 

./sge/nps_check.sh decor testdata/npsdat/ 0 

qsub -cwd -t 1-22 sge/nps_prune.job testdata/npsdat/ 0
qsub -cwd -t 1-22 sge/nps_prune.job testdata/npsdat/ 20
qsub -cwd -t 1-22 sge/nps_prune.job testdata/npsdat/ 40
qsub -cwd -t 1-22 sge/nps_prune.job testdata/npsdat/ 60

qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/npsdat/ 0
qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/npsdat/ 20
qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/npsdat/ 40
qsub -cwd -t 1-22 sge/nps_gwassig.job testdata/npsdat/ 60

Rscript npsR/nps_prep_part.R testdata/npsdat/ 0 10 10 
Rscript npsR/nps_prep_part.R testdata/npsdat/ 20 10 10 
Rscript npsR/nps_prep_part.R testdata/npsdat/ 40 10 10 
Rscript npsR/nps_prep_part.R testdata/npsdat/ 60 10 10 

for chrom in `seq 1 22`; do env SGE_TASK_ID=$chrom sge/nps_part.job testdata/npsdat/ 0; done 

Rscript npsR/nps_weight.R testdata/npsdat/ 0 
Rscript npsR/nps_weight.R testdata/npsdat/ 20 
Rscript npsR/nps_weight.R testdata/npsdat/ 40 
Rscript npsR/nps_weight.R testdata/npsdat/ 60 

./sge/nps_check.sh weight

Rscript npsR/nps_train_AUC.R testdata/npsdat/ 0 20 40 60
Data: 2500 controls < 2500 cases.
Area under the curve: 0.8799
95% CI: 0.8707-0.8891 (DeLong)

Rscript npsR/nps_plot_shrinkage.R testdata/npsdat/ Test1.nps.pdf 0 20 40 60

for chrom in `seq 1 22`; do echo $chrom; env SGE_TASK_ID=$chrom sge/nps_back2snpeff.job testdata/npsdat/ 0; done 

./sge/nps_check.sh back2snpeff testdata/npsdat/ 0 

for chrom in `seq 1 22`; do echo $chrom; env SGE_TASK_ID=$chrom sge/nps_score.job testdata/npsdat/ Test1.train testdata/ Test1.val ; done 
for chrom in `seq 1 22`; do echo $chrom; env SGE_TASK_ID=$chrom sge/nps_score.job testdata/npsdat/ Test1.train.win_20 testdata/ Test1.val ; done 
for chrom in `seq 1 22`; do echo $chrom; env SGE_TASK_ID=$chrom sge/nps_score.job testdata/npsdat/ Test1.train.win_40 testdata/ Test1.val ; done 
for chrom in `seq 1 22`; do echo $chrom; env SGE_TASK_ID=$chrom sge/nps_score.job testdata/npsdat/ Test1.train.win_60 testdata/ Test1.val ; done 

Rscript npsR/nps_val.R testdata/npsdat/ testdata/ testdata/Test1.val.5K.fam testdata/Test1.val.5K.phen 0 20 40 60 

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
### Example data

### Running JLIM on 

## File formats

## Using UK Biobank as a training cohort

