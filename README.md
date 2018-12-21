
# Non-Parametric Shrinkage (NPS)
NPS is a non-parametric polygenic risk prediction algorithm described in [Chun et al. (2018) BioRxiv](https://www.biorxiv.org/content/early/2018/07/16/370064). 


## How to install

1. The core NPS module is implemented in R. R can be downloaded from [here](https://www.r-project.org/)(v3.0 or later is required). Although NPS can run on a plain-vanilla version of R, we strongly recommend to use R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [R open](https://mran.microsoft.com/open). 

2. (Optional) NPS relies on R modules, **pROC** and **DescTools**, to calculate AUC and Nagelkerke's R2 statistics. These modules are optional; if they are not installed, AUC and Nagelkerke's R2 will not be reported. To enable this feature, please install these packages by running the following lines on command line: 

```
$ Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
$ Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
```

3. Download and unpack NPS package as below. Some of NPS codes are optimized in C++ and need to be compiled. For this, you need GNU C++ compiler with C++0x support (GCC 4.4 or later). This step will create two executable binaries, **stdgt** and **grs**, in the top-level folder. stdgt converts an imputed dosage file to standardized genotypes with the mean of 0 and variance of 1. grs calculates genetic risk scores for an imputed dosage file using NPS-generated per-SNP weighted effects. 

```
$ tar -zxvf nps-1.0.0.tar.gz
$ cd nps-1.0.0/
$ make
```

4. We recommend to run NPS with computer clusters, processing chromosomes in parallel. We provide sample job scripts in **sge** and **lsf** directories to use for SGE and LSF clusters. You may still need to modify the provided job scripts, e.g. to load necessary modules. Note that some steps are memory-intensive or take long computing time; you may need to request resources to your job scheduler. 

5. (Optional) We provide SGE job scripts to prepare UK Biobank data as NPS training genotypes. To use these scripts, you need [bgenix](https://bitbucket.org/gavinband/bgen/wiki/bgenix) and [QCTOOL v2](https://www.well.ox.ac.uk/~gav/qctool/). 

## Input files for NPS
To run NPS, you need the following set of input files: 
1. GWAS summary statistics 
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
2. Training genotypes in [QCTOOL dosage format](https://www.well.ox.ac.uk/~gav/qctool/)
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
3. Training sample information in [PLINK .fam format](https://www.cog-genomics.org/plink2/formats#fam)
```
trainF2 trainI2 0 0 0 -9
trainF3 trainI3 0 0 0 -9
trainF39 trainI39 0 0 0 -9
trainF41 trainI41 0 0 0 -9
trainF58 trainI58 0 0 0 -9
...
```

4. Training phenotypes [in PLINK phenotype format](http://zzz.bwh.harvard.edu/plink/data.shtml#pheno)
```
FID	IID	TotalLiability	Outcome
trainF68266	trainI68266	2.00037013482937	1
trainF77481	trainI77481	2.19146817328717	1
trainF39184	trainI39184	0.833227631267031	0
trainF76746	trainI76746	0.168893429255002	0
trainF65453	trainI65453	0.817881328612934	0
...
```
5. Validation genotypes in QCTOOL dosage format
6. Validation sample IDs in PLINK .fam format
7. Validation phenotypes in PLINK phenotype format


## Test cases

```
$ cd nps-1.0.0/testdata/
$ tar -zxvf NPS.Test1.tar.gz 
# This will create test data in nps-1.0.0/testdata/Test1
$ tar -zxvf NPS.Test2.tar.gz 
# This will create test data in nps-1.0.0/testdata/Test2
```

### Running NPS on Test1

1. Standardize genotypes 

```
$ cd nps-1.0.0/

# SGE cluster
$ qsub -cwd -t 1-22 sge/nps_stdgt.job testdata/Test1 Test1.train 5000

# Batch processing
$ ./batch_all_chroms.sh sge/nps_stdgt.job testdata/Test1 Test1.train 5000
```

```
$ ./nps_check.sh stdgt testdata/Test1 Test1.train 
Verifying nps_stdgt:
Checking testdata/Test1/chrom1.Test2.train ...OK
Checking testdata/Test1/chrom2.Test2.train ...OK
Checking testdata/Test1/chrom3.Test2.train ...OK
...
```

2. Set up NPS

```
$ Rscript npsR/nps_init.R testdata/Test1/Test1.summstats.txt testdata/Test1 testdata/Test1/Test1.train.2.5K_2.5K.fam testdata/Test1/Test1.train.2.5K_2.5K.phen Test1.train 80 testdata/Test1/npsdat

$ ./nps_check.sh init testdata/Test1/npsdat/
```

3. Project data to the decorrelated "eigenlocus" space
```
qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 0 
qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 20 
qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 40 
qsub -cwd -t 1-22 sge/nps_decor.job testdata/Test1/npsdat/ 60 

```
