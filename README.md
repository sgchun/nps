
# Non-Parametric Shrinkage (NPS)
NPS is a non-parametric polygenic risk prediction algorithm described in Chun et al. [BioRxiv](https://www.biorxiv.org/content/early/2018/07/16/370064). 


## How to install

1. The core NPS module is implemented in R. R can be downloaded from [here](https://www.r-project.org/). Although NPS can run on a plain-vanilla version of R, we strongly recommend to use R linked with a linear algebra acceleration library, such as [OpenBLAS](https://www.openblas.net/), [Intel Math Kernel Library (MKL)](https://software.intel.com/en-us/articles/using-intel-mkl-with-r) or [R open](https://mran.microsoft.com/open). 

2. In addition, NPS relies on R modules, **pROC** and **DescTools** to calculate AUC and Nagelkerke's R2 statistics. These modules are optional; if they are not installed, AUC and Nagelkerke's R2 will not be reported. To enable this feature, please install these packages by running the following lines on command line: 

```
$ Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
$ Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
```

3. Please download and unpack NPS package as below. Some of NPS codes are optimized in C++ and need to be compiled. For this, you need GNU C++ compiler with C++0x support (GCC 4.4 or later). This step will create two executable binaries, **stdgt** and **grs**, in the top-level folder. stdgt converts an imputed dosage file to standardized genotypes with the mean of 0 and variance of 1. grs calculates genetic risk scores for an imputed dosage file using NPS-generated per-SNP weighted effects. 

```
$ tar -zxvf nps-1.0.0.tar.gz
$ cd nps-1.0.0/
$ make
```

4. We recommend to run NPS with computer clusters, processing chromosomes in parallel. We provide sample job scripts in *sge/* and *lsf/* directories to use for SGE and LSF clusters, respectively. You may still need to modify the provided job scripts. Also note that some steps are memory-intensive or take longer computing time; you may need to specify the resource constraints for your job scheduler. 

## Input files for NPS
To run NPS, you need the following set of input files: 
1. GWAS summary statistics 
2. Training genotypes in [QCTOOL dosage format](https://www.well.ox.ac.uk/~gav/qctool/)
3. Training sample information in [PLINK .fam format](https://www.cog-genomics.org/plink2/formats#fam)
4. Training phenotypes [in PLINK phenotype format](http://zzz.bwh.harvard.edu/plink/data.shtml#pheno)
5. Validation genotypes in QCTOOL dosage format
6. Validation sample IDs in PLINK .fam format
7. Validation phenotypes in PLINK phenotype format


## Running JLIM on provided examples  

