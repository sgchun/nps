
# Non-Parametric Shrinkage (NPS)
NPS is a non-parametric polygenic risk prediction algorithm, which is described in Chun et al. ([BioRxiv] https://www.biorxiv.org/content/early/2018/07/16/370064). 


## How to install

The core NPS module is implemented in R ([Download] https://www.r-project.org/). We strongly recommend to use R linked with a linear algebra acceleration library, such as [OpenBLAS] https://www.openblas.net/ or [R open] https://mran.microsoft.com/open. Altough NPS can run on a plain-vanilla version of R, BLAS accerlation can significantly reduce the running time of NPS. We do not provide expertise on how to install BLAS libraries for R. 

To calculate AUC and Nagelkerke's R2, we use R modules, **pROC** and **DescTools**. These modules are optional; if not installed, AUC and Nagelkerke's R2 calculation will be skipped. To enable this, please install the packages by the following on command line: 

```
Rscript -e 'install.packages("pROC", repos="http://cran.r-project.org")' 
Rscript -e 'install.packages("DescTools", repos="http://cran.r-project.org")' 
```

To install NPS, please unpack the provided NPS package as below. We optimized certains steps using C++ codes. To compile these codes, you need GNU C++ compiler with C++0x support. If this is successful, you will see two executable binaries: **stdgt** and **grs** in the top-level directory. stdgt converts dosage files to standardized genotypes. grs calculate risk scores based on dosage and NPS-generated per-SNP effect sizes. 

```
tar -zxvf nps-1.0.0.tar.gz
cd nps-1.0.0/
make
```

We recommend to run NPS in a compute cluster, processing chromosomes in parallel. We provide sample job scripts in *sge/* and *lsf/* directories for SGE and LSF clusters, respectively. Depending on the system, you may need to modify the provided job scripts. Also note that some steps are memory-intensive or consume longer computing time; you may need to specify the resource constraints to your job scheduler. 

## Running JLIM on provided examples  

