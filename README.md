# GWAS tutorial
Tutorial on running basic GWAS

[![Binder](https://mybinder.org/badge_logo.svg)]( https://mybinder.org/v2/gh/dpaudel/gwas_tutorial/main?urlpath=rstudio )

## Get source codes for GWAS analysis

```
source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
```


## Import phenotype

```
myY <- read.table("http://zzlab.net/GAPIT/data/mdp_traits.txt", header = T)
```
## Import genotype

```
myG <- read.table("http://zzlab.net/GAPIT/data/mdp_genotype_test.hmp.txt", head = FALSE)
```

