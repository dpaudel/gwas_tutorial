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

# A simplified GWAS script for running our in-class exercise
# Data for this example comes from Romay et al 2013, Genome Biology201314:R55 DOI: 10.1186/gb-2013-14-6-r55

# Load data and confirm samples are the same
phenos=read.delim('ames_phenos.txt', row.names=1)
genos=read.delim('ames_color_chr6.numeric.txt.gz', skip=1, row.names=1)

phenos=read.delim('10_ames_data.txt', row.names=1)
genos=read.delim('09_ames_chr6.numeric.txt.gz', skip=1, row.names=1)
identical(rownames(genos), rownames(phenos))

# Basic GWAS helper function
#   data = data frame containing phenotype and any covariates
#   genos = genotype matrix with samples in rows and sites in columns. (Rows should match those in data)
#   base_model = String specifying the linear regression model that SNPs will be added to. 
basic.gwas=function(data, genos, base_model){
  pvals = sapply(names(genos), function(site){
	mymodel = paste(base_model, site, sep="+")
	print(mymodel)
	mydata=cbind(data, genos[,site])
	names(mydata)[ncol(mydata)] = site	# Reset name
	fit = lm(mymodel, data=mydata)
	pvals = summary(fit)$coefficients[,4] 
	return(pvals[site])	# return the p-value corresponding to the site
  })
  return(pvals)
}

# Run a naive model (no covariates) for kernel color
phenos$color <- as.numeric(phenos$color)

# GWAS model
model1 = "color~1" # model with only intercept; no covariates
results1 = basic.gwas(data=phenos, genos=genos, base_model=model1)
# plot results
plot(positions, -log10(results1)) # initial manhattan plot
#

model2 = "color~PC1+PC2+PC3"
results2=basic.gwas(data=phenos, genos=genos, base_model = model2)
plot(positions, -log10(results2))
