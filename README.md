# GWAS tutorial
Tutorial on running basic GWAS

[![Binder](https://mybinder.org/badge_logo.svg)]( https://mybinder.org/v2/gh/dpaudel/gwas_tutorial/main?urlpath=rstudio )

# Height example

```
height <- read.table("https://raw.githubusercontent.com/dpaudel/gwas_tutorial/main/data/height.txt", header=T)
head(height)
summary(lm(height~m1))
```

# Simplified GWAS script for running in-class exercise

Data for this example comes from Romay et al 2013, Genome Biology201314:R55 [DOI: 10.1186/gb-2013-14-6-r55](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2013-14-6-r55)

## Load data and confirm samples are the same

```
phenos=read.delim('https://raw.githubusercontent.com/dpaudel/gwas_tutorial/main/data/ames_data.txt', row.names=1)
genos=read.delim('https://raw.githubusercontent.com/dpaudel/gwas_tutorial/main/data/ames_chr6.numeric.txt.gz', skip=1, row.names=1)
identical(rownames(genos), rownames(phenos))
phenos$color <- as.numeric(phenos$color)
# Extract marker positions (are part of the marker names) 
positions = as.numeric(sub(names(genos), pattern="S._", repl=""))
```
## Basic GWAS helper function

data = data frame containing phenotype and any covariates

genos = genotype matrix with samples in rows and sites in columns. (Rows should match those in data)

base_model = String specifying the linear regression model that SNPs will be added to. 

```
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
```

## Run a naive model (no covariates) for kernel color

```
model1 = "color~1" # model with only intercept; no covariates
results1 = basic.gwas(data=phenos, genos=genos, base_model=model1)
```

## plot results

```
head(results1)
plot(positions, -log10(results1)) # initial manhattan plot
```
## Model with PCs

```
model2 = "color~PC1+PC2+PC3"
results2=basic.gwas(data=phenos, genos=genos, base_model = model2)
plot(positions, -log10(results2))
```

# More advanced GWAS
## Get source codes for GWAS analysis

```
source("http://zzlab.net/GAPIT/GAPIT.library.R")
source("http://zzlab.net/GAPIT/gapit_functions.txt")
```
## Import data from Zhiwu Zhang Lab

Import phenotypes:

```
myY <- read.table("http://zzlab.net/GAPIT/data/mdp_traits.txt", head = TRUE)
```

Import genotypes:

```
myGD=read.table(file="http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
```

Import genetic map:

```
myGM=read.table(file="http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
```

Run GAPIT

```
myGAPIT <- GAPIT(
Y=myY,
GD=myGD,
GM=myGM,
PCA.total=3,
model="MLM"   # can use any of "GLM" "CMLM" "FarmCPU" "Blink"
)
```
