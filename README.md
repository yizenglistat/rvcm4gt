## GVCM4GT

This repository contains R codes along with simulation results for "**Bayesian varying coefficient mixed models for group testing data**". To reproduce the results in the paper, we provide implementation details as follows. 

```sh
username@login001 ~$ git clone git@github.com:yizenglistat/GVCM4GT.git
username@login001 ~$ cd GVCM4GT
```

### Configuration


```r
# install required packagaes below in alphabetical order
# R >= 3.6.1
> install.packages(...)
```

- BayesLogit
- coda
- extraDistr
- geoR, [v1.8-1](https://cran.r-project.org/src/contrib/Archive/geoR/geoR_1.8-1.tar.gz)
- glue
- hdf5r
- Hmisc
- ltsa
- Matrix
- mvtnorm
- Rcpp
- RcppEigen


### Arguments

>> `task_id`

```r
# task_id will be distributed as 1,...,100 for example if running on the cluster. 
# In this way, we will run 5 simulations independently on 100 nodes to have a total of 500 repetitions. 
task_id 		<- 1 						
nreps 			<- 500
Ns 				<- c(5000)
pool_sizes 		<- c(5, 10)
model_names 	<- c("m1", "m2")
testings 		<- c("AT", "DT", "IT")
N_test 			<- 600
sigma 			<- 0.5
folder 			<- 'output'

packages <- c("BayesLogit", "RcppEigen", "extraDistr", "coda", "geoR", "ltsa", "mvtnorm", "Matrix","hdf5r","Rcpp","glue", "Hmisc")
temp <- suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
temp <- sapply(list.files("src", pattern="*.cpp$", full.names=TRUE, ignore.case=TRUE), function(x) sourceCpp(x))
rm(temp)
```

### Demo

### Thanks to:

* Carl Boettiger and his [template package](https://github.com/cboettig/template)
* Jeff Hollister and his [manuscriptPackage](https://github.com/jhollist/manuscriptPackage)
* Robert Flight: http://rmflight.github.io/posts/2014/07/analyses_as_packages.html
* Hadley Wickham: http://r-pkgs.had.co.nz/
* Yihui Xie: http://yihui.name/knitr/
* Rstudio


### Links

* https://github.com/ropensci/rrrpkg
* https://github.com/Reproducible-Science-Curriculum/rr-init
* http://ropensci.github.io/reproducibility-guide/
* https://github.com/jdblischak/r-project-workflows
* https://github.com/benmarwick/rrtools

