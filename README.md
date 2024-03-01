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

- `task_id`
> The machine id. For example, 1,...,100 if running on the cluster. In this way, we will run 5 simulations independently on 100 nodes to have a total of 500 repetitions. 

- `nreps`
> The repetitions.

- `Ns`
> A vector of sample sizes.

- `pool_sizes`
> A vector of pool sizes.

- `model_names`
> A vector of model names. Different model names corresponds to different varying function sets.

- `testings`
> A vector of testing protocols such as AT (array testing), DT (Dorfman Testing) or IT (Individual Testing).

- `N_test`
> Number of knots values in inference for estimated varying functions. 

- `sigma`
> True random effect standard deviation


```r
# A demo example to run 500 repetitions in one machine.
task_id 		<- 1 						
nreps 			<- 500
Ns 				<- c(3000, 5000)
pool_sizes 		<- c(5, 10)
model_names 	<- c("m1", "m2")
testings 		<- c("AT", "DT", "IT")
N_test 			<- 600
sigma 			<- 0.5
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

