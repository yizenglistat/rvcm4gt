## GVCM4GT

This repository contains R codes along with simulation results for "**Bayesian varying coefficient mixed models for group testing data**". Our model is try to estimate an individual-level regression model based on
group testing data that can capture the age-varying impact on
the Chlamydia risk with selection. To relate available information, we consider

$$
\text{logit}(\text{Pr}(\widetilde Y_i=1\mid \boldsymbol x_i, u_i))=\underbrace{\sum_{d=0}^p x_{id}\psi_d(u_i)}_{\text{Age-varying Effects}} + \underbrace{\sum_{\ell=1}^L r_\ell(i)\gamma_\ell}_{\text{Random effect}} \quad\text{for }i=1,\ldots,N.
$$


To reproduce the results in the paper, we provide implementation details as follows. 

```sh
username@login001 ~$ git clone git@github.com:yizenglistat/GVCM4GT.git
username@login001 ~$ cd GVCM4GT
```

### Arguments

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

