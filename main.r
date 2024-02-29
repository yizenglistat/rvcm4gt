rm(list=ls(all=TRUE))
graphics.off()

# running locally which might be time-consuming to run reps=500
setwd("set-working-directory") 
task_id <- 1
nreps <- 500
knowns <- c(FALSE)
Ns <- c(5000)
pool_sizes <- c(5, 10)
model_names <- c("m1", "m2")
testings <- c("AT", "DT", "IT")
N_test <- 600
sigma <- 0.5
folder <- 'output'

packages <- c("BayesLogit", "RcppEigen", "extraDistr", "coda", "geoR", "ltsa", "mvtnorm", "Matrix","hdf5r","Rcpp","glue", "Hmisc")
temp <- suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
temp <- sapply(list.files("src", pattern="*.cpp$", full.names=TRUE, ignore.case=TRUE), function(x) sourceCpp(x))
rm(temp)

for(rep in 1:nreps){
	for(known in knowns){
		for(N in Ns){
			for(model_name in model_names){
				for(pool_size in pool_sizes){

					#set.seed(2222)

					sim_options <- syn_options(
						N=N, 
						pool_size=pool_size, 
						dist="uniform", 
						u_lower=-3, 
						u_upper=3, 
						N_test=N_test, 
						nsites=64, 
						sigma=sigma, 
						se=c(0.95,0.98), 
						sp=c(0.98,0.99)
					)

					full_data <- synthetic(
						model_name=model_name, 
						options=sim_options
					)

					for(testing in testings){
						data <- full_data[[testing]]
						print(data$prevalence)
						options <- mcmc_options(
							task_id=nreps*(task_id-1)+rep,
							model_name=model_name, 
							outdir=glue('output/{ifelse(known, "known", "unknown")}/{N}/{model_name}/cj{pool_size}/{testing}/'), 
							nchain=1, 
							nburn=1000, 
							nkeep=2000, 
							nmem=2000, 
							nknots=100,
							nthin=5, 
							ndisp=-1, 
							a_se=1/2, b_se=1/2, 
							a_sp=1/2, b_sp=1/2,
							a_sigma2=2.0, b_sigma2=1.0, 
							known=known,
							dirac=TRUE, 
							phi_sd=0.1, kappa=2,
							delete=TRUE, 
							seed=FALSE
						)

						fit <- try(gpp_estimate(data, options), silent=TRUE)
						
						if(all(class(fit)=="try-error")){
							rdata_file <- glue(options$outdir, "fitted_task", options$task_id, ".ERROR")
							if(!file.exists(rdata_file)) file.create(rdata_file)
						}else{
							rdata_file <- glue(options$outdir, "fitted_task", options$task_id, ".RData")
							saveRDS(fit, file=rdata_file)
						}

					}
				}
			}
		}
	}
}
