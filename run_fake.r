rm(list=ls(all=TRUE))

SEED <- 414 # ssvs
#SEED <- 415 # novs

set.seed(SEED)

# ..........................................................................................
# prerequisites
# ..........................................................................................
packages <- c("BayesLogit", "extraDistr", "coda", "geoR", "ltsa", "mvtnorm", "Matrix","hdf5r","Rcpp","glue")
temp <- suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE))
temp <- sapply(list.files("R", pattern="*.R$", full.names=TRUE, ignore.case=TRUE), function(x) source(x))
temp <- sapply(list.files("src", pattern="*.cpp$", full.names=TRUE, ignore.case=TRUE), function(x) sourceCpp(x))
rm(temp)

# ..........................................................................................
# settings
# ..........................................................................................
task_id <- SEED
model_name <- "fake"
N_test <- 600
dirac <- TRUE
datadir <- "./data/fake"
outdir <- glue("output/{model_name}/")
if(dirac){
	savefile <- glue("output/{model_name}/fitted_ssvs_fake.RData")
}else{
	savefile <- glue("output/{model_name}/fitted_novs_fake.RData")
}
# ..........................................................................................
# create and format simulated fake data
# ..........................................................................................
create_fake(datadir)
data <- preprocess_fake(folder=datadir, N_test=N_test)

# ..........................................................................................
# mcmc options
# ..........................................................................................
options <- mcmc_options(
	task_id=task_id,							# task id number
	model_name=model_name, 						# model name
	outdir=outdir, 								# output folder
	nchain=1, 									# number chain running
	nburn=5000, 								# number burn-in samples before keeping
	nkeep=2000, 								# number of kept samples
	nmem=1000, 									# number of memory in storing
	nknots=100, 								# number of knots for GPP
	nthin=10, 									# number of thinning
	ndisp=1, 									# display for debugging
	a_se=0.5, b_se=0.5, 						# prior for Se
	a_sp=0.5, b_sp=0.5,							# prior for Sp
	a_sigma2=2, b_sigma2=1, 					# prior for sigma2
	known=FALSE,								# Se and Sp unknown
	phi_sd=0.1, kappa=2,						# hyper-parameter for GPP
	dirac=dirac, 								# SSVS or not
	hpd95=TRUE, 								# hpd95 or equal-tail
	delete=TRUE,  								# FALSE means keep the whole chain for convergence diagnosis
)

# ..........................................................................................
# mcmc
# ..........................................................................................
fit <- gpp_estimate(data, options)
saveRDS(fit, file=savefile)

# ..........................................................................................
# create figures for fake data in the Supplementary Materials
# ..........................................................................................
verbose <- TRUE
if(dirac){
	savefile <- glue("output/results/app/fitted_ssvs_fake.RData")
}else{
	savefile <- glue("output/{model_name}/fitted_novs_fake.RData")
}
fit <- readRDS(savefile)

knots <- rep(seq(-3, 3, length.out=fit$N_test), times=length(fit$alpha_beta_hat)/fit$N_test)
knots <- do_unnormalization(knots)

mean <- fit$alpha_beta_hat
lower <- fit$alpha_beta_hat_lower
upper <- fit$alpha_beta_hat_upper
label <- rep(c("Intercept", "Race", "New", "Multiple", "Contact", "Symptoms", "Cervical", "Cervicitis", "PID"), each=fit$N_test)
df <- data.frame(knots=knots, mean=mean, lower=lower, upper=upper, label=label)
	
if(verbose){
par(mfrow=c(3,3))
	if(dirac){
		ylims <- list(c(-5.5,0), c(-1,1), c(-0.5,1), 
				c(-0.5,1), c(-0.5,2), c(-0.1,0.1),
				c(-0.5,1), c(-0.1,0.1), c(-0.1,0.1))
	}else{
		ylims <- list(c(-5.5,0), c(-1.8,1), c(-0.5,1.8), 
				c(-0.5,1.5), c(-0.5,2), c(-0.5,0.5),
				c(-0.5,1.5), c(-1.5,1), c(-1.5,1.0))
	}
	ylabs <- c(expression(psi[1](u)), expression(psi[2](u)), expression(psi[3](u)),
			expression(psi[4](u)), expression(psi[5](u)), expression(psi[6](u)),
			expression(psi[7](u)), expression(psi[8](u)), expression(psi[9](u)))
	for(idx in 1:9){
		xx <- df[df$label==unique(label)[idx], "knots"]
		yy <- df[df$label==unique(label)[idx], "mean"]
		lw <- df[df$label==unique(label)[idx], "lower"]
		up <- df[df$label==unique(label)[idx], "upper"]
		plot(xx, yy, 'l', lwd=2, lty=1, ylim=ylims[[idx]], col="black",
			xlab='u=Age', ylab=ylabs[idx], main=unique(label)[idx])
		abline(h=0, col='red', lwd=1.5, lty=2)
		lines(xx, lw, lty=1, lwd=1.0, col="darkgray")
		lines(xx, up, lty=1, lwd=1.0, col="darkgray")
		polygon(c(xx, rev(xx)), c(up, rev(lw)), 
			col=rgb(0.1, 0.1, 0.1, 0.1), border=NA)
	}
}
if(dirac){
	saveRDS(df, "output/results/app/df_ssvs_fake.RData")
}else{
	saveRDS(df, "output/results/app/df_novs_fake.RData")
}
# ..........................................................................................
# create the assay accuracy probabilities table for fake data in the Supplementary Materials
# ..........................................................................................
if(dirac){
	savefile <- glue("output/{model_name}/fitted_ssvs_fake.RData")
}else{
	savefile <- glue("output/{model_name}/fitted_novs_fake.RData")
}
fit <- readRDS(savefile)

sesp_tab <- data.frame(
	se_mean = fit$se_hat, 
	se_lower = fit$se_hat_lower, 
	se_upper = fit$se_hat_upper,
	sp_mean = fit$sp_hat, 
	sp_lower = fit$sp_hat_lower, 
	sp_upper = fit$sp_hat_upper,
	row.names=c("swab individual", "urine individual", "swab pool")
)
print(sesp_tab)
# ..........................................................................................
# end
# ..........................................................................................