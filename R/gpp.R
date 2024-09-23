# nmem must be divisible by nkeep
mcmc_options <- function(
	task_id=NULL, model_name='m1', outdir='output/', 
	nchain=1, nburn=500, nkeep=500, nmem=500, nthin=2, ndisp=100, nknots=100,
	a_se=0.5, b_se=0.5, a_sp=0.5, b_sp=0.5,
	a_sigma2=5, b_sigma2=50, known=FALSE, 
	phi_sd=0.1, kappa=2, hpd95=FALSE, dirac=TRUE, cdc=FALSE, 
	delete=TRUE)
{

	if(nmem>nkeep|!(nkeep%%nmem)) nmem <- nkeep
	
	return(list(task_id=task_id, model_name=model_name, outdir=outdir, 
		nchain=nchain, nburn=nburn, nkeep=nkeep, nmem=nmem, nthin=nthin, ndisp=ndisp, nknots=nknots, 
		known=known, a_se=a_se, b_se=b_se, a_sp=a_sp, b_sp=b_sp, a_sigma2=a_sigma2, b_sigma2=b_sigma2,
		phi_sd=phi_sd, kappa=kappa, hpd95=hpd95, dirac=dirac, cdc=cdc,
		delete=delete))
}

gpp_mcmc <- function(chain_id=1, data, options=mcmc_options()){

	# data
	name 	<- data$name
	testing <- data$testing
	u 		<- data$u
	X 		<- data$X
	X_lin 	<- data$X
	u_seq 	<- data$u_seq
	race 	<- data$race
	S 		<- data$S
	Z 		<- data$Z
	Y		<- data$Y
	pools 	<- data$pools
	npools 	<- data$npools
	nsites 	<- data$nsites
	nassay 	<- length(unique(Z[,3]))
	nbeta 	<- ncol(X)
	N 		<- data$N
	N_test 	<- data$N_test
	pool_size <- data$pool_size
	
	SE 		<- data$se
	SP 		<- data$sp

	# mcmc options
	nburn 	<- options$nburn
	nkeep 	<- options$nkeep
	nthin 	<- options$nthin
	nknots 	<- options$nknots
	ndisp 	<- options$ndisp
	nmem 	<- options$nmem
	known 	<- if(testing %in% c("DT", "AT")) options$known else TRUE
	phi_sd 	<- options$phi_sd
	kappa 	<- options$kappa
	dirac 	<- options$dirac
	cdc 	<- options$cdc

	a_se 	<- options$a_se
	b_se 	<- options$b_se
	a_sp 	<- options$a_sp
	b_sp 	<- options$b_sp
	a_sigma2 <- options$a_sigma2
	b_sigma2 <- options$b_sigma2

	task_id <- options$task_id
	outdir 	<- options$outdir
	model_name <- options$model_name

	# config gpp
	config 			<- gpp_config(t=u, 					# t sequence
								  t_new=u_seq,
							 	  nknots=nknots, 		# number of selected knots
							 	  nbeta=nbeta,			# number of beta including intercept
							 	  phi_sd=phi_sd,		# hyperparameter phi_sd, default is 0.1
							 	  kappa=kappa)			# hyperparameter kappa, default is 2

	nnew 			<- config$nnew 					# *note* number of t_new = number of t_unique 

	t_unique 		<- config$t_unique 				# unique t values in t sequence
	nunique 		<- config$nunique 				# number of unique t values
	t_knots  		<- config$t_knots 				# initial selected knots
	nknots 			<- config$nknots				# number of selected knots
	
	# in order to calculate correlation matrix of GPP easily, compute distance matrix in advance.
	unique_dist 	<- config$unique_dist			# nunique x nunique, L1 dist matrix
	knots_dist  	<- config$knots_dist			# nknot x nknot, L1 dist matrix
	cross_dist  	<- config$cross_dist 			# nunique x nknot, L1 dist matrix
	cross_dist_new 	<- config$cross_dist_new 		# nnew x nknot, L1 dist matrix

	# ---------------------- beta_config matrix description ----------------------- #
			# each columns means a beta function											#
			# row1: gamma distribution shape, prior distribution for tau 					#
			# row2: gamma distribution scale, prior distribution for tau 					#
			# row3: sampled tau from gamma prior distribution based on row1 and row2 		#
			# row4: uniform distribution lower bound for phi 								#
			# row5: uniform distribution upper bound for phi 								#
			# row6: phi averaged value based on uniform distribution based on row4 and row5	#
			# row7: proposed phi for M-H sampler, which made a transformation on row6 		#
			# row8: hyperparameter phi_sd, default 0.1										#
			# row9: hyperparameter kappa, default 2											#
			# ----------------------------------------------------------------------------- #

	G_list <- list()
	for(d in 1:nbeta) G_list[[d]] <- cbind(X[,d],Diagonal(x=X[,d]))

	beta_config 	<- config$beta_config 			# matrix of beta configuration for gpp 

	R_knots_list 	<- config$R_knots_list 			# list of correlation matrix, nknots x nknots
	R_knots_inv_list<- config$R_knots_inv_list		# list of inverse correlation matrix, nknots x nknots
	R_cross_list 	<- config$R_cross_list 			# list of cross correlation matrix, nunique x nknots
	Q_list 			<- config$Q_list 				# convert beta(t_knots) vector back to beta(t_unique) vector
	P_list 			<- config$P_list
	delta_list 		<- config$delta_list

	beta			<- matrix(0, N, nbeta)				# initial beta(t) matrix
	beta_unique 	<- matrix(0, nunique, nbeta)		# initial beta(t_unique) matrix
	beta_knots 	 	<- matrix(0, nknots, nbeta)			# initial beta(t_knots) matrix

	# initialization
	sigma <- runif(1)
	gamma <- rnorm(nsites,0,sigma)

	if(testing=='IT'){
		Y[,1] <- Z[,1]
		
	}else if(testing=='MPT'){
		for(idx in 1:npools){
	  		Y[seq((idx-1)*pool_size+1,idx*pool_size),1] <- Z[idx,1]
		}

	}else if(testing=='DT'){
		if(model_name!="app"){
			Y[(Y[,2]==1),1] <- 0
			Y[(Y[,2]==2),1] <- Z[Y[(Y[,2]==2),4]]
		}
	}else if(testing=='AT'){
		Y[Y[,2]==2,1]	<- 0
		Y[(Y[,2]==3),1]	<- Z[Y[(Y[,2]==3),5]]
	}
	
	omega <- rpg(N, 1, 0)
	y_latent <- Y[,1]
	h <- (y_latent - 0.5)/omega
	rands <- 0
	nalpha <- if(is.null(X_lin)) 1 else ncol(X_lin)
	alpha <- rep(0, nalpha)
	lins <- if(is.null(X_lin)) rep(0, N) else X_lin %*% alpha

	# spike and slab initialization
	delta1 <- rep(1, nalpha)
	delta2 <- rep(1, nbeta)
	large_vars <- rep(50/4, nalpha)

	# storage
	if(!dir.exists(outdir)) dir.create(outdir, recursive=TRUE, showWarnings = FALSE)

	# create h5 chain to save matrices	
	file_name <- paste0(outdir, name, 'chain', chain_id, '_task', task_id, '.h5')
	if(file.exists(file_name)) file.remove(file_name)
	chain <- H5File$new(file_name, mode = "a")	

	# create y_store
	y_store <- chain$create_dataset(name="y_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep,N), maxdims=c(nkeep,N)),chunk_dims=c(nkeep,1))
	y_store_tmp <- matrix(NA, nmem, N)

	# create beta_store
	beta_store <- chain$create_dataset(name="beta_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep, N*nbeta), maxdims=c(nkeep,N*nbeta)),chunk_dims=c(nkeep,1))
	beta_store_tmp <- matrix(NA, nmem, N*nbeta)

	# create beta_test_store
	beta_test_store <- chain$create_dataset(name="beta_test_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep, N_test*nbeta), maxdims=c(nkeep,N_test*nbeta)),chunk_dims=c(nkeep,1))
	beta_test_store_tmp <- matrix(NA, nmem, N_test*nbeta)

	# create alpha_beta_store
	alpha_beta_store <- chain$create_dataset(name="alpha_beta_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep, N_test*nbeta), maxdims=c(nkeep,N_test*nbeta)),chunk_dims=c(nkeep,1))
	alpha_beta_store_tmp <- matrix(NA, nmem, N_test*nbeta)

	# create gamma_store
	gamma_store <- chain$create_dataset(name="gamma_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nsites), maxdims=c(nkeep,nsites)),chunk_dims=c(nkeep,1))
	gamma_store_tmp <- matrix(NA, nmem, nsites)

	# create alpha_store
	alpha_store <- chain$create_dataset(name="alpha_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nalpha), maxdims=c(nkeep,nalpha)),chunk_dims=c(nkeep,1))
	alpha_store_tmp <- matrix(NA, nmem, nalpha)

	# create se and sp store
	se_store <- chain$create_dataset(name="se_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nassay), maxdims=c(nkeep,nassay)),chunk_dims=c(nkeep,1))
	se_store_tmp <- matrix(NA, nmem, nassay)

	sp_store <- chain$create_dataset(name="sp_store", dtype=h5types$H5T_NATIVE_FLOAT,
		space=H5S$new(dims=c(nkeep,nassay), maxdims=c(nkeep,nassay)),chunk_dims=c(nkeep,1))
	sp_store_tmp <- matrix(NA, nmem, nassay)

	# create sigma_store
	sigma_store <- chain$create_dataset(name="sigma_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep,1), maxdims=c(nkeep,1)),chunk_dims=c(nkeep,1))
	sigma_store_tmp <- matrix(NA, nmem, 1)

	# create delta1_store
	delta1_store <- chain$create_dataset(name="delta1_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep,nalpha), maxdims=c(nkeep,nalpha)),chunk_dims=c(nkeep,1))
	delta1_store_tmp <- matrix(NA, nmem, nalpha)

	# create delta2_store
	delta2_store <- chain$create_dataset(name="delta2_store", dtype=h5types$H5T_NATIVE_FLOAT, 
		space=H5S$new(dims=c(nkeep, nbeta), maxdims=c(nkeep, nbeta)),chunk_dims=c(nkeep,1))
	delta2_store_tmp <- matrix(NA, nmem, nbeta)

	# ----- summary statistics store required ------ 

	# create prevalence store
	prevalence_store <- chain$create_dataset(name="prevalence_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,1), maxdims=c(nkeep,1)),chunk_dims=c(nkeep,1))
	prevalence_store_tmp <- matrix(0, nmem, 1)

	# create number of tests stores
	ntests_store <- chain$create_dataset(name="ntests_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,1), maxdims=c(nkeep,1)),chunk_dims=c(nkeep,1))
	ntests_store_tmp <- matrix(0, nmem, 1)
	
	# create posterior inclusion probability for varying effect
	pip0_store <- chain$create_dataset(name="pip0_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,nbeta), maxdims=c(nkeep,nbeta)),chunk_dims=c(nkeep,1))
	pip0_store_tmp <- matrix(0, nmem, nbeta)

	# create posterior inclusion probability for fixed effect
	pip1_store <- chain$create_dataset(name="pip1_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,nalpha), maxdims=c(nkeep,nalpha)),chunk_dims=c(nkeep,1))
	pip1_store_tmp <- matrix(0, nmem, nalpha)

	# create posterior inclusion probability for varying effect
	pip2_store <- chain$create_dataset(name="pip2_store", dtype=h5types$H5T_NATIVE_FLOAT, 
	space=H5S$new(dims=c(nkeep,nbeta), maxdims=c(nkeep,nbeta)),chunk_dims=c(nkeep,1))
	pip2_store_tmp <- matrix(0, nmem, nbeta)

	niter 	<- nburn + nkeep * nthin

	# mcmc
	ikeep <- 0
	slice_end <- 0
	sample_state <- 'burn in'
		
	do_loglikeli_d <- function(case, C_knots_d, C_knots_inv_d=NULL, GdP_d=NULL,omega=NULL, h_delta_d=NULL, large_var=50/4)
	{

		nknots <- dim(C_knots_d)[1]
		# (delta_1d, delta_2d) = (0, 0)
		if(case==1) out <- 0.5 * (log(large_var) + determinant(C_knots_d, logarithm=TRUE)$modulus)
		# (delta_1d, delta_2d) = (1, 0)
		if(case==2){
			Sigma_delta_d <- sum(GdP_d^2*omega)+1/large_var
			Mu_delta_d <- (1/Sigma_delta_d)*sum(GdP_d*omega*h_delta_d)
			out <- -0.5*determinant(C_knots_inv_d, logarithm=TRUE)$modulus - 0.5*log(Sigma_delta_d) + 0.5*Sigma_delta_d%*%Mu_delta_d^2
		}

		# (delta_1d, delta_2d) = (1, 0)
		if(case==3){
			Sigma_delta_d <- t(GdP_d*omega)%*%GdP_d + rbind(c(1/large_var,rep(0,nknots)), cbind(0, C_knots_inv_d))
			Mu_delta_d <- solve(Sigma_delta_d)%*%t(GdP_d*omega)%*%h_delta_d
			out <- -0.5*determinant(Sigma_delta_d, logarithm=TRUE)$modulus + 0.5*t(Mu_delta_d)%*%Sigma_delta_d%*%Mu_delta_d
		}

		return(as.numeric(out))
	}

	running_avg1 <- c()
	running_avg2 <- c()
	running_up <- c()
	flag <- TRUE
	pb <- progress::progress_bar$new(format="[:bar] :current/:total[:percent] elapsed :elapsed eta :eta",
	                       total = niter)

	for(iter in 1:niter){

		# update Se, Sp
		if(!known){
			sesp_config	<- matrix(0,nrow=nassay,ncol=4)										# col1, col2 updated a_Se, b_Se; col3, col4 updated a_Sp, b_Sp					
			sesp_output <- update_sesp(N, npools, Y, Z, sesp_config, nassay) 				# update error prior
			updated_se 	<- matrix(sesp_output[,1:2],nassay, 2)								# updated part of Se prior
			updated_sp 	<- matrix(sesp_output[,3:4],nassay, 2)								# updated part of Sp prior
			se 			<- rbeta(nassay, a_se+updated_se[,1], b_se+updated_se[,2])
			sp 			<- rbeta(nassay, a_sp+updated_sp[,1], b_sp+updated_sp[,2])
		}else{
			se 			<- if(testing=="MPT") SE[1] else if(testing=="IT") SE[2] else SE
			sp 			<- if(testing=="MPT") SP[1] else if(testing=="IT") SP[2] else SP
		}

		if(iter>nburn&iter%%nthin==0&flag){
			running_avg1 <- c()
			running_avg2 <- c()
			running_up <- c()
			flag <- FALSE
		}

		running_avg1 <- rbind(running_avg1, delta1)
		running_avg2 <- rbind(running_avg2, delta2)
		running_up <- rbind(running_up, se)
		avg1 <- round(apply(running_avg1, 2, mean),2)
		avg2 <- round(apply(running_avg2, 2, mean),2)
		up1 <- apply(running_up, 2, quantile, 0.975, type=1)
		med1 <- apply(running_up, 2, quantile, 0.5, type=1)


		print(glue(
			"

			iter:{iter} | sigma:{round(sigma,2)}
			alpha1:{round(alpha[1],2)} | alpha2:{round(alpha[2],2)}
			alpha3:{round(alpha[3],2)} | alpha4:{round(alpha[4],2)}
			alpha5:{round(alpha[5],2)} | alpha6:{round(alpha[6],2)}
			
			delta1: {delta1[1]}, {delta1[2]}, {delta1[3]}, {delta1[4]}, {delta1[5]}, {delta1[6]}, {delta1[7]}, {delta1[8]}, {delta1[9]}
			avg1  : {avg1[1]}, {avg1[2]}, {avg1[3]}, {avg1[4]}, {avg1[5]}, {avg1[6]}, {avg1[7]}, {avg1[8]}, {avg1[9]}
			delta2: {delta2[1]}, {delta2[2]}, {delta2[3]}, {delta2[4]}, {delta2[5]}, {delta2[6]}, {delta2[7]}, {delta2[8]}, {delta2[9]}
			avg2  : {avg2[1]}, {avg2[2]}, {avg2[3]}, {avg2[4]}, {avg2[5]}, {avg2[6]}, {avg2[7]}, {avg2[8]}, {avg2[9]}

			se1:{round(se[1],4)} | se2:{round(se[2],4)} | se3:{round(se[3],4)}
			sp1:{round(sp[1],4)} | sp2:{round(sp[2],4)} | sp3:{round(sp[3],4)}

			se1_upr:{round(up1[1],4)} | se2_upr:{round(up1[2],4)} | se3_upr:{round(up1[3],4)}
			se1_med:{round(med1[1],4)} | se2_med:{round(med1[2],4)} | se3_med:{round(med1[3],4)}"))

		# update auxiliary h (or omega) 
		fixs 			<- do_rowsums(X, beta)
		rands 			<- sapply(S, function(idx) gamma[idx])
		omega 			<- rpg(N, 1, fixs + rands + lins)
		h 				<- (y_latent - 0.5)/omega 	

		# update gamma
		for(site in 1:nsites){
			gamma_sd 	<- sqrt((sigma^2)/(1+(sigma^2)*sum(omega[(S==site)]))) 
			gamma_mu 	<- (gamma_sd^2)*sum(omega[(S==site)]*(h[S==site]-fixs[S==site]-lins[S==site]))
			gamma[site] <- rnorm(1,gamma_mu,gamma_sd)
		}
		rands 			<- sapply(S, function(idx) gamma[idx])

		# update sigma 
		sigma 			<- sqrt(1/rgamma(1,a_sigma2+nsites/2, b_sigma2+sum(gamma^2)/2))

		# update y_latent
		prob  		<- do_logit_inv(fixs+rands+lins) 					
		uniform_var <- runif(N)
		newY 		<- rep(0,N)
		
		Y[,1] 		<- update_y(N, prob, Y, Z, newY, uniform_var, se, sp)
		y_latent 	<- Y[,1]
	
		# --------- alpha - spike & slab
		# update theta1
		theta1 <- rbeta(nalpha, delta1 + 1, 2 - delta1)
		# update theta2
		theta2 <- rbeta(nbeta, delta1*delta2 + 1, 1 + delta1*(1 - delta2))

		# update delta1 and delta2
		nmixed <- sample(1:nbeta,nbeta,replace=FALSE)
		for(d in nmixed){
			if(d!=1&dirac&iter>0){
				G_d 			<- G_list[[d]]
				P_d 			<- config$E_mat%*%Q_list[[d]]
				C_knots_d 		<- R_knots_list[[d]]/beta_config[3,d]
				C_knots_inv_d 	<- R_knots_inv_list[[d]]*beta_config[3,d]

				# (0, 0)
				delta_00_d 	<- delta_list[[1]]
				loglikeli_00 <- do_loglikeli_d(case=1, 
					C_knots_d=C_knots_d, 
					C_knots_inv_d=C_knots_inv_d, 
					GdP_d=0,
					omega=omega,
					h_delta_d=0,
					large_var=large_vars[d])

				# (1, 0)
				delta_10_d 	<- delta_list[[2]]
				h_10_delta 	<- h - fixs - lins - rands + X[,d]*alpha[d] 
				loglikeli_10 <- do_loglikeli_d(case=2, 
					C_knots_d=C_knots_d, 
					C_knots_inv_d=C_knots_inv_d, 
					GdP_d=X[,d],
					omega=omega,
					h_delta_d=h_10_delta,
					large_var=large_vars[d])

				# (1, 1)
				delta_11_d	<- delta_list[[3]]
				GdP_11_d 	<- cbind(X[,d], X[,d]*P_d)
				h_11_delta 	<- h - fixs - lins - rands + X[,d]*beta[,d] + X[,d]*alpha[d] 
				loglikeli_11 <- do_loglikeli_d(case=3, 
					C_knots_d=C_knots_d, 
					C_knots_inv_d=C_knots_inv_d, 
					GdP_d=GdP_11_d,
					omega=omega,
					h_delta_d=h_11_delta,
					large_var=large_vars[d])

				p00 <- 1/(1+theta1[d]*(1-theta2[d])/(1-theta1[d])*exp(loglikeli_10-loglikeli_00) + theta1[d]*theta2[d]/(1-theta1[d])*exp(loglikeli_11-loglikeli_00))
				p10 <- 1/((1-theta1[d])/(theta1[d]*(1-theta2[d]))*exp(loglikeli_00-loglikeli_10)+1+theta2[d]/(1-theta2[d])*exp(loglikeli_11-loglikeli_10))
				p11 <- 1/((1-theta1[d])/(theta1[d]*theta2[d])*exp(loglikeli_00-loglikeli_11) + (1-theta2[d])/theta2[d]*exp(loglikeli_10-loglikeli_11) + 1)
				
				draw <- rcat(1, c(p00, p10, p11), c(0,1,2))
			}else{
				
				if(cdc){
					#"Intercept ", "Race ", "New Partner", "Multiple Partners", "Contact ", "Symptoms ", "Cervical Friability ", "Cervicitis ", "PID "
					if(d==1) draw <- 2 # intercept
					if(d==2) draw <- 2 # race
					if(d==3) draw <- 1 # new
					if(d==4) draw <- 1 # multiple
					if(d==5) draw <- 1 # contact
					if(d==6) draw <- 0 # symptom
					if(d==7) draw <- 1 # cf
					if(d==8) draw <- 0 # cerv
					if(d==9) draw <- 0 # pid
				}else{
					draw <- 2
				}

			}
		
			if(draw==2){

				delta1[d] <- 1
				delta2[d] <- 1

				# update beta_d
				Xd				<- matrix(X[,d], N, 1)												
				fixs_d 			<- do_rowsums(matrix(X[,-d],N,nbeta-1), matrix(beta[,-d],N,nbeta-1))

				h_beta			<- h - fixs_d - rands - lins									# construct h_beta = h - h without beta_d
				E_mat 			<- config$E_mat 												# beta_d indicator mat, from t_unique to t seq
				Q_mat 			<- Q_list[[d]] 													# beta_d transform mat, from t_knots to t_unique
				R_knots_mat 	<- R_knots_list[[d]] 											# beta_d knots correlation matrix 
				R_knots_inv_mat <- R_knots_inv_list[[d]]										# beta_d inverse knots correlation matrix
				R_cross_mat 	<- R_cross_list[[d]]											# beta_d nunique x nknots cross corr matrix
				beta_d_config   <- beta_config[,d]												# beta_d config matrix

				# update beta_d fun and pour ouput into a list
				beta_d_output 	<- update_beta(center=TRUE, 
												Xd=Xd, omega=omega, h=h, h_beta=h_beta, 
												E_mat=E_mat, Q_mat=Q_mat, 
												knots_dist=knots_dist, cross_dist=cross_dist, 
												R_knots_mat=R_knots_mat, 
												R_knots_inv_mat=R_knots_inv_mat, 
												R_cross_mat=R_cross_mat,
												beta_d_config=beta_d_config)

				beta[,d] 				<- beta_d_output$beta_d									# update beta_d(t)
				beta_unique[,d] 		<- beta_d_output$beta_d_unique							# update beta_d(t_unique)
				beta_knots[,d] 			<- beta_d_output$beta_d_knots 							# update beta_d(t_knots)
				beta_config[3,d] 		<- beta_d_output$tau 									# update tau in beta_config (d_th)
				beta_config[6,d] 		<- beta_d_output$phi 									# update phi in beta_config (d_th)
				beta_config[7,d] 		<- beta_d_output$trphi									# update trphi in beta_config (d_th)
				R_knots_list[[d]] 		<- beta_d_output$R_knots_mat 							# update knots corr mat (d_th)
				R_knots_inv_list[[d]] 	<- beta_d_output$R_knots_inv_mat						# update inverse knots corr mat (d_th)
				R_cross_list[[d]] 		<- beta_d_output$R_cross_mat							# update nunique x nknots cross corr mat (d_th)
				Q_list[[d]] 			<- beta_d_output$Q_mat	
				fixs <- do_rowsums(X, beta)

				# update alpha_d
				large_vars[d] <- 50#1/rgamma(1, 2+1/2, 0.08+alpha[d]^2/2)
				Sigma_alpha <- solve(t(Xd <- matrix(X[,d], N, 1))%*%Diagonal(x=omega)%*%Xd+1/large_vars[d])
				mu_alpha 	<- Sigma_alpha%*%t(Xd <- matrix(X[,d], N, 1))%*%Diagonal(x=omega)%*%(h-fixs-rands-lins+Xd*alpha[d])
				alpha[d] <- as.vector(mvtnorm::rmvnorm(1,mu_alpha, as.matrix(Sigma_alpha), method = "svd"))
				lins 		<- X_lin%*%alpha 

			}else if(draw==1){
				delta1[d] <- 1
				delta2[d] <- 0

				# update beta_d
				beta_unique[,d] <- 0
				beta_knots[,d] <- 0
				beta[,d] <- 0
				fixs <- do_rowsums(X, beta)

				# update alpha_d
				large_vars[d] <- 50#1/rgamma(1, 2+1/2, 0.08+alpha[d]^2/2)
				Xd <- matrix(X[,d], N, 1)
				Sigma_alpha <- solve(t(Xd)%*%Diagonal(x=omega)%*%Xd+1/large_vars[d])
				mu_alpha 	<- Sigma_alpha%*%t(Xd)%*%Diagonal(x=omega)%*%(h-fixs-rands-lins+Xd*alpha[d])
				alpha[d] 	<- as.vector(mvtnorm::rmvnorm(1,mu_alpha, as.matrix(Sigma_alpha), method = "svd"))
				lins 		<- X_lin%*%alpha 

			}else{
				delta1[d] <- 0
				delta2[d] <- 0

				# update beta_d
				beta_unique[,d] <- 0
				beta_knots[,d] <- 0
				beta[,d] <- 0
				fixs <- do_rowsums(X, beta)

				# update alpha_d
				alpha[d] <- 0
				lins <- X_lin%*%alpha
			}
		}
		# display
		if(iter%%ndisp==0&ndisp>0){
			if(iter > nburn){
				sample_state = 'sampling'
			}
			
			verbose = paste0('iteration ',iter, ' (', sample_state,')')
			print(verbose)
		}



		# save
		if(iter > nburn){
			if(iter %% nthin==0){
				ikeep <- ikeep + 1

				beta_test <- c()
				for (d in 1:nbeta){											
					if(delta2[d]==1){
						R_knots_mat_new 				<- matern(cross_dist_new, phi=beta_config[6,d], kappa=beta_config[9,d])
						Q_mat_new 						<- R_knots_mat_new%*%R_knots_inv_list[[d]]
						Sigma_mat 						<- rep(1,nnew)-diag(tcrossprod(Q_mat_new,R_knots_mat_new))
						Sigma_mat[Sigma_mat<=0] 		<- 0
						beta_test 						<- c(beta_test, as.vector(rnorm(nnew, Q_mat_new%*%beta_knots[,d], Sigma_mat/beta_config[3,d])))
					}else{
						beta_test 						<- c(beta_test, rep(0, N_test))
					}
				}

				y_store_tmp[ikeep, ] 			<- y_latent
				beta_store_tmp[ikeep, ]			<- as.vector(beta)
				beta_test_store_tmp[ikeep, ]	<- beta_test
				alpha_beta_store_tmp[ikeep, ]	<- rep(alpha, each=N_test)+beta_test
				gamma_store_tmp[ikeep, ]		<- gamma
				alpha_store_tmp[ikeep, ]		<- alpha
				sigma_store_tmp[ikeep, ] 		<- sigma
				prevalence_store_tmp[ikeep, ]	<- mean(y_latent)
				ntests_store_tmp[ikeep, ]		<- nrow(Z)

				se_store_tmp[ikeep, ] 			<- se
				sp_store_tmp[ikeep, ] 			<- sp

				delta1_store_tmp[ikeep, ]		<- delta1
				delta2_store_tmp[ikeep, ]		<- delta2

				pip0_store_tmp[ikeep, ] 		<- 1.0*I(delta1==0&delta2==0)
				pip1_store_tmp[ikeep, ] 		<- 1.0*I(delta1==1&delta2==0)
				pip2_store_tmp[ikeep, ] 		<- 1.0*I(delta1==1&delta2==1)
				
			}
			if(ikeep==nmem){
				slice_start <- slice_end + 1
				slice_end <- slice_start + nmem - 1
				slice <- slice_start:slice_end
				print(paste0('Storing task ', options$task_id, ' chain ', chain_id, ': ', round(100*slice_end/nkeep,2),'%'))

				# hdf5r
				y_store[slice, ] 			<- y_store_tmp
				gamma_store[slice, ] 		<- gamma_store_tmp
				alpha_store[slice, ]		<- alpha_store_tmp
				sigma_store[slice, ] 		<- sigma_store_tmp
				prevalence_store[slice, ]	<- prevalence_store_tmp
				ntests_store[slice, ]		<- ntests_store_tmp 
				beta_store[slice, ]			<- beta_store_tmp
				beta_test_store[slice, ]	<- beta_test_store_tmp
				alpha_beta_store[slice, ]	<- alpha_beta_store_tmp
				se_store[slice, ] 			<- se_store_tmp
				sp_store[slice, ] 			<- sp_store_tmp

				delta1_store[slice, ]		<- delta1_store_tmp
				delta2_store[slice, ]		<- delta2_store_tmp 
				pip0_store[slice, ]			<- pip0_store_tmp
				pip1_store[slice, ]			<- pip1_store_tmp
				pip2_store[slice, ]			<- pip2_store_tmp

				ikeep <- 0
			}
		}
		pb$tick()
	}
	chain$close_all()
}

post_mean <- function(data_name=NULL, options=mcmc_options(), param_name, nparam, binary=0){
	
	draws <- array(NA,c(options$nkeep, nparam, options$nchain))

	for(chain_id in 1:options$nchain){
		file_name <- paste0(options$outdir, data_name, 'chain', chain_id, '_task', options$task_id, '.h5')
		chain <- H5File$new(file_name, mode="r+")
		#draws[,,chain_id] <- h5read(file_name, paste0(param_name,'_store'))
		draws[,,chain_id] <- chain[[paste0(param_name,'_store')]][,]
		if(binary>0){
			for(binary_idx in 1:binary){
				zero <- 2*binary_idx-1
				ones <- 2*binary_idx
				difference <- draws[,ones,chain_id] - draws[,zero,chain_id]
				draws[,zero,chain_id] <- difference
				draws[,ones,chain_id] <- difference
			}
		}
		chain$close_all()
	}


	post_mean <- apply(draws, 2, mean)
	post_med <- apply(draws, 2, median)
	post_std <- apply(draws, 2, sd)

	if(options$hpd95){
		post_bands <- apply(draws, 2, function(x) HPDinterval(as.mcmc(x), prob=0.95))
		post_lower <- post_bands[1,]
		post_upper <- post_bands[2,]
	}else{
		post_lower <- apply(draws, 2, quantile, 0.025, type=6)
		post_upper <- apply(draws, 2, quantile, 0.975, type=6)
	}
	
	return(list(mean=post_mean,median=post_med,std=post_std,lower=post_lower,upper=post_upper))
}

gpp_estimate <- function(data, options=mcmc_options()){
	
	tic <- proc.time()

	# paralell computing
	for(chain_id in 1:options$nchain){
		gpp_mcmc(chain_id, data, options)
	}

	toc <- as.vector(proc.time() - tic)[3]

	print(glue('TASK [{options$task_id}]; estimation time [min]: {round(toc/60,2)}'))

	post_y_hat <- post_mean(data$name, options, 'y', data$N)
	y_hat <- post_y_hat$mean
	y_hat_med <- post_y_hat$median
	y_hat_std <- post_y_hat$std
	y_hat_lower <- post_y_hat$lower
	y_hat_upper <- post_y_hat$upper

	post_beta_hat <- post_mean(data$name, options, 'beta', data$N*ncol(data$X))
	beta_hat <- post_beta_hat$mean
	beta_hat_med <- post_beta_hat$median
	beta_hat_std <- post_beta_hat$std
	beta_hat_lower <- post_beta_hat$lower
	beta_hat_upper <- post_beta_hat$upper

	post_beta_test_hat <- post_mean(data$name, options, 'beta_test', data$N_test*ncol(data$X))
	beta_test_hat <- post_beta_test_hat$mean
	beta_test_hat_med <- post_beta_test_hat$median
	beta_test_hat_std <- post_beta_test_hat$std
	beta_test_hat_lower <- post_beta_test_hat$lower
	beta_test_hat_upper <- post_beta_test_hat$upper

	post_alpha_beta_hat <- post_mean(data$name, options, 'alpha_beta', data$N_test*ncol(data$X))
	alpha_beta_hat <- post_alpha_beta_hat$mean
	alpha_beta_hat_med <- post_alpha_beta_hat$median
	alpha_beta_hat_std <- post_alpha_beta_hat$std
	alpha_beta_hat_lower <- post_alpha_beta_hat$lower
	alpha_beta_hat_upper <- post_alpha_beta_hat$upper

	post_se_hat <- post_mean(data$name, options, 'se', length(unique(data$Z[,3])))
	se_hat <- post_se_hat$mean
	se_hat_med <- post_se_hat$median
	se_hat_std <- post_se_hat$std
	se_hat_lower <- post_se_hat$lower
	se_hat_upper <- post_se_hat$upper

	post_sp_hat <- post_mean(data$name, options, 'sp', length(unique(data$Z[,3])))
	sp_hat <- post_sp_hat$mean
	sp_hat_med <- post_sp_hat$median
	sp_hat_std <- post_sp_hat$std
	sp_hat_lower <- post_sp_hat$lower
	sp_hat_upper <- post_sp_hat$upper

	post_sigma_hat <- post_mean(data$name, options, 'sigma', 1)
	sigma_hat <- post_sigma_hat$mean
	sigma_hat_med <- post_sigma_hat$median
	sigma_hat_std <- post_sigma_hat$std
	sigma_hat_lower <- post_sigma_hat$lower
	sigma_hat_upper <- post_sigma_hat$upper

	post_alpha_hat <- post_mean(data$name, options, 'alpha', ncol(data$X))
	alpha_hat <- post_alpha_hat$mean
	alpha_hat_med <- post_alpha_hat$median
	alpha_hat_std <- post_alpha_hat$std
	alpha_hat_lower <- post_alpha_hat$lower
	alpha_hat_upper <- post_alpha_hat$upper

	post_gamma_hat <- post_mean(data$name, options, 'gamma', data$nsites)
	gamma_hat <- post_gamma_hat$mean
	gamma_hat_med <- post_gamma_hat$median
	gamma_hat_std <- post_gamma_hat$std
	gamma_hat_lower <- post_gamma_hat$lower
	gamma_hat_upper <- post_gamma_hat$upper

	post_delta1_hat <- post_mean(data$name, options, 'delta1', ncol(data$X))
	delta1_hat <- post_delta1_hat$mean
	delta1_hat_med <- post_delta1_hat$median
	delta1_hat_std <- post_delta1_hat$std
	delta1_hat_lower <- post_delta1_hat$lower
	delta1_hat_upper <- post_delta1_hat$upper

	post_delta2_hat <- post_mean(data$name, options, 'delta2', ncol(data$X))
	delta2_hat <- post_delta2_hat$mean
	delta2_hat_med <- post_delta2_hat$median
	delta2_hat_std <- post_delta2_hat$std
	delta2_hat_lower <- post_delta2_hat$lower
	delta2_hat_upper <- post_delta2_hat$upper

	post_pip0_hat <- post_mean(data$name, options, 'pip0', ncol(data$X))
	pip0_hat <- post_pip0_hat$mean
	pip0_hat_med <- post_pip0_hat$median
	pip0_hat_std <- post_pip0_hat$std
	pip0_hat_lower <- post_pip0_hat$lower
	pip0_hat_upper <- post_pip0_hat$upper

	post_pip1_hat <- post_mean(data$name, options, 'pip1', ncol(data$X))
	pip1_hat <- post_pip1_hat$mean
	pip1_hat_med <- post_pip1_hat$median
	pip1_hat_std <- post_pip1_hat$std
	pip1_hat_lower <- post_pip1_hat$lower
	pip1_hat_upper <- post_pip1_hat$upper

	post_pip2_hat <- post_mean(data$name, options, 'pip2', ncol(data$X))
	pip2_hat <- post_pip2_hat$mean
	pip2_hat_med <- post_pip2_hat$median
	pip2_hat_std <- post_pip2_hat$std
	pip2_hat_lower <- post_pip2_hat$lower
	pip2_hat_upper <- post_pip2_hat$upper

	post_prevalence_hat <- post_mean(data$name, options, 'prevalence', 1)
	prevalence_hat <- post_prevalence_hat$mean
	prevalence_hat_med <- post_prevalence_hat$median
	prevalence_hat_std <- post_prevalence_hat$std
	prevalence_hat_lower <- post_prevalence_hat$lower
	prevalence_hat_upper <- post_prevalence_hat$upper

	post_ntests_hat <- post_mean(data$name, options, 'ntests', 1)
	ntests_hat <- post_ntests_hat$mean
	ntests_hat_med <- post_ntests_hat$median
	ntests_hat_std <- post_ntests_hat$std
	ntests_hat_lower <- post_ntests_hat$lower
	ntests_hat_upper <- post_ntests_hat$upper

	if(options$delete){
		for(chain_id in 1:options$nchain){
			file_name <- paste0(options$outdir, data$name, 'chain', chain_id, '_task', options$task_id, '.h5')
			file.remove(file_name)
		}
	}

	return(list(time=toc/options$nchain,N_test=data$N_test,
		y_hat=y_hat, y_hat_med=y_hat_med, y_hat_std=y_hat_std, y_hat_lower=y_hat_lower, y_hat_upper=y_hat_upper,
		sigma_hat=sigma_hat, sigma_hat_med=sigma_hat_med, sigma_hat_std=sigma_hat_std, sigma_hat_lower=sigma_hat_lower, sigma_hat_upper=sigma_hat_upper,
		gamma_hat=gamma_hat, gamma_hat_med=gamma_hat_med, gamma_hat_std=gamma_hat_std, gamma_hat_lower=gamma_hat_lower, gamma_hat_upper=gamma_hat_upper,
		alpha_hat=alpha_hat, alpha_hat_med=alpha_hat_med, alpha_hat_std=alpha_hat_std, alpha_hat_lower=alpha_hat_lower, alpha_hat_upper=alpha_hat_upper,
		beta_test_hat=beta_test_hat, beta_test_hat_med=beta_test_hat_med, beta_test_hat_std=beta_test_hat_std, beta_test_hat_lower=beta_test_hat_lower, beta_test_hat_upper=beta_test_hat_upper,
		alpha_beta_hat=alpha_beta_hat, alpha_beta_hat_med=alpha_beta_hat_med, alpha_beta_hat_std=alpha_beta_hat_std, alpha_beta_hat_lower=alpha_beta_hat_lower, alpha_beta_hat_upper=alpha_beta_hat_upper,
		beta_hat=beta_hat, beta_hat_med=beta_hat_med, beta_hat_std=beta_hat_std, beta_hat_lower=beta_hat_lower, beta_hat_upper=beta_hat_upper,
		se_hat=se_hat, se_hat_med=se_hat_med, se_hat_std=se_hat_std, se_hat_lower=se_hat_lower, se_hat_upper=se_hat_upper,
		sp_hat=sp_hat, sp_hat_med=sp_hat_med, sp_hat_std=sp_hat_std, sp_hat_lower=sp_hat_lower, sp_hat_upper=sp_hat_upper,
		prevalence_hat=prevalence_hat, prevalence_hat_med=prevalence_hat_med, prevalence_hat_std=prevalence_hat_std, prevalence_hat_lower=prevalence_hat_lower, prevalence_hat_upper=prevalence_hat_upper,
		ntests_hat=ntests_hat, ntests_hat_med=ntests_hat_med, ntests_hat_std=ntests_hat_std, ntests_hat_lower=ntests_hat_lower, ntests_hat_upper=ntests_hat_upper,
		delta1_hat=delta1_hat, delta1_hat_med=delta1_hat_med, delta1_hat_std=delta1_hat_std, delta1_hat_lower=delta1_hat_lower, delta1_hat_upper=delta1_hat_upper,
		delta2_hat=delta2_hat, delta2_hat_med=delta2_hat_med, delta2_hat_std=delta2_hat_std, delta2_hat_lower=delta2_hat_lower, delta2_hat_upper=delta2_hat_upper,
		pip0_hat=pip0_hat, pip0_hat_med=pip0_hat_med, pip0_hat_std=pip0_hat_std, pip0_hat_lower=pip0_hat_lower, pip0_hat_upper=pip0_hat_upper,
		pip1_hat=pip1_hat, pip1_hat_med=pip1_hat_med, pip1_hat_std=pip1_hat_std, pip1_hat_lower=pip1_hat_lower, pip1_hat_upper=pip1_hat_upper,
		pip2_hat=pip2_hat, pip2_hat_med=pip2_hat_med, pip2_hat_std=pip2_hat_std, pip2_hat_lower=pip2_hat_lower, pip2_hat_upper=pip2_hat_upper
		))
}