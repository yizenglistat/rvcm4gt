gpp_config<-function(t, t_new, nknots=100, nbeta, phi_sd=0.1, kappa=2, step=0.01){
	N 				<- length(t)						# number of individuals
	eta_idx 		<- rep(NA, N)						# linear predictor eta initialized

	t_lower 		<- min(t)							# lower bound of variable t
	t_upper 		<- max(t)							# upper bound of variable t
	
	t_unique 		<- sort(unique(t))					# sorted extracted t_unique
	nunique 		<- length(t_unique)					# number of unique in t sequence 
	
	nnew 			<- length(t_new)					# number of unique for new t sequence
														# generate indicator matrix for tunqiue for N individuals 
	for(idx in 1:nunique){
		eta_idx[ t == t_unique[idx] ] = idx 			# store those individuals for each t_unique
	}
	E_mat = Matrix(0, N, nunique)						# initialize E_mat indicator matrix N x nunique
	E_mat[cbind(1:N, eta_idx)] = 1 						# assign those individuals to 1 if appears in eta_idx

	# distance matrix for computing correlation matrix later in convenience
	t_knots = seq(t_lower, t_upper, length.out=nknots)
	unique_dist = as.matrix(dist(t_unique, t_unique, method = "manhattan", diag = FALSE, upper = FALSE))
	knots_dist  = as.matrix(dist(t_knots, t_knots, method = "manhattan", diag = FALSE, upper = FALSE))
	cross_dist  = matrix(NA, nunique, nknots)
	cross_dist_new = matrix(NA, nnew, nknots)
	for(idx in 1:nknots){
	  cross_dist[,idx] = abs(t_unique-t_knots[idx])
	  cross_dist_new[,idx] = abs(t_new-t_knots[idx])
	}

	# configure the range of phi
	t_range = t_upper-t_lower
	phi_seq = seq(0.01,100,by=0.01)
	phi_upper = phi_lower = rep(NA,length(phi_seq))
	for(idx in 1:length(phi_seq)){
	  phi_lower[idx] = matern(t_range*2/30,phi=phi_seq[idx],kappa=kappa)
	  phi_upper[idx] = matern(t_range*2/3,phi=phi_seq[idx],kappa=kappa)
	}
	phi_min = phi_seq[which.min(abs(phi_lower-0.05))]
	phi_max = phi_seq[which.min(abs(phi_upper-0.05))]

	# configure required beta functional parameters 

			# ---------------------- beta_config matrix description ----------------------- #
			# each columns means a beta function											#
			# row1: gamma distribution shape, prior distribution for tau 					#
			# row2: gamma distribution scale, prior distribution for tau 					#
			# row3: sampled tau from gamma prior distribution based on row1 and row2 		#
			# row4: uniform distribution lower bound for phi 								#
			# row5: uniform distribution upper bound for phi 								#
			# row6: phi averaged value based on uniform distribution based on row4 and row5	#
			# row7: transformed phi for M-H sampler, which made a transformation on row6 	#
			# row8: hyperparameter phi_sd, default 0.1										#
			# row9: hyperparameter kappa, default 2											#
			# ----------------------------------------------------------------------------- #

	beta_config <- matrix(NA, 9, nbeta)
	for (d in 1:nbeta){
		beta_config[1,d] <- 2
		beta_config[2,d] <- 1
		beta_config[3,d] <- rgamma(1, beta_config[1,d], beta_config[2,d])
		beta_config[4,d] <- phi_min
		beta_config[5,d] <- phi_max
		beta_config[6,d] <- (beta_config[4,d]+beta_config[5,d])/2
		beta_config[7,d] <- do_logit((beta_config[6,d]-beta_config[4,d])/(beta_config[5,d]-beta_config[4,d]))
		beta_config[8,d] <- phi_sd
		beta_config[9,d] <- kappa
	}

	R_knots_list 	 	<- list()						# initial list knots corr mat
	R_knots_inv_list	<- list() 						# initial list inverse knots corr mat
	R_cross_list 		<- list()						# initial list nunique x nknots cross corr mat
	Q_list 				<- list()						# initial list conditional transform mat
	P_list 				<- list()						# initial list E_mat x conditional transform mat 
	for (d in 1:nbeta){
		R_knots_list[[d]]		<- toeplitz(matern(knots_dist[1,], phi= beta_config[6,d], kappa = kappa))
		R_knots_inv_list[[d]] 	<- TrenchInverse(R_knots_list[[d]])
		R_cross_list[[d]] 		<- matern(cross_dist, phi=beta_config[6,d], kappa=kappa)
		Q_list[[d]] 			<- R_cross_list[[d]]%*%R_knots_inv_list[[d]]
		P_list[[d]] 			<- form_P_mat(E_mat, Q_list[[d]])
	}

	delta_list <- list()
	delta_list[[1]] <- Diagonal(x=c(0, rep(0, N)))
	delta_list[[2]] <- Diagonal(x=c(1, rep(0, N)))
	delta_list[[3]] <- Diagonal(x=c(1, rep(1, N)))

	# output list
	output_list = list(E_mat=E_mat,kappa=kappa, nnew=nnew,
					   t_unique=t_unique, nunique=nunique, t_knots=t_knots, nknots=nknots,
					   unique_dist=unique_dist, 
					   knots_dist=knots_dist, 
					   cross_dist=cross_dist,
					   cross_dist_new=cross_dist_new, 
					   phi_init=c(phi_min,phi_max),
					   beta_config=beta_config,
					   R_knots_list=R_knots_list, 
					   R_knots_inv_list=R_knots_inv_list,
					   R_cross_list=R_cross_list,
					   Q_list=Q_list,
					   P_list=P_list,
					   delta_list=delta_list
					   )
	return(output_list)
}

update_beta <- function(center=FALSE, Xd, omega=1, sigma=1.0, epsilon=1.0,
	h, h_beta, E_mat, Q_mat, 
	knots_dist, cross_dist, 
	R_knots_mat, R_knots_inv_mat, R_cross_mat, 
	beta_d_config)
{
	
	# ------------------------------------------ Update beta_d(t), beta_d(t_unique), beta_d(t_knots) ------------------------------------------ #

	tau  			 <- beta_d_config[3]													# tau parameter for beta_d
	Omega_d2		 <- Diagonal(x=Xd^2*omega)												# Omega_d2. diagonal matrix
	EOmegah_beta  	 <- crossprod(E_mat, Xd*omega*h_beta)									# part of updated mean of GPP
	D_mat 			 <- crossprod(E_mat,Omega_d2)%*%E_mat
	presicion_mat 	 <- forceSymmetric(crossprod(Q_mat, D_mat)%*%Q_mat + R_knots_inv_mat*(tau/epsilon)*(sigma**2))
	presicion_decomp <- chol(presicion_mat)													# chol decomposition to make inverse faster
	Sigma_beta_d  	 <- chol2inv(presicion_decomp) 											# updated covariance matrix of GPP
	mu_beta_d  		 <- Sigma_beta_d %*% crossprod(Q_mat,EOmegah_beta)   					# updated mean of GPP

	# use override old matrix to make name short, otherwise, updated_beta_d_knots, etc
	nknots 		  	 <- nrow(R_knots_mat)													# number of knots
	beta_d_knots  	 <- as.numeric(backsolve(presicion_decomp/sigma,rnorm(nknots))+mu_beta_d)		# sample from presicion matrix,solve upper triangle system
	
	beta_d_unique 	 <- Q_mat%*%beta_d_knots 												# convert from knots to unique values in t
	beta_d_knots  	 <- beta_d_knots  - mean(beta_d_unique)*center							# center beta_d(t_knots) or not 
	beta_d_unique 	 <- beta_d_unique  - mean(beta_d_unique)*center							# center beta_d(t_unique) or not
	beta_d        	 <- as.vector(E_mat%*%beta_d_unique)									# updated beta_d(t)

	# -------------------------------------------------------------- Update tau -------------------------------------------------------------- #
	
	a 				 <- beta_d_config[1]													# shape parameter for gamma prior
	b 				 <- beta_d_config[2]													# scale parameter for gamma prior
	updated_a 		 <- a + nknots/2														# updated shape and scale for gamma prior
	updated_b 		 <- as.numeric(b + crossprod(beta_d_knots, R_knots_inv_mat)%*%beta_d_knots/(2*epsilon))
	updated_tau 	 <- rgamma(1, updated_a, updated_b)										# sample a new tau from updated gamma prior
	
	# -------------------------------------------------------- Update GPP correlation matrix ------------------------------------------------------- #
	
	# -$-$-$-$-$-$-$-$- configure required arguments -$-$-$-$-$-$-$-$-

	phi_min 							<- beta_d_config[4]									# phi lower bound
	phi_max 		 					<- beta_d_config[5]									# phi upper bound
	phi     						 	<- beta_d_config[6]									# phi parameter in covariance matrix
	trphi   						  	<- beta_d_config[7]									# trphi value for M-H algorithm.
	phi_sd 								<- beta_d_config[8]									# hyperparameter, phi_sd
	kappa 								<- beta_d_config[9]									# hyperparameter, kappa

	# -$-$-$-$-$-$-$-$- proposed point from simple random walk -$-$-$-$-$-$-$-$-

	proposed_trphi 						<- rnorm(1, trphi, phi_sd)										# propose a new phi from a Gaussian random walk		
	proposed_phi   						<- (exp(proposed_trphi)*phi_max+phi_min)/(1+exp(proposed_trphi)) 	# transformation on the proposed phi
	proposed_R_knots_mat 				<- toeplitz(matern(knots_dist[1,], phi=proposed_phi, kappa=kappa)) 	# updated knots mat for beta_d
	proposed_R_knots_inv_mat 			<- TrenchInverse(proposed_R_knots_mat)								# updated inverse knots mat for beta_d
	proposed_R_cross_mat 				<- matern(cross_dist, phi=proposed_phi, kappa=kappa)				# updated nunique x nknots cross mat for beta_d
	proposed_Q_mat 						<- proposed_R_cross_mat%*%proposed_R_knots_inv_mat					# updated Q_mat for beta_d

	h_minus_eta 						<- h_beta - Xd * beta_d 												# h - eta where eta is linear predictor
	proposed_h_minus_eta 				<- h_beta - Xd * (E_mat%*%proposed_Q_mat%*%beta_d_knots)				# proposed (h - eta): N x 1 vector

	proposed_one 						<-  -determinant(proposed_R_knots_mat)$modulus/2-(updated_tau/epsilon)*crossprod(beta_d_knots, proposed_R_knots_inv_mat)%*%beta_d_knots/2-sum(proposed_h_minus_eta*omega*proposed_h_minus_eta)/(2*sigma**2)
	current_one 						<-  -determinant(R_knots_mat)$modulus/2-(updated_tau/epsilon)*crossprod(beta_d_knots, R_knots_inv_mat)%*%beta_d_knots/2-sum(h_minus_eta*omega*h_minus_eta)/(2*sigma**2)

	rate 								<- (proposed_phi-phi_min)*(phi_max-proposed_phi)/((phi_max-phi)*(phi-phi_min))	# to control the accept probability								
	accept_prob 						<- min(1, as.numeric(exp(proposed_one-current_one)*rate))						# accept probability to proposed one or reject

	bernoulli_var 						<- rbinom(1,1,accept_prob)														# generate bernoulli variable with accept prob
	updated_phi  						<- bernoulli_var*proposed_phi   + (1-bernoulli_var)*phi 						# updated phi (could be no change)
	updated_trphi						<- bernoulli_var*proposed_trphi + (1-bernoulli_var)*trphi 						# updated trphi (could be no change)

	# use override old matrix to be consistent with the input matrix name
	if(bernoulli_var==1){ 																	# if proposed accepted, then update all matrics
	  R_knots_mat 						<- proposed_R_knots_mat 							# updated knots correlation matrix
	  R_knots_inv_mat 					<- proposed_R_knots_inv_mat							# updated inverse knots correlation matrix			
	  R_cross_mat 						<- proposed_R_cross_mat 							# updated (nunique x nknots) cross unique knots correlation matrix
	  Q_mat 							<- proposed_Q_mat 									# updated Q_mat converting knots to unique values in t.
	}

	output_list <- list(beta_d=beta_d, beta_d_knots=beta_d_knots, beta_d_unique=beta_d_unique, 
						tau=updated_tau, phi=updated_phi, trphi=updated_trphi,
						R_knots_mat=R_knots_mat, R_knots_inv_mat=R_knots_inv_mat, R_cross_mat=R_cross_mat,Q_mat=Q_mat)

	# output_list <- list(beta_d=beta_d, beta_d_knots=beta_d_knots, beta_d_unique=beta_d_unique, 
	# 					tau=updated_tau, phi=phi, trphi=trphi,
	# 					R_knots_mat=R_knots_mat, R_knots_inv_mat=R_knots_inv_mat, R_cross_mat=R_cross_mat,Q_mat=Q_mat)


	return(output_list)
}

