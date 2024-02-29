syn_options <- function(N=5000, pool_size=4, dist="uniform", u_lower=-3, u_upper=3, 
	N_test=600, nsites=64, sigma=0.5, se=c(0.95,0.98), sp=c(0.98,0.99))
{

	out <- list(N=N, npools=N/pool_size, pool_size=pool_size, 
		dist=dist,
		u_lower=u_lower, u_upper=u_upper, 
		N_test=N_test, nsites=nsites, sigma=sigma,
		se=se, sp=sp)
	return(out)
}

synthetic <- function(model_name="m1", options=syn_options()){

	N 				<- options$N
	npools 			<- options$npools
	pool_size 		<- options$pool_size
	u_lower 		<- options$u_lower
	u_upper 		<- options$u_upper
	pools 			<- rep(pool_size, npools)
	N_test 			<- options$N_test
	nsites 			<- options$nsites
	sigma 			<- options$sigma
	se 				<- options$se
	sp 				<- options$sp
	dist 			<- options$dist

	f_list <- model_set(model_name)

	if(dist=="normal"){
		u 		<- rnorm(N, 0, 1.5)
		u[u>=3] <- 3
		u[u<=-3]<- -3
		u 		<- round(u, 2)
		u_seq 	<- seq(u_lower, u_upper, length.out=N_test)
	}

	if(dist=="uniform"){
		u 		<- round(runif(N, u_lower, u_upper), 2)
		u_seq 	<- seq(u_lower, u_upper, length.out=N_test)
	}

	
	
	# fixed effect
	# X 		<- cbind(1, runif(N), rnorm(N), rbinom(N,1,0.5)-1/2, rbinom(N,1,0.5)-1/2, rbinom(N,1,0.5)-1/2, rbinom(N,1,0.5)-1/2)
	# beta 	<- cbind(f_list$f0(u), 0, f_list$f1(u), 0, f_list$f2(u), 0, 0)
	# alpha 	<- c(-2.9, -0.8, -0.5, 0.5, 0.3, 0, 0)

	# fixed effect
	X		<- cbind(1, rnorm(N), rbinom(N,1,0.5)*2-1, rbinom(N,1,0.5)*2-1, rbinom(N,1,0.5)*2-1, rbinom(N,1,0.5)*2-1, rbinom(N,1,0.5)*2-1)
	beta 	<- cbind(f_list$f0(u), 0, f_list$f1(u), 0, f_list$f2(u), 0, 0)
	alpha	<- c(-3.3, -1.0, 0.5, -0.5, 0.5, 0, 0)

	# random effect
	S		<- sample(1:nsites, N, replace=TRUE)
	gamma 	<- rnorm(nsites, 0, sigma)
	rands 	<- sapply(S, function(idx) gamma[idx])
	fixs 	<- do_rowsums(X, beta) + X%*%alpha
	f 		<- fixs + rands
	prob 	<- do_logit_inv(f)
	Y_true 	<- rbinom(N, 1, prob)
	
	res_IT 	<- do_IT(Y_true, se[2], sp[2], 1)
	Z_IT 	<- res_IT$Z
	Y_IT 	<- res_IT$Y

	res_MPT <- do_MPT(Y_true, se[1], sp[1], pool_size)
	Z_MPT 	<- res_MPT$Z
	Y_MPT 	<- res_MPT$Y

	res_DT 	<- do_DT(Y_true, se, sp, pool_size)
	Z_DT 	<- res_DT$Z
	Y_DT 	<- res_DT$Y

	res_AT 	<- do_AT(Y_true, se, sp, pool_size)
	Z_AT 	<- res_AT$Z
	Y_AT 	<- res_AT$Y

	data_IT <- list(testing="IT", X=X, u=u,
		beta=beta, prevalence=mean(Y_true),
		N=N, npools=npools, pool_size=pool_size, pools=pools, N_test=N_test, 
		u_seq=u_seq, se=se, sp=sp, S=S, nsites=nsites, Z=Z_IT, Y=Y_IT)

	data_MPT <- list(testing="MPT", X=X, u=u,
		beta=beta, prevalence=mean(Y_true),
		N=N, npools=npools, pool_size=pool_size, pools=pools, N_test=N_test, 
		u_seq=u_seq, se=se, sp=sp, S=S, nsites=nsites, Z=Z_MPT, Y=Y_MPT)

	data_DT <- list(testing="DT", X=X, u=u,
		beta=beta, prevalence=mean(Y_true),
		N=N, npools=npools, pool_size=pool_size, pools=pools, N_test=N_test, 
		u_seq=u_seq, se=se, sp=sp, S=S, nsites=nsites, Z=Z_DT, Y=Y_DT)

	data_AT <- list(testing="AT", X=X, u=u,
		beta=beta, prevalence=mean(Y_true),
		N=N, npools=npools, pool_size=pool_size, pools=pools, N_test=N_test, 
		u_seq=u_seq, se=se, sp=sp, S=S, nsites=nsites, Z=Z_AT, Y=Y_AT)

	return(list(IT=data_IT, MPT=data_MPT, DT=data_DT, AT=data_AT))

}

