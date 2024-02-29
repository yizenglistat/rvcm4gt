rename_csv_files <- function(folder_path) {
    # List all .RData files matching the pattern in the directory
    rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)
    
    # Sort files to maintain some order based on existing numbering
    rdata_files <- sort(rdata_files)

    # Rename each .RData file sequentially
    for (i in seq_along(rdata_files)) {
        new_file_name <- file.path(folder_path, sprintf("fitted_task-%d.RData", i))
        file.rename(rdata_files[i], new_file_name)
    }

    # List all .RData files matching the pattern in the directory
    rdata_files <- list.files(folder_path, pattern = "\\.RData$", full.names = TRUE)
    
    # Sort files to maintain some order based on existing numbering
    rdata_files <- sort(rdata_files)

    # Rename each .RData file sequentially
    for (i in seq_along(rdata_files)) {
        new_file_name <- file.path(folder_path, sprintf("fitted_task%d.RData", i))
        file.rename(rdata_files[i], new_file_name)
    }

    # List and remove all .h5 files
    h5_files <- list.files(folder_path, pattern = "\\.h5$", full.names = TRUE)
    sapply(h5_files, file.remove)

    # Return the count of .RData files renamed
    return(length(rdata_files))
}

combine_vectors_slash <- function(vec_a, vec_b, vec_c) {
  combined_result <- sprintf("%.3f/%.3f/%.3f", vec_a, vec_b, vec_c)
  return(combined_result)
}

combine_vectors <- function(vectorA, vectorB) {
  # Check if both vectors are of the same length
  if (length(vectorA) != length(vectorB)) {
    stop("Vectors must be of the same length.")
  }
  if(any(round(vectorA,3)==0)) vectorA[round(vectorA,3)==0] <- abs(vectorA[round(vectorA,3)==0])
  if(any(round(vectorB,3)==0)) vectorB[round(vectorB,3)==0] <- abs(vectorB[round(vectorB,3)==0])
  # Combine vectors with the specified format and round to 2 decimal places
  result <- sprintf("%.3f(%.3f)", round(vectorA, 3), round(vectorB, 3))

  return(result)
}

combine_alternate <- function(vec1, vec2, vec3, vec4, vec5) {
  # Check if all vectors are of the same length
  if (length(vec1) != length(vec2) || length(vec1) != length(vec3) ||
      length(vec1) != length(vec4) || length(vec1) != length(vec5)) {
    stop("All vectors must be of the same length.")
  }

  # Combine vectors by alternating elements
  result <- c(rbind(
  	format(vec1),
    format(vec2),
    format(round(vec3, 2), nsmall = 2),
    format(round(vec4, 2), nsmall = 2),
    format(round(vec5, 2), nsmall = 2)
  ))

  return(result)
}

safe_divide <- function(x, y) {
  # Check if both vectors are of the same length
  if (length(x) != length(y)) {
    stop("Vectors must be of the same length.")
  }
  
  # Perform division, treating 0/0 as 0
  result <- ifelse(y == 0, ifelse(x == 0, 0, Inf), x / y)
  
  return(result)
}

form_P_mat <- function(E_mat, Q_mat)
{
	nknots <- ncol(Q_mat)
	P_mat <- rbind(c(1,rep(0,nknots)), cbind(0, E_mat%*%Q_mat))
	return(P_mat)
}


misc_str_pad <- function(string) return(format(round(as.numeric(string),3), nsmall=3))

misc_str <- function(string){
      if(is.na(string)) return(NA)
      if(grepl("\\(", string)){
            slice_start <- gregexpr("\\(", string)[[1]][1] + 1
            slice_end <- nchar(string) - 1

            return(glue(misc_str_pad(substr(string, 1, slice_start-2)), 
                  "(",
                  misc_str_pad(substr(string, slice_start, slice_end)),
                  ")"))
      }else{
            return(misc_str_pad(string))
      }
}

join_tab_item <- function(x, y){

	out <- rep(NA, 17)
	out[seq(1,16,3)] <- x
	out[seq(2,17,3)] <- y
	return(out)
}

glue_tab_item <- function(x, y, digits=3)
{
	out <- unlist(lapply(glue("{round(x,digits)}({round(y,digits)})"), misc_str))

	out[y==0] <- ""

	return(out)
}

do_logit <- function(x) log(x/(1-x))
do_logit_inv <- function(x) return(1/(1+exp(-x)))

do_normalization <- function(x, log=FALSE)
{
	if(log){
		out <- (log(x)-3)
	}else{
		out <- (x-5.629612)/(69.949588-5.629612) 
		out <- out*6 - 3
	}

	return(out)
}

do_unnormalization <- function(x, x_min=5.629612, x_max=69.949588, log=FALSE)
{
	if(log){
		out <- exp(x+3)
	}else{
		out <- (x + 3)/6
		out <- (x_max-x_min)*out+x_min
	}

	return(out)
}

do_slice_add <- function(vec, p)
{
	nvec <- length(vec)
	nknots <- nvec/p
	slices <- seq(1, (p+1)*nknots, by=nknots)
	slice_base <- slices[1]:(slices[2]-1)
	base <- vec[slice_base]
	out <- c(base)
	for(idx in 2:p){
		slice <- slices[idx]:(slices[idx+1]-1)
		out <- c(out, base + vec[slice])
	}
	return(out)
}

do_twoway <- function(X){
	p <- ncol(X)
	names <- colnames(X)
	for(i in 2:(p-1)){
		for(j in (i+1):p){
			names <- c(names, glue::glue("{names[i]}_{names[j]}"))
			X <- cbind(X, as.vector(X[,i]*X[,j]))
		}
	}
	colnames(X) <- names
	return(X)
}

replace_id <- function(old, new, vec){
	vec[vec %in% old] <- new[match(vec, old, nomatch = 0)]
	return(vec)
}

remove_na <- function(df, col)
{
	ids <- df[is.na(df[,col]),c("seq_id","pool_id")]
	remove_pools <- ids$pool_id[!is.na(ids$pool_id)]
	df <- df[!(df$pool_id %in% remove_pools | df$seq_id %in% ids$seq_id),]

	return(df)
}


model_set <- function(model_name="m1", lambda=1.0){

	model_name <- substr(model_name,1,2)

	if(model_name=="m1"){
		
		f0 <- function(x, offset=0) -sin(pi*x/3) + offset
		f1 <- function(x, offset=0) x^3/16 + offset
		f2 <- function(x, offset=0.75025) -x^2/4 + offset

		f_list <- list(f0=f0, f1=f1, f2=f2)
	}

	if(model_name=="m2"){
		
		f0 <- function(x, offset=0.6392108) -0.5*exp(-sin(x)) + offset
		#f0 <- function(x, offset=0) (x/6)*cos(pi*x/3) + offset
		f1 <- function(x, offset=-0.9003) 0.3*x^2+sin(x/3) + offset
		#f2 <- function(x, offset=-0.75025) x^2/4 + offset
		f2 <- function(x, offset=-0.5) pnorm(x) + offset

		f_list <- list(f0=f0, f1=f1,f2=f2)
	
	}

	if(model_name=="m3"){
		
		f0 <- function(x, offset=-0.4729678) pnorm(x/2+cos(pi*x/2)) + offset
		f1 <- function(x, offset=0) x*cos(pi*x/3) + offset
		f2 <- function(x, offset=-0.9469494) 2^(x-1) + offset

		f_list <- list(f0=f0, f1=f1, f2=f2)
	
	}

	# x=seq(-3,3,0.001)
	# par(mfrow=c(1,3), pty='s')
	# plot(x,f0(x),'l')
	# plot(x,f1(x),'l')
	# plot(x,f2(x),'l')
	# mean(f0(x))
	# mean(f1(x))
	# mean(f2(x))

	return(f_list)

}

do_IT<-function(Y_true,Se,Sp,cj=1){
   N<-length(Y_true)
   Jmax<-N/cj 
   J<-1

   Y<-matrix(-99,nrow=N, ncol=3) 
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


   for(j in 1:(N/cj)){
	   prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
	   Z[J,1]<-rbinom(n=1,size=1,prob=prob)
	   Z[J,2]<-cj
	   Z[J,3]<-2
	   Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
	   Y[((j-1)*cj+1):(j*cj),2]<-1
	   Y[((j-1)*cj+1):(j*cj),3]<-J
	   J<-J+1
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

do_MPT<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N/cj 
   J<-1

   Y<-matrix(-99,nrow=N, ncol=3) 
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


   for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

do_DT<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N+N/cj
   J<-1

   Y<-matrix(-99,nrow=N, ncol=4)
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3)


   for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se[1],1-Sp[1])
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1   ## swab pool used
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
      if(Z[J-1,1]==1){
         for(k in ((j-1)*cj+1):(j*cj)){
            prob<-ifelse(Y_true[k]>0,Se[2],1-Sp[2])
            Z[J,1]<- rbinom(n=1,size=1,prob=prob)
            Z[J,2]<-1
            Z[J,3]<-2   ## swab ind used
            Z[J,4]<-k
            Y[k,2]<-2
            Y[k,4]<-J
            J<-J+1
         }
      }
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

do_AT<-function(Y_true, Se, Sp, cj){

	N<-length(Y_true)
	Jmax<-2*N/cj +N
	J<-1
	AT<-N/(cj^2)

	Y<-matrix(-99,nrow=N, ncol=5) 
	Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 

	Y.A<-array(-99,c(cj,cj,AT))
	ID.A<-array(-99,c(cj,cj,AT))
	ind<-1
	for(i in 1:AT){
	for(j in 1:cj){
	for(m in 1:cj){
	Y.A[m,j,i]<-Y_true[ind]
	ID.A[m,j,i]<-ind
	ind<-ind+1
	}}}



	for(s in 1:AT){

	array.yk<-Y.A[,,s]
	array.id<-ID.A[,,s]

	a<-rep(0,nrow(array.yk))
	b<-rep(0,ncol(array.yk))

	for(i in 1:cj){
	   prob<- ifelse(sum(array.yk[i,])>0, Se[1], 1-Sp[1])
	   g<- rbinom(n=1,size=1,prob=prob)
	   a[i]<-g
	   Z[J,1]<-g 
	   Z[J,2]<-cj 
	   Z[J,3]<-1
	   Z[J,4:(cj+3)]<-array.id[i,]
	   Y[array.id[i,],2]<-2
	   Y[array.id[i,],3]<-J
	   J<-J+1
	}
	for(j in 1:cj){
	   prob<- ifelse(sum(array.yk[,j])>0, Se[1], 1-Sp[1])
	   g<- rbinom(n=1,size=1,prob=prob)
	   b[j]<-g
	   Z[J,1]<-g 
	   Z[J,2]<-cj 
	   Z[J,3]<-1
	   Z[J,4:(cj+3)]<-array.id[,j]
	   Y[array.id[,j],4]<-J
	   J<-J+1
	}


	if(sum(a)>0 & sum(b)>0){
	array.yk1<-as.matrix(array.yk[(a==1),(b==1)])
	array.id1<-as.matrix(array.id[(a==1),(b==1)])
	for(i in 1:nrow(array.yk1)){
	for(j in 1:ncol(array.yk1)){
	   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
	   g<- rbinom(n=1,size=1,prob=prob)
	   Z[J,1]<-g 
	   Z[J,2]<-1 
	   Z[J,3]<-2
	   Z[J,4]<-array.id1[i,j]
	   Y[array.id1[i,j],2]<-3
	   Y[array.id1[i,j],5]<-J
	   J<-J+1
	}}}



	if(sum(a)>0 & sum(b)==0){
	array.yk1<-as.matrix(array.yk[(a==1),])
	array.id1<-as.matrix(array.id[(a==1),])
	for(i in 1:nrow(array.yk1)){
	for(j in 1:ncol(array.yk1)){
	   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
	   g<- rbinom(n=1,size=1,prob=prob)
	   Z[J,1]<-g 
	   Z[J,2]<-1 
	   Z[J,3]<-2
	   Z[J,4]<-array.id1[i,j]
	   Y[array.id1[i,j],2]<-3
	   Y[array.id1[i,j],5]<-J
	   J<-J+1
	}}}

	if(sum(a)==0 & sum(b)>0){
	array.yk1<-as.matrix(array.yk[,(b==1)])
	array.id1<-as.matrix(array.id[,(b==1)])
	for(i in 1:nrow(array.yk1)){
	for(j in 1:ncol(array.yk1)){
	   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
	   g<- rbinom(n=1,size=1,prob=prob)
	   Z[J,1]<-g 
	   Z[J,2]<-1 
	   Z[J,3]<-2
	   Z[J,4]<-array.id1[i,j]
	   Y[array.id1[i,j],2]<-3
	   Y[array.id1[i,j],5]<-J
	   J<-J+1
	}}}

	} 

	J<-J-1
	Z<-Z[1:J,]

	return(list("Z"=Z, "Y"=Y))
}

do_stats <- function(est, std, true){

	est[is.na(est)] <- 0
	std[is.na(std)] <- 0


	bias <- apply(est-true, 1, mean)
	rmse <- apply(est-true, 1, function(err) rmse(err,0))
	ssd <- apply(est, 1, sd)
	ese <- apply(std, 1, mean)

	med <- apply(est, 1, quantile, probs=c(0.5))
	lower <- apply(est, 1, quantile, probs=c(0.025))
	upper <- apply(est, 1, quantile, probs=c(1-0.025))

	return(list(bias=bias, rmse=rmse, ese=ese, ssd=ssd, med=med, lower=lower, upper=upper))
}

post_summary <- function(
	folder='./output', 
	isknown='unknown', 
	sample_size=5000, 
	model_name='m1', 
	pool_size=5, 
	testing='DT', 
	sigma=0.5,
	cache=TRUE, 
	replace=FALSE)
{
	cache_name <- glue('{folder}/{isknown}/{sample_size}/{model_name}/cj{pool_size}/{testing}_summary.RData')
	if(cache&(!replace)&file.exists(cache_name)){
		return(readRDS(cache_name))
	}else{
		nreps <- rename_csv_files(glue("{folder}/{isknown}/{sample_size}/{model_name}/cj{pool_size}/{testing}"))

		file_name <- glue('{folder}/{isknown}/{sample_size}/{model_name}/cj{pool_size}/{testing}_summary.h5')
		
		if(file.exists(file_name)) file.remove(file_name)
		
		param <- H5File$new(file_name, mode = "a")		
		
		# create store matrix
		param_med <- param$create_dataset(name="param_med", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))
		param_std <- param$create_dataset(name="param_std", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))
		param_ext <- param$create_dataset(name="param_ext", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))
		param_eci <- param$create_dataset(name="param_eci", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))

		SE <- c(0.95, 0.98)
		SP <- c(0.98, 0.99)
		#alpha <- c(-3.5, -1.0, 0.5, -0.5, 0.3, 0, 0)
		alpha	<- c(-3.3, -1.0, 0.5, -0.5, 0.5, 0, 0)
		nalpha <- length(alpha)
		delta1 <- c(1,1,1,1,1,0,0)
		delta2 <- c(1,0,1,0,1,0,0) 
		pip0 <- 1.0*I(delta1==0&delta2==0)
		pip1 <- 1.0*I(delta1==1&delta2==0)
		pip2 <- 1.0*I(delta1==1&delta2==1)

		pb <- progress::progress_bar$new(format="[:bar] :current/:total[:percent] elapsed :elapsed eta :eta",
	                       total = nreps)

		for(rep in 1:nreps){
			
			rdata_file <- glue("{folder}/{isknown}/{sample_size}/{model_name}/cj{pool_size}/{testing}/fitted_task{rep}.RData")
			fitted <- readRDS(rdata_file)

			# sigma
			param_med[1,rep] <- fitted$sigma_hat_med
			param_std[1,rep] <- fitted$sigma_hat_std
			param_eci[1,rep] <- 1.0*I(fitted$sigma_hat_lower<sigma && sigma<fitted$sigma_hat_upper)

			# SE1
			param_med[2,rep] <- fitted$se_hat_med[1]
			param_std[2,rep] <- fitted$se_hat_std[1]	
			param_eci[2,rep] <- 1.0*I(round(fitted$se_hat_lower[1],2)<=SE[1] && SE[1]<=round(fitted$se_hat_upper[1],2))

			# SE2
			param_med[3,rep] <- fitted$se_hat_med[2]
			param_std[3,rep] <- fitted$se_hat_std[2]	
			param_eci[3,rep] <- 1.0*I(fitted$se_hat_lower[2]<SE[2] && SE[2]<fitted$se_hat_upper[2])

			# SP1
			param_med[4,rep] <- fitted$sp_hat_med[1]
			param_std[4,rep] <- fitted$sp_hat_std[1]	
			param_eci[4,rep] <- 1.0*I(fitted$sp_hat_lower[1]<SP[1] && SP[1]<fitted$sp_hat_upper[1])

			# SP2
			param_med[5,rep] <- fitted$sp_hat_med[2]
			param_std[5,rep] <- fitted$sp_hat_std[2]	
			param_eci[5,rep] <- 1.0*I(fitted$sp_hat_lower[2]<SP[2] && SP[2]<fitted$sp_hat_upper[2])

			# alpha
			param_med[6:(5+nalpha),rep] <- fitted$alpha_hat_med
			param_std[6:(5+nalpha),rep] <- fitted$alpha_hat_std	
			param_eci[6:(5+nalpha),rep] <- 1.0*I(fitted$alpha_hat_lower<=alpha & alpha<=fitted$alpha_hat_upper)

			# pip1
			param_med[(6+nalpha):(5+2*nalpha),rep] <- fitted$pip1_hat
			param_std[(6+nalpha):(5+2*nalpha),rep] <- fitted$pip1_hat_std	
			param_eci[(6+nalpha):(5+2*nalpha),rep] <- 1.0*I(fitted$pip1_hat_lower<=pip1 & pip1<=fitted$pip1_hat_upper)

			# pip2
			param_med[(6+2*nalpha):(5+3*nalpha),rep] <- fitted$pip2_hat
			param_std[(6+2*nalpha):(5+3*nalpha),rep] <- fitted$pip2_hat_std	
			param_eci[(6+2*nalpha):(5+3*nalpha),rep] <- 1.0*I(fitted$pip2_hat_lower<=pip2 & pip2<=fitted$pip2_hat_upper)

			if(model_name %in% c("m1", "m2")){
				f_list <- model_set(model_name=model_name)
				true <- c(sigma, SE, SP, alpha, pip1, pip2)
			}

			npts <- 5 + 3*nalpha

			# beta(s) & f(s)
			N_test <- fitted$N_test
			varying <- c(1,3,5)
			nvarying <- 3
			fitted$levs_test_hat_med <- rep(seq(-3,3,length.out=N_test), nalpha)

			for(var_id in 1:nalpha){
				
				slice_start <- 1 + N_test*(var_id-1)
				slice_end <- N_test*var_id
				slice <- slice_start:slice_end
				if(var_id %in% varying){
					true <- c(true, f_list[[glue("f{(var_id-1)/2}")]](fitted$levs_test_hat_med[slice]) + rep(alpha[var_id], N_test))
				}else{
					true <- c(true, rep(0, N_test) + rep(alpha[var_id], N_test))
				}
				
				param_med[slice+npts,rep] <- fitted$beta_test_hat_med[slice] + rep(fitted$alpha_hat_med[var_id], N_test)
				param_std[slice+npts,rep] <- fitted$beta_test_hat_std[slice] 
				param_ext[var_id,rep] <- rmse(fitted$beta_test_hat_med[slice], true[slice])

			}

			param_ext[nvarying+1,rep] <- fitted$time
			param_ext[nvarying+2,rep] <- fitted$ntests_hat_med
			# conditional pip(delta2 | delta1=1)
			param_ext[(nvarying+3):(nvarying+2+nalpha),rep] <- safe_divide(fitted$delta2_hat, fitted$delta1_hat)
			
			pb$tick()

		}

		results <- do_stats(param_med[,], param_std[,], true)

		results$levs <- c(rep(NA, npts), fitted$levs_test_hat_med)

		results$true <- true
		results$param_med <- param_med[,]
		results$param_std <- param_std[,]
		results$param_eci <- param_eci[,]
		results$param_ext <- param_ext[,]

		param$close_all()
		
		saveRDS(results, cache_name)

		return(results)
	}
}

rmse <- function(y_hat, y_true) sqrt(mean((y_hat-y_true)^2))
















