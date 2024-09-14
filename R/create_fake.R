# Fake Chlamydia Group Testing Generation

# 	age.txt 			store the fake chlamydia subjects age. N x 1 vector
# 	covariates.txt: 	store the fake chlamydia subjects along with attributes. N x p matrix
# 	S.txt: 				store the fake clinic locations for each subject. N x 1 vector
# 	alpha.txt: 			store the age-independent effect regression coefficients. p x 1 vector  
# 	beta.txt: 			store the age-varying effect regression coefficients. N x p matrix
# 	gamma.txt: 			store the random effect for each site

# Output: Simulated_Chlamydia_GT_Data.csv and formated data

create_fake <- function(folder="./data/fake"){

	# load settings for generating fake data
	u <- scan(glue("{folder}/age.txt"), quiet=TRUE)
	X <- as.matrix(read.table(glue("{folder}/covariates.txt"), header = TRUE))
	S <- scan(glue("{folder}/sites.txt"), quiet=TRUE)
	nsites <- length(unique(S))
	alpha <- scan(glue("{folder}/alpha.txt"), quiet=TRUE)
	beta <- as.matrix(read.table(glue("{folder}/beta.txt"), header = FALSE))
	gamma <- scan(glue("{folder}/gamma.txt"), quiet=TRUE)

	# generate 13862 individual-level status
	X[,-1] <- apply(X[,-1], 2, function(x) (x-mean(x))/sd(x))
	u <- do_normalization(u) # min=6, max=70, round to integer already for simplicity
	N_test <- 600
	u_seq <- seq(-3, 3, length.out=N_test)

	N <- nrow(X)
	p <- ncol(X)

	fixs <- X%*%alpha
	varys <- do_rowsums(X, beta)
	rands <- sapply(S, function(idx) gamma[idx])

	logits 	<- fixs + varys + rands
	probs 	<- do_logit_inv(logits)
	Y_true 	<- rbinom(N, 1, probs)

	print(mean(Y_true))

	# generate group testing responses

	# 	- Total N1 = 9288 Swab specimens
	# 		
	# 	- Total N2 = 4574 Urine specimens

	N1 <- 9288
	N2 <- N - N1

	Y_true_swab <- Y_true[1:N1]
	Y_true_urine <- Y_true[(N1+1):N]

	se_swap <- c(0.91, 0.99)
	sp_swap <- c(0.99, 0.98)
	pool_size <- 4
	res_DT 	<- do_DT(Y_true_swab, se_swap, sp_swap, pool_size)
	Z_DT 	<- res_DT$Z
	Y_DT 	<- res_DT$Y

	Z_DT[Z_DT[,3]==1,3] <- 3
	Z_DT[Z_DT[,3]==2,3] <- 1

	se_urine <- 0.98
	sp_urine <- 0.99

	res_IT 	<- do_IT(Y_true_urine, se_urine, sp_urine, 1)
	Z_IT 	<- res_IT$Z
	Z_IT[,4] <- Z_IT[,4] + length(Y_true_swab)
	Z_IT <- cbind(Z_IT, matrix(-99, nrow(Z_IT), 3))

	Y_IT 	<- res_IT$Y
	Y_IT[,3] <- Y_IT[,3] + nrow(Z_DT)
	Y_IT <- cbind(Y_IT, matrix(-99, nrow(Y_IT), 1))

	Z <- rbind(Z_DT, Z_IT)
	Y <- rbind(Y_DT, Y_IT)

	data <- list(testing="DT", 
				u=u, X=X, X_lin=X,
				N=N, npools=nrow(Z), pools=Z[,2], 
				N_test=N_test, 
				u_seq=u_seq, 
				se=NULL, sp=NULL, 
				S=S, 
				nsites=nsites, 
				Z=Z, 
				Y=Y
	)


	# mimic the structure to Iowa chlamydia group testing data
	X <- as.matrix(read.table(glue("{folder}/covariates.txt"), header = TRUE))
	X <- apply(X, 2, function(x) ifelse(x>0,"Y","N"))
	X[,1] <- scan(glue("{folder}/age.txt"), quiet=TRUE)
	X[,2] <- ifelse(X[,2]=="Y", "W", "O")
	colnames(X) <- c("Age", "Race", "Risk.New", "Risk.Mutiple", "Risk.Contact", "Symptoms", 
		"Sign.Cervical.friability", "Sign.Cervicitis", "Sign.PID")
	specimens <- c(rep("Swab", N1), rep("Urine", N2))
	df <- data.frame(
		Study.ID=1:N,
		Clinic.ID=S, 
		Pool.ID=c(Y[1:N1,3], rep(NA, N2)), 
		Specimen=specimens, 
		X,
		Ind.C.Result=c(rep(NA, N1), Z_IT[,1]),
		Pool.C.Result=NA)

	for(idx in 1:N1){
		if(-99 %in% Y_DT[idx,3:4]){
			pool_id = Y_DT[idx, 3]
			df$Pool.C.Result[idx] = Z_DT[pool_id, 1]
		}else if(!(-99 %in% Y_DT[idx,3:4])){
			pool_id = Y_DT[idx,3]
			df$Pool.C.Result[idx] = Z_DT[pool_id, 1]
			pool_id2 = Y_DT[idx,4]
			df$Ind.C.Result[idx] = Z_DT[pool_id2,1]
		}
	}
	
	write.csv(df, file = glue("{folder}/Simulated_Chlamydia_GT_Data.csv"), row.names=FALSE)

	#return(data)
}







