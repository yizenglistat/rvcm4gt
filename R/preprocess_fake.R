preprocess_fake <- function(folder="data/fake", N_test){

	csv_file <- glue("{folder}/simulated_chlamydia_gt_data.csv")
	df <- read.csv(csv_file)

	colnames(df) <- c(
		"seq_id", 
		"site_id", 
		"pool_id", 
		"specimen", # swab or urine  
		"age", # age
		"race", # nonhispanic-white = 1
		"risk_new", # new partner = 1 
		"risk_mult", # multiple partners = 1 
		"risk_contact", # contact with STD partner = 1
		"symptoms", # symptoms = 1
		"sign_cf", # Sign.Cervical.friability = 1
		"sign_c", # Sign.Cervicitis = 1
		"sign_pid", # sign PID = 1
		"ct", 
		"ct_pool")

	df$race <- ifelse(df$race=='W', 1, 0)
	df$race <- as.numeric(df$race)
	cols_Y_or_N <- c("risk_new", "risk_mult", "risk_contact", "symptoms", "sign_cf", "sign_c", "sign_pid")

	for(col in cols_Y_or_N) df[,col] <- ifelse(df[,col]=='Y', 1, 0)


	# Divide the data into Swab and Urine data sets
	# Note that Urine samples are tested individually (pool size = 1)
	df_indv <- df[is.na(df$pool_id),]
	df_pool <- df[!is.na(df$pool_id),]

	# Creating the design matrices and age
	u_indv <- df_indv[,5]
	u_pool <- df_pool[,5]

	X_indv <- cbind(1,df_indv[,6:13])
	X_pool <- cbind(1,df_pool[,6:13])

	# Building the Z and Y matrices
	Z_indv <- matrix(-99, nrow=nrow(df_indv), ncol=7)
	Y_indv <- matrix(-99, nrow=nrow(df_indv), ncol=4)

	Z_pool <- matrix(-99, nrow=nrow(df_pool), ncol=7)
	Y_pool <- matrix(-99, nrow=nrow(df_pool), ncol=4)


	# assay
	# 1: swab indv
	# 2: urine indv
	# 3: swab pool

	# For Swab testing
	pids <- unique(df_pool$pool_id)

	trackid <- 1

	for(idx in 1:length(pids)){

		temp <- df_pool[df_pool$pool_id==pids[idx],]
		temp_id <- temp$seq_id

		ct <- temp$ct

		pool_size <- length(ct)

		retesting <- ifelse(all(is.na(ct)), 0, 1)

		# pool-level
		Z_pool[trackid, 1] <- retesting
		Z_pool[trackid, 2] <- pool_size
		Z_pool[trackid, 3] <- 3		## swab pool assay used
		Z_pool[trackid, 4:(pool_size+3)] <- temp_id

		# indv-level
		Y_pool[temp_id, 1] <- -99
		Y_pool[temp_id, 2] <- 1
		Y_pool[temp_id, 3] <- trackid

		if(retesting==0){
			trackid <- trackid + 1
		}

		if(retesting>0){
			tid <- (trackid+1):(trackid+pool_size)
			Z_pool[tid,1] <- ct
			Z_pool[tid,2] <- 1
			Z_pool[tid,3] <- 1		## swab individual assay used
			Z_pool[tid,4] <- temp_id

			Y_pool[temp_id,1] <- -99
			Y_pool[temp_id,2] <- 2
			Y_pool[temp_id,4] <- tid
			trackid <- max(tid)+1
		}
	}

	Z_pool <- Z_pool[1:(trackid-1),]


	# For Urine testing
	Z_indv[,1] <- df_indv$ct 	## Test outcome for Chlamydia
	Z_indv[,2] <- 1                          ## Pool size    
	Z_indv[,3] <- 2				## urine assay used
	Z_indv[,4] <- (1:(nrow(df_indv))) + nrow(Y_pool)	

	Y_indv[,1] <- -99	## Diagnosed status
	Y_indv[,2] <- 1
	Y_indv[,3] <- 1:(nrow(df_indv)) + nrow(Z_pool)

	# Putting everything together
	S <- df$site_id
	nsites <- length(unique(S))
	u <- c(u_pool, u_indv)
	u <- do_normalization(u)
	u_seq <- seq(-3, 3, length.out=N_test)
	X <- rbind(X_pool, X_indv)
	X[,-1] <- apply(X[,-1], 2, function(x) (x-mean(x))/sd(x))
	X <- as.matrix(X)
	Z <- rbind(Z_pool, Z_indv)
	Y <- rbind(Y_pool, Y_indv)

	data <- list(testing="DT", 
				u=u, X=X, X_lin=X,
				N=nrow(Y), npools=nrow(Z), pools=Z[,2], 
				N_test=N_test, 
				u_seq=u_seq, 
				se=NULL, sp=NULL, 
				S=S, 
				nsites=nsites, 
				Z=Z, 
				Y=Y)

	return(data)
}



# if(FALSE){

# library(xgboost)

# # Define the function
# xgboost_predict <- function(x, y, xtest) {
#   # Convert x and y into a DMatrix format for xgboost
#   dtrain <- xgb.DMatrix(data = as.matrix(x), label = y)
  
#   # Define xgboost parameters
#   params <- list(
#     objective = "reg:squarederror",  # For regression tasks
#     max_depth = 6,                   # Maximum depth of a tree
#     eta = 0.1,                       # Learning rate
#     nthread = 2,                     # Number of threads to use
#     verbosity = 0                    # Silent mode
#   )
  
#   # Train the model using xgboost
#   model <- xgboost(data = dtrain, params = params, nrounds = 100, verbose = 0)
  
#   # Prepare xtest for prediction
#   dtest <- xgb.DMatrix(data = as.matrix(xtest))
  
#   # Predict ytest using the trained model
#   ytest <- predict(model, dtest)
  
#   return(ytest)
# }

# u_seq <- seq(-3,3,length.out=600)
# u <- data$u
# y <- df[df$label=="Race","med"]
# race_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="Intercept","med"]
# intercept_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="New","med"]
# New_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="Multiple","med"]
# Multiple_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="Contact","med"]
# Contact_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="Symptoms","med"]
# Symptoms_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="Cervical","med"]
# Cervical_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="Cervicitis","med"]
# Cervicitis_function <- xgboost_predict(u_seq, y, u)
# y <- df[df$label=="PID","med"]
# PID_function <- xgboost_predict(u_seq, y, u)
# #"Intercept", "Race", "New", "Multiple", "Contact", "Symptoms", "Cervical", "Cervicitis", "PID")

# #u <- data$u
# #plot(u, 0.5*race_function, ylim=c(-1,1))



# data <- preprocess(folder='./data', vars=c(1:9), linear=NULL, twoway=FALSE, N_test=600, reformat=FALSE, fake=FALSE)

# folder="./data"


# u <- data$u
# write.table(round(do_unnormalization(u),2), file = "data/fake/age.txt", row.names = FALSE, col.names = FALSE)
# u <- scan("data/fake/age.txt", quiet=TRUE)


# S <- sample(data$S)
# write.table(S, file = "data/fake/sites.txt", row.names = FALSE, col.names = FALSE)
# S <- scan("data/fake/sites.txt", quiet=TRUE)

# N <- length(S)
# X <- data$X

# intercept <- 1
# X <- cbind(intercept, X[,-1])
# write.table(X, file = "data/fake/covariates.txt", row.names = FALSE, col.names = TRUE)
# X <- as.matrix(read.table("data/fake/covariates.txt", header = TRUE))
# # X[,-1] <- apply(X[,-1], 2, function(x) (x-mean(x))/sd(x))

# beta <- matrix(0, ncol=9, nrow=13862)
# beta[,1] <- intercept_function - mean(intercept_function)
# beta[,2] <- 1.5*(race_function - mean(race_function))
# beta[,3] <- 1.2*(New_function - mean(New_function))
# beta[,4] <- 1.2*(Multiple_function - mean(Multiple_function))
# beta[,5] <- 1.0*(Contact_function - mean(Contact_function))
# beta[,6] <- 1.0*(Symptoms_function - mean(Symptoms_function))
# beta[,7] <- 1.2*(Cervical_function - mean(Cervical_function))
# beta[,8] <- Cervicitis_function - mean(Cervicitis_function)
# beta[,9] <- PID_function - mean(PID_function)
# write.table(beta, file = "data/fake/beta.txt", row.names = FALSE, col.names = FALSE)
# beta <- as.matrix(read.table("data/fake/beta.txt", header = FALSE))

# alpha	<- fitted$alpha_hat
# write.table(alpha, file = "data/fake/alpha.txt", row.names = FALSE, col.names = FALSE)
# alpha <- scan("data/fake/alpha.txt", quiet=TRUE)


# # random effect
# gamma 	<- fitted$gamma_hat
# write.table(gamma, file = "data/fake/gamma.txt", row.names = FALSE, col.names = FALSE)
# gamma <- scan("data/fake/gamma.txt", quiet=TRUE)

# rands 	<- sapply(S, function(idx) gamma[idx]) 
# logits 	<- do_rowsums(X, beta) + X%*%alpha + rands
# probs 	<- do_logit_inv(logits)
# Y_true 	<- rbinom(N, 1, probs)
# print(mean(Y_true))

# }






