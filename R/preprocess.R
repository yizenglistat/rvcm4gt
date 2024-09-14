preprocess <- function(folder="./data", vars=NULL, linear=NULL, twoway=FALSE, N_test=600, reformat=FALSE, fake=FALSE)
{

	if(reformat){

		csv_file <- glue("{folder}/raw.csv")
		df <- read.csv(csv_file, stringsAsFactors=F, na.strings = "NaN")
		df <- df[, 3:ncol(df)]
		df$Study.ID <- as.numeric(substring(df$Study.ID,4,10))

		colnames(df) <- c(
			"seq_id", 
			"date", 
			"site_id", 
			"age", # age
			"male", # male = 1
			"race", # nonhispanic-white = 1
			"ethnicity", # H: 1900, N: 12270, U: 207 不知道是啥 hispanic or non-hispanic
			"method", # character variable:"tma" 不知道是啥
			"pool_id", 
			"specimen", # swab or urine  
			"reason", # F: 5413   I:323    P:221    R:350   S:8210 不知道是啥
			"risk_new", # new partner = 1 
			"risk_mult", # multiple partners = 1 
			"risk_contact", # contact with STD partner = 1
			"risk_msm", # msm = 1 men sex with mex 
			"risk_none", # none = 1
			"symptoms", # symptoms = 1
			"sign_cf", # Sign.Cervical.friability = 1
			"sign_c", # Sign.Cervicitis = 1
			"sign_pid", # sign PID = 1
			"sign_u", # sign Urethritis = 1
			"sign_no", # sign No Exam Performed = 1
			"sign_none", # Sign None = 1
			"insurance", 
			"ct", 
			"gc", 
			"ct_pool", 
			"gc_pool")

		filter_colnames <- c(
			"seq_id", "pool_id", "site_id", "specimen", "reason",
			"age", "ethnicity", "race", "risk_new", "risk_mult", "risk_contact", "risk_msm",
			"risk_none", "symptoms", "sign_cf", "sign_c", "sign_pid", 
			"sign_u", "sign_no", "sign_none", "ct", "ct_pool")
		df <- df[, filter_colnames]
		

		# replace site id
		df$site_id <- replace_id(sort(unique(df$site_id)), 1:length(unique(df$site_id)), df$site_id)

		# age, remove 110455 that contains 
		ids <- df[(df$age=="#VALUE!"),c("seq_id","pool_id")]
		remove_pools <- ids$pool_id[!is.na(ids$pool_id)]
		df <- df[!(df$pool_id %in% remove_pools | df$seq_id %in% ids$seq_id),]
		df$age <- as.numeric(df$age)

		ids <- df[(df$age<2)|(df$age>80),c("seq_id","pool_id")]
		#ids <- df[(df$age<10)|(df$age>60),c("seq_id","pool_id")]
		remove_pools <- ids$pool_id[!is.na(ids$pool_id)]
		df <- df[!(df$pool_id %in% remove_pools | df$seq_id %in% ids$seq_id),]

		df$ethnicity[is.na(df$ethnicity)] <- "U"

		# remove NA from race, ..., sign_none
		for(col in c("ethnicity", "race", "risk_new", "risk_mult", "risk_contact", "risk_msm",
			"risk_none", "symptoms", "sign_cf", "sign_c", "sign_pid", 
			"sign_u", "sign_no", "sign_none", "ct")){
			df <- remove_na(df, col)
		}

		# ct, remove "E" status
		ids <- df[df$ct=="E", c("seq_id","pool_id")]
		remove_pools <- ids$pool_id[!is.na(ids$pool_id)]
		df <- df[!(df$pool_id %in% remove_pools | df$seq_id %in% ids$seq_id),]

		# fix specimen error
		df[df[,1]==20677,'specimen'] <- 'Swab' # df[(!is.na(df$pool_id))&df$pool_id==110644,]

		# replace pool id
		df[!is.na(df$pool_id),]$pool_id <- replace_id(unique(df[!is.na(df$pool_id),]$pool_id), 1:length(unique(df[!is.na(df$pool_id),]$pool_id)), df[!is.na(df$pool_id),]$pool_id)

		#df$race <- ifelse((df$race=='B' | df$race=='A' | df$race=='BA')&(df$ethnicity=='N'), 1, 0)
		#NW <- ifelse((df$race=='W')&(df$ethnicity=='N'), 1, 0)
		#NB <- ifelse((df$race=='B')&(df$ethnicity=='N'), 1, 0)
		#df$race[NB==1] <- 0
		#df$race[NW==1] <- 1
		#df$race[NB==0 & NW==0] <- 0
		#df$race <- as.numeric(df$race)
		df$race <- ifelse((df$race=='W')&(df$ethnicity=='N'), 1, 0) #ifelse(df$race=='W', 1, 0)
		df$race <- as.numeric(df$race)
		df <- df[,-7]# remove ethnicity column
		df <- df[,-5]# remove reason column
		cols_Y_or_N <- c("risk_new", "risk_mult", "risk_contact", "risk_msm",
			"risk_none", "symptoms", "sign_cf", "sign_c", "sign_pid", 
			"sign_u", "sign_no", "sign_none")

		for(col in cols_Y_or_N) df[,col] <- ifelse(df[,col]=='Y', 1, 0)

		df$ct <- ifelse(df$ct=='P',1,0)
		df$ct_pool <- ifelse(df$ct_pool=="P", 1, 0)

		#df$ct <- fitted$y_hat_med

		# indv swab: 1
		indv_swab <- df[is.na(df$pool_id)&df$specimen=='Swab',]
		indv_swab$seq_id <- replace_id(sort(unique(indv_swab$seq_id)), 1:length(unique(indv_swab$seq_id)),indv_swab$seq_id)
		indv_swab$pool_size <- 1
		indv_swab$assay <- 1
		indv_swab_ids <- matrix(-99, nrow(indv_swab), 4)
		indv_swab_ids[,1] <- indv_swab$seq_id
		indv_swab$pool_id <- indv_swab$seq_id
		Z_indv_swab <- cbind(indv_swab$ct, indv_swab$pool_size, indv_swab$assay, indv_swab_ids)
		X_indv_swab <- indv_swab[,5:18]
		S_indv_swab <- indv_swab$site_id
		Y_indv_swab <- cbind(indv_swab$ct, 1, indv_swab$pool_id, -99)

		# indv urine: 2
		indv_urine <- df[is.na(df$pool_id)&df$specimen=='Urine',]
		indv_urine$seq_id <- replace_id(sort(unique(indv_urine$seq_id)), 1:length(unique(indv_urine$seq_id)),indv_urine$seq_id)
		indv_urine$seq_id <- indv_urine$seq_id + max(indv_swab$seq_id)
		indv_urine$pool_size <- 1
		indv_urine$assay <- 2
		indv_urine_ids <- matrix(-99, nrow(indv_urine), 4)
		indv_urine_ids[,1] <- indv_urine$seq_id
		indv_urine$pool_id <- indv_urine$seq_id
		Z_indv_urine <- cbind(indv_urine$ct, indv_urine$pool_size, indv_urine$assay, indv_urine_ids)
		X_indv_urine <- indv_urine[,5:18]
		S_indv_urine <- indv_urine$site_id
		Y_indv_urine <- cbind(indv_urine$ct, 1, indv_urine$pool_id, -99)

		# pool swab: 3
		pool_swab <- df[!is.na(df$pool_id),]
		pool_swab$seq_id <- 1:length(unique(pool_swab$seq_id))
		pool_swab$seq_id <- pool_swab$seq_id + max(indv_urine$seq_id)

		for(pool_id in as.numeric(names(table(pool_swab$pool_id)))){
			pool_swab[pool_swab$pool_id==pool_id,"pool_size"] <- table(pool_swab$pool_id)[pool_id]
		}

		pool_swab$assay <- 3
		npools_swab <- length(unique(pool_swab$pool_id))
		Z_pool_swab <- cbind(rep(-99, npools_swab*4), 4, 3, matrix(-99, npools_swab*4, 4))

		idx <- 1
		for(pool_id in unique(pool_swab$pool_id)){
			tmp <- pool_swab[pool_swab$pool_id==pool_id,]
			pool_size <- max(tmp$pool_size)
			Z_pool_swab[idx,2] <- pool_size
			Z_pool_swab[idx,4:(3+pool_size)][1:pool_size] <- tmp$seq_id

			#Y_pool_swab[pool_swab$seq_id %in% tmp$seq_id, 3] <- idx

			if(max(tmp$ct_pool)==0){
				Z_pool_swab[idx,1] <- 0
			}else{
				
				Z_pool_swab[idx,1] <- 1

				# resolve from pool
				Z_pool_swab[(idx+1):(idx+pool_size),1] <- tmp$ct
				# pool size
				Z_pool_swab[(idx+1):(idx+pool_size),2] <- 1
				# assay indv swab 
				Z_pool_swab[(idx+1):(idx+pool_size),3] <- 1
				Z_pool_swab[(idx+1):(idx+pool_size),4] <- tmp$seq_id

				# # Y matrix
				# Y_pool_swab[pool_swab$seq_id %in% tmp$seq_id, 4] <- (idx+1):(idx+pool_size)
				# Y_pool_swab[pool_swab$seq_id %in% tmp$seq_id, 2] <- 2

				idx <- idx + pool_size

			}

			idx <- idx + 1

		}

		Z_pool_swab <- Z_pool_swab[Z_pool_swab[,1]>=0,]

		Y_pool_swab <- cbind(pool_swab$seq_id, 0, matrix(-99, nrow(pool_swab), 2))

		i <- 0
		for(seq_id in pool_swab$seq_id){
			print(round(i/length(pool_swab$seq_id), 4)*100)
			for(pool_id in 1:nrow(Z_pool_swab)){
				indvs <- Z_pool_swab[pool_id, -c(1:3)]
				if(seq_id %in% indvs[indvs>0]){
					cond <- Y_pool_swab[,1]==seq_id
					Y_pool_swab[cond, (3+Y_pool_swab[cond,2])] <- pool_id + max(indv_urine$pool_id)
					Y_pool_swab[cond, 2] <- Y_pool_swab[cond,2] + 1
				}
			}
			i <- i + 1
		}
		Y_pool_swab[,1] <- -99

		Y_pool_swab[(Y_pool_swab[,2]==1),1] <- 0
		Y_pool_swab[(Y_pool_swab[,2]==2),1] <- Z_pool_swab[Y_pool_swab[(Y_pool_swab[,2]==2),4]-max(indv_urine$seq_id)]

		X_pool_swab <- pool_swab[,5:18]
		S_pool_swab <- pool_swab$site_id

		# together
		Z <- rbind(Z_indv_swab, Z_indv_urine, Z_pool_swab)
		Y <- rbind(Y_indv_swab, Y_indv_urine, Y_pool_swab)
		S <- c(S_indv_swab, S_indv_urine, S_pool_swab)
		X <- rbind(X_indv_swab, X_indv_urine, X_pool_swab)

		N <- nrow(X)

		data <- list(Z=Z, Y=Y, S=S, X=X)
		saveRDS(data, 'data/formated.RData')

	}else{
		
		formated <- readRDS(glue("{folder}/formated.RData"))

		S <- formated$S
		X <- formated$X
		Z <- formated$Z
		Y <- formated$Y

		u <- do_normalization(X[,1], log=FALSE)
		
		race <- X[,2]
		#X[,2] <- ifelse(X[,2]==1, 1, 0)

		#X <- cbind(model.matrix(~race-1)[,-1], X[,-c(1:2)])

		u_seq <- seq(-3, 3, length.out=N_test)
		#X <- as.matrix(cbind(1, apply(X[,-1], 2, function(x) (x-mean(x))/sd(x))))
		X <- as.matrix(cbind(1, X[,-1]-0.5))
		X <- X[,c(1:5, 8:12, 6,7,13:14)]

		if(!is.null(vars)){
			X <- X[,vars]
		}

		if(twoway) X <- do_twoway(X)

		# if(!is.null(linear)){
		# 	X_lin <- X[, linear]
		# 	X <- as.matrix(X[, -linear])
		# }else{
		# 	X_lin <- NULL
		# }

		X_lin <- X
		
		N <- length(u)
		
		npools <- nrow(Z)
		pools <- Z[,2]
	
		# no varying intercept
		# X <- X[,-1]

		data <- list(testing="DT", X=as.matrix(X), u=u, X_lin=as.matrix(X_lin),
			N=N, npools=npools, pools=pools, N_test=N_test, 
			u_seq=u_seq, se=NULL, sp=NULL, S=S, nsites=length(unique(S)), Z=Z, Y=Y, race=race)

		return(data)
	}	
}






























