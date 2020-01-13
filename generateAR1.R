#### ---- Project: APHRC Wash Data ----
#### ---- Task: Simulation ----
#### ---- Sub-task: Simulate multivariate gausian responses ----
#### ---- By: Steve and Jonathan ----
#### ---- Date: 2019 May 28 (Tue) ----

library(mvtnorm)
library(data.table)
library(dplyr)
options(dplyr.width = Inf)

library(tidyr)
library(tibble)
library(nlme)

set.seed(7777)

# Simulation parameters
nsims <- 1		# Number of simulations to run
nHH <- 100		# Number of HH (primary units) per year

nyrs <- 30	# Number of years to simulate
yrs <- 2000 + c(1:nyrs) # Years to simulate
N <- nyrs * nHH

# Generate dataset template
temp_df <- (data.frame(hhid = rep(c(1:nHH), each = nyrs)
		, years = rep(yrs, nHH)
		, wealthindex = rnorm(n = N)
	)
	%>% group_by(hhid)
	%>% ungroup()
)

# Beta values
y1_beta0 <- 0.3
y1_beta1 <- 0.4 
y2_beta0 <- 0.3
y2_beta1 <- 0.8
y3_beta0 <- 0.4
y3_beta1 <- 0.5

# Correlation matrix
cor_y1y2 <- 0.20
cor_y1y3 <- 0.30
cor_y2y3 <- 0.50
corMat <- matrix(
	c(1, cor_y1y2, cor_y1y3
		, cor_y1y2, 1, cor_y2y3
		, cor_y1y3, cor_y2y3, 1
	), 3, 3
)

# Sd
y1_sd <- 0.5
y2_sd <- 0.3
y3_sd <- 0.7
sdVec <- c(y1_sd, y2_sd, y3_sd)
varMat <- sdVec %*% t(sdVec)
varMat
corMat
# varcov matrix
covMat <- varMat * corMat
covMat

# Simulate AR1 process: 
## https://newonlinecourses.science.psu.edu/stat502/lesson/10/10.3
## https://stats.stackexchange.com/questions/295102/how-to-write-variance-covariance-matrix-of-ar1-process-in-r

phi <- 0.8

simAr1 <- function(phi
		, sdar1
		, ngrp
		, npergrp){
	cmat <- sdar1*phi^abs(outer(0:(npergrp-1), 0:(npergrp-1), "-"))
	errs <- MASS::mvrnorm(ngrp, mu = rep(0, npergrp), Sigma = cmat)
	errs <- c(t(errs)) ## unpack errors by row
	return(errs)
}

# Generate dataset
sim_dflist <- list() # Simulated datasets
betas0_dflist <- list() # Simulated random effects estimates for Year 
hhRE_dflist <- list() # Simulated random effects datasets estimates for HH
ar1Errors_dflist <- list() # Simulated yearly ar1

for (i in 1:nsims){
	
	# AR1 process errors
	y1_ar1 <- simAr1(phi = phi, sdar1 = 1, ngrp = nHH, npergrp = nyrs)
	y2_ar1 <- simAr1(phi = phi, sdar1 = 1, ngrp = nHH, npergrp = nyrs)
	y3_ar1 <- simAr1(phi = phi, sdar1 = 1, ngrp = nHH, npergrp = nyrs)
	
	# Save the ar1 df
	ar1Errors <- (data.frame(X1 = y1_ar1, X2 = y2_ar1, X3 = y3_ar1)
		%>% setnames(c("X1", "X2", "X3"), c("y1", "y2", "y3"))
	)
	ar1Errors_dflist[[i]] <- ar1Errors
	
	# Simulate B0 for each year and then merge to HH data. Different HH has same B0 for same year
	betas0 <- (MASS::mvrnorm(nyrs
			, mu = c(y1_beta0, y2_beta0, y3_beta0)
			, Sigma = covMat
			, empirical = TRUE
		)
   	%>% data.frame()
		%>% mutate(years = yrs)
		%>% right_join(temp_df)
		%>% select(c("years", "X1", "X2", "X3"))
	)

	# Save the B0 RE estimates
	betas0_dflist[[i]] <- (betas0
		%>% distinct()
		%>% setnames(c("X1", "X2", "X3"), c("y1", "y2", "y3"))
	)

	# Remove year from betas0
	betas0 <- (betas0
		%>% select(-years)
	)

	# Simulate HH-level random effects (residual error)
	hhRE <- MASS::mvrnorm(nHH
		, mu = c(0, 0, 0)
		, Sigma = covMat
		, empirical = TRUE
	)
	hhRE <- hhRE[temp_df$hhid, ]
	
	# Save HH RE simulation estimates
	hhRE_dflist[[i]] <- (hhRE
		%>% data.frame()
		%>% mutate(hhid = pull(temp_df, hhid))
		%>% distinct()
	)
	
	dat <- (temp_df
		%>% mutate(y1 = betas0[,1] + y1_beta1*wealthindex + hhRE[,1] + ar1Errors[,1]
			, y2 = betas0[,2] + y2_beta1*wealthindex + hhRE[,2] + ar1Errors[,2]
			, y3 = betas0[,3] + y3_beta1*wealthindex + hhRE[,3] + ar1Errors[,3]
			, y1bin = rbinom(N, 1, plogis(y1))
			, y2bin = rbinom(N, 1, plogis(y2))
			, y3bin = rbinom(N, 1, plogis(y3))
		)
	)
	sim_dflist[[i]] <- dat
}

print((sim_dflist[[1]]))

# Extract beta values assigned in the simulation
betas <- sapply(grep("y[1-3]_beta|cor_y[1-3]|y[1-3]_sd$", ls(), value = TRUE), get)
#betas <- betas[!names(betas) %in% grep("_sd", names(betas), value = TRUE)]
betas_df <- (data.frame(betas) 
	%>% rownames_to_column("coef")
	%>% mutate(n = extract_numeric(coef)
		, coef = ifelse(grepl("_0$", coef)
			, paste0("y", n)
				, ifelse(grepl("_beta", coef) 
					, paste0("wealthindex:y", n)
				, coef
			)
		)
	)
)


# Create a dataframe of covariance matrix
covmat_df <- (
	data.frame(coef = c("y1y1", "y2y2", "y3y3", "y2y1", "y3y1", "y3y2")
		, values = c(diag(covMat), covMat[lower.tri(covMat)])
	)
	%>% mutate(n = extract_numeric(coef)
		, coef_clean = paste0("Sigma[years:y", substr(n, 1, 1), ",y", substr(n, 2, 2), "]")
	) 
)

# Maybe there is a better way but I am lazy!!!
ar1_objs <- list(phi = phi
	, sdar1 = c(y1_sd, y2_sd, y3_sd)
)

save(file = "generateAR1.rda"
	, sim_dflist
	, betas_df
	, covmat_df
	, betas0_dflist
	, hhRE_dflist
	, betas
	, corMat
	, covMat
	, ar1_objs
	, ar1Errors_dflist
)

