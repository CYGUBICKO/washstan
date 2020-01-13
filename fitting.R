library(brms)
library(dplyr)
library(mvtnorm)
library(rstan)

set.seed(7777)

if (file.exists("generateAR1.rda")){
   load("generateAR1.rda")
} else{
   source("generateAR1.R")
}

sim_data <- sim_dflist[[1]]
print(sim_data, n = 100)

## Long format for non-bayesian methods
services <- c("y1", "y2", "y3")
long_df <- (sim_data
	%>% select(c("hhid", "years", "wealthindex", services))
	%>% gather(service, status, services)
)


# Fitting on continous outcome
brms1 <- brm(y1 ~ wealthindex + (1|hhid) + (1|years)
	, data = sim_data
	, chains = 2
	, cores = 4
	, autocor=cor_ar(formula = ~1, p = 1, cov = FALSE)
)

brms2 <- brm(mvbind(y1, y2, y3) ~ wealthindex + (1|hhid) + (1|years)
	, data = sim_data
	, chains = 2
	, cores = 4
	, autocor=cor_ar(formula = ~1, p = 1, cov = FALSE)
)

## Brms binary response

priors <- c(prior(normal(0, 1), class = b, resp = y1bin)
   , prior(normal(0, 1), class = b, resp = y2bin)
   , prior(normal(0, 1), class = b, resp = y3bin)
   , prior(normal(0, 1), class = b, coef = intercept, resp = y1bin)
   , prior(normal(0, 1), class = b, coef = intercept, resp = y2bin)
   , prior(normal(0, 1), class = b, coef = intercept, resp = y3bin)
   , prior(normal(0, 1), class = b, coef = wealthindex, resp = y1bin)
   , prior(normal(0, 1), class = b, coef = wealthindex, resp = y2bin)
   , prior(normal(0, 1), class = b, coef = wealthindex, resp = y3bin)
   , prior(cauchy(0, 5), class = sd, resp = y1bin)
   , prior(cauchy(0, 5), class = sd, resp = y2bin)
   , prior(cauchy(0, 5), class = sd, resp = y3bin)
   , prior(cauchy(0, 5), class = sd, group = hhid, resp = y1bin)
   , prior(cauchy(0, 5), class = sd, group = hhid, resp = y2bin)
   , prior(cauchy(0, 5), class = sd, group = hhid, resp = y3bin)
   , prior(cauchy(0, 5), class = sd, coef = Intercept, group = hhid, resp = y1bin)
   , prior(cauchy(0, 5), class = sd, coef = Intercept, group = hhid, resp = y2bin)
   , prior(cauchy(0, 5), class = sd, coef = Intercept, group = hhid, resp = y3bin)
   , prior(cauchy(0, 5), class = sd, group = years, resp = y1bin)
   , prior(cauchy(0, 5), class = sd, group = years, resp = y2bin)
   , prior(cauchy(0, 5), class = sd, group = years, resp = y3bin)
   , prior(cauchy(0, 5), class = sd, coef = Intercept, group = years, resp = y1bin)
   , prior(cauchy(0, 5), class = sd, coef = Intercept, group = years, resp = y2bin)
   , prior(cauchy(0, 5), class = sd, coef = Intercept, group = years, resp = y3bin)
)

brmsbin <- brm(mvbind(y1bin, y2bin, y3bin) ~ 0 + intercept + wealthindex + (1|hhid) + (1|years)
	, data = sim_data
	, chains = 2
	, cores = 4
	, family = list(bernoulli(link = "logit")
		, bernoulli(link = "logit")
		, bernoulli(link = "logit")
	)
	, prior = priors
	, autocor=cor_ar(formula = ~hhid|years, p = 1, cov = TRUE)
)

### Stan model
stan_df <- (sim_data
#	%>% filter(years==2007)
	%>% select(wealthindex, y1bin, y2bin, y3bin)
)

standata <- list(
	N=nrow(stan_df)
	, K=2
	, D=ncol(stan_df)-1
	, y=cbind(stan_df$y1bin, stan_df$y2bin, stan_df$y3bin)
	, x=cbind(1, stan_df$wealthindex)
)

rt <- stanc("mvt_probit.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)
stan1 <- sampling(sm
	, data = standata
	, chains = 2
	, cores = 4
	, iter = 2000
	, thin = 1
	, seed = 101
)
print(stan1, pars=c("beta"))

save(file = "fitting.rda"
	, brms1
	, brms2
	, brmsbin
	, stan1
	, betas_df
	, covmat_df
)
