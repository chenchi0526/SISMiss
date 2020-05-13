rm(list=ls())
#.libPaths(c(.libPaths(), "~/R/lib"))

#ccr.task.id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

source("z_u_fit.R")
source("pair_trans.R")
source("NFMiss_se.R")
source("G_matrix.R")
source("AB_matrix.R")
source("se_asymp.R")
source("perturbate.R")
source("NFMiss_fit.R")

library(devtools)
library(glmnet)


icu_data <- read.csv("data3.csv")

## select columns
# demo: 8, 9, 10, 11
# chart: 39, 42, 45, 48, 51, 54
# lab: 15, 18, 21, 24, 27, 30, 32, 36, 56
# metrics: 58, 59

ind_demo <- c(8, 9, 10)
ind_chart <- c(39, 42, 45, 48, 51, 54)
ind_lab <- c(15, 18, 20, 24, 27, 30, 32, 36, 56)
ind_metrics <- c(58, 59)

# only include "ind" columns
icu_v1 <- icu_data[, c(ind_demo, ind_chart, ind_lab, ind_metrics)]
# remove NAs
icu_v2 <- icu_v1[complete.cases(icu_v1[, -which(colnames(icu_v1) == "albumin_mean")]), ]
# age 20-40 and insurance in (medicare & medicaid)
if_age <- icu_v2$age>=25 & icu_v2$age<=40
if_marital <- icu_v2$marital_status=="MARRIED"
#if_insurance <- icu_v2$insurance %in% c("Medicare", "Medicaid", "Government")
#icu_v3 <- icu_v2[if_age & if_marital & if_insurance, ]
icu_v3 <- icu_v2[if_age & if_marital, -which(colnames(icu_v2) == "marital_status")]

# scale "score" variables
icu_v3$albumin_mean <- scale(icu_v3$albumin_mean)
icu_v3$sapsii <- scale(icu_v3$sapsii)
icu_v3$sofa <- scale(icu_v3$sofa)

# factorize
icu_v3$gender <- as.numeric(icu_v3$gender == "F")

z_name <- "calcium_mean"

y <- icu_v3$albumin_mean

z <- icu_v3[ , which(colnames(icu_v3) == z_name)]
# remove column of y and z in u
u <- as.matrix(icu_v3[ , -which(colnames(icu_v3) %in% c("albumin_mean", z_name))])

## End of import data
######################################################################################

N <- length(y)
r <- as.numeric(!is.na(y))
n <- sum(is.na(y) == FALSE)
p <- ncol(u)


regularize = FALSE
cov.names = c("albumin_mean", colnames(u))
se.method = "asymp"
M = 500
seed_num = 123

if(is.null(cov.names) == TRUE){
  prefix <- "u."
  suffix <- seq(1:p)
  cov_names <- c("z", paste(prefix, suffix, sep = ""))
}else{
  cov_names <- cov.names
}

res_out <- as.data.frame(matrix(data = NA, ncol = 2, nrow = p+1))
rownames(res_out) <- cov_names
colnames(res_out) <- c("Estimates", "Standard Error")

z_u_model <- z_u_fit(z, u)
eta_est <- z_u_model$eta_est
z.sigma_est <- z_u_model$z.sigma_est


pair_trans_res <- pair_trans(y, z, u, eta_est, z.sigma_est)
data.pair <- pair_trans_res$data_pair
newcov.pair <- pair_trans_res$newcov_pair
obs.if <- pair_trans_res$obs_id
w.hat <- pair_trans_res$w_hat

log.fit <- glm(newcov.pair[obs.if, 1] ~ 0+newcov.pair[obs.if, 2:(p+2)],
               offset = newcov.pair[obs.if,p+3],
               family = binomial(link = "logit"))
theta.init_est <- as.numeric(log.fit$coefficients)
gamma_tilde.init_est <- theta.init_est[1]
beta.init_est <- theta.init_est[-1]/gamma_tilde.init_est

## Result of unregularized (se.method=asymp)
res_se_1 <- NFMiss_se(y, z, u,
                      gamma_tilde_hat = gamma_tilde.init_est,
                      beta_hat = beta.init_est,
                      se.method = "asymp",
                      M = NULL, seed_num = NULL, regularize = NULL,
                      data_pair = NULL, newcov_pair = NULL,
                      obs_id = NULL, w_hat = w.hat)

## Result of unregularized (se.method=perturb)
res_se_2 <- NFMiss_se(y, z, u,
                      gamma_tilde_hat = gamma_tilde.init_est,
                      beta_hat = beta.init_est,
                      se.method = "perturb",
                      M = 500, seed_num = 1234, regularize = FALSE,
                      data_pair = data.pair, newcov_pair = newcov.pair,
                      obs_id = obs.if, w_hat = NULL)
res_se_2 <- apply(res_se_2, )

## Result of regularized (se.method=perturb)

res_se_3 <- NFMiss_fit(y, z, u,
                       regularize = TRUE, cov.names = NULL,
                       se.method = "perturb", M = 500, seed_num = 1234)


########################################################################
## One step
## Result of unregularized (se.method=asymp)
res_out_1 <- NFMiss_fit(y, z, u, regularize = FALSE, cov.names = NULL,
                        se.method = "asymp", CI.alpha = 0.05,
                        M = NULL, seed_num = NULL)

## Result of unregularized (se.method=perturb)
res_out_2 <- NFMiss_fit(y, z, u, regularize = FALSE, cov.names = NULL,
                        se.method = "perturb", CI.alpha = 0.05,
                        M = 500, seed_num = 1234)

## Result of regularized (se.method=perturb)
res_out_3 <- NFMiss_fit(y, z, u, regularize = TRUE, cov.names = NULL,
                        se.method = "perturb", CI.alpha = 0.05,
                        M = 10, seed_num = 1234)
