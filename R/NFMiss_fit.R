#' Fit the primary linear model with Nuissance-Free conditional
#' likelihood method.
#'
#' @import glmnet
#'
#' @param y Vector. The response variable subject to missing data.
#' @param z Vector. The shadow variable (fully observed).
#' @param u Matrix. The covariate matrix excluding
#' shadow variable \eqn{z} (fully observed).
#' @param cov.names Vector. The vector of names for all the covariates.
#' The first element is the name of variable \eqn{z},
#' followed by the name of variable \eqn{u}.
#' @param regularize Logical. If \code{regularize=TRUE}, a regularized linear
#' model is going to be fitted. Adaptive LASSO penalty is adopt. Default
#' is \code{regularize=FALSE}.
#' @param se.method Charactor.
#' The method for estimating standard error of parameters.
#' Three options are available.
#' \describe{
#' \item{\code{"asymp"}} Obtain estimates of standard error via asymptotic
#' distribution. This option is applicable to \code{regularize=FALSE} only.
#' \item{\code{"perturb"}} Obtain estimates of standard error via
#' perturbation method. This option is applicable to both
#' \code{regularize=FALSE} and \code{regularize=TRUE}.
#' \item{\code{NULL}} Standard error is not to be estimated.
#' }
#' @param CI.alpha \eqn{\alpha} value for confidence interval of each parameter.
#' \code{CI.alpha} must be a value between (0, 1).
#' If \code{CI.alpha=NULL}, confidence interval will not be returned.
#' @param M Number of resampling in perturbation.
#' @param seed_num Number of seed to control randomness of perturbation
#' term generated.
#'
#' @return A table. Summary of estimates and standard error presented in the table
#' with \eqn{p+1} rows and 2 columns.

NFMiss_fit <- function(y, z, u, regularize = FALSE, cov.names = NULL,
                       se.method = NULL, CI.alpha = NULL,
                       M = 500, seed_num = 123){
  regularize_TF <- regularize
  if(regularize_TF != TRUE & regularize_TF!= FALSE){
    stop("Argument regularize is invalid.")
  }

  ## Check conflict of SE estimation method and regularization option
  if(se.method == "asymp" & regularize == TRUE){
    stop("Estimation of standard error via asymptotic approach is
         applicable to unregularized regression only!")
  }

  ## Check conflict of SE estimation method and CI option
  if(is.null(se.method) == TRUE & is.null(CI.alpha) == FALSE){
    stop("Specify se.method to enable confidence interval estimation.")
  }



  N <- length(y)
  n <- sum(is.na(y) == FALSE)
  p <- ncol(u)
  r <- as.numeric(is.na(y))

  ## Default setting of variable names
  if(is.null(cov.names) == TRUE){
    prefix <- "u."
    suffix <- seq(1:p)
    cov_names <- c("z", paste(prefix, suffix, sep = ""))
  }else{
    cov_names <- cov.names
  }

  z_u_model <- z_u_fit(z, u)
  eta_est <- z_u_model$eta_est
  z.sigma_est <- z_u_model$z.sigma_est

  pair_trans_res <- pair_trans(y, z, u, eta_est, z.sigma_est)
  data.pair <- pair_trans_res$data_pair
  newcov.pair <- pair_trans_res$newcov_pair
  obs.if <- pair_trans_res$obs_id
  w.hat <- pair_trans_res$w_hat

  print("pair_trans() completed.")

  log.fit <- glm(newcov.pair[obs.if, 1] ~ 0+newcov.pair[obs.if, 2:(p+2)],
                 offset = newcov.pair[obs.if,p+3],
                 family = binomial(link = "logit"))
  theta.init <- as.numeric(log.fit$coefficients)
  gamma_tilde.init <- theta.init[1]
  beta.init <- theta.init[-1]/gamma_tilde.init

  ## Obtain estimates (unregularize / regularize)
  if(regularize_TF == FALSE){
    ## Output estimates
    res_est <- c(gamma_tilde.init, beta.init)
    print("Estimates obtained.")
    ## Output standard error
    if(se.method == "asymp"){
      res_se <- NFMiss_se(y, z, u,
                          gamma_tilde_hat = gamma_tilde.init,
                          beta_hat = beta.init,
                          se.method,
                          M = NULL, seed_num = NULL, regularize_TF,
                          data_pair = NULL, newcov_pair = NULL,
                          obs_id = NULL, w_hat = w.hat)
      res_lb <- res_est - qnorm(1-CI.alpha/2, 0, 1)*res_se
      res_ub <- res_est + qnorm(1-CI.alpha/2, 0, 1)*res_se
      print("Standard error obtained")
    }
    else if(se.method == "perturb"){
      res_pert <- NFMiss_se(y, z, u, gamma_tilde_hat = NULL, beta_hat = NULL,
                            se.method,
                            M, seed_num, regularize_TF,
                            data.pair, newcov.pair,
                            obs.if, w_hat = NULL)
      res_se <- apply(res_pert, 1, function(x){sd(x, na.rm = TRUE)})
      res_lb <- apply(res_pert, 1, function(x){quantile(x, CI.alpha/2, na.rm = TRUE)})
      res_ub <- apply(res_pert, 1, function(x){quantile(x, 1-CI.alpha/2, na.rm = TRUE)})
      print("Standard error obtained")
    }
    else if(is.null(se.method) == TRUE){
      res_se <- rep(".", p+1)
    }
    else{stop("Invalid option for se.method")}
  }else{
    p.fac <- c(0,abs(theta.init[-1])^(-1))

    ## Adaptive LASSO
    print("Start to fit adaptive LASSO.")
    fit.adalasso <- glmnet(x = newcov.pair[obs.if, 2:(p+2)],
                           y = newcov.pair[obs.if, 1],
                           intercept = FALSE,
                           offset = newcov.pair[obs.if, p+3],
                           family = 'binomial',
                           penalty.factor = p.fac
                           )
    lambda.seq <- fit.adalasso$lambda
    df.seq <- fit.adalasso$df
    lambda.length <- length(lambda.seq)

    ## Compute BIC
    BIC <- rep(NA, lambda.length)
    print("Start to evaluate BIC.")
    for(b in 1:lambda.length){
      #print(b)
      ## Remove the intercept term
      fit.theta <- as.numeric(coef(fit.adalasso, s = lambda.seq[b]))[-1]
      gamma.tilde <- fit.theta[1]
      beta.tilde <- fit.theta[-1]

      ##BIC
      loglik.vector=apply(data.pair[obs.if,],
                          1,
                          function(x){
                            t <- -gamma.tilde*x[2]*x[1]+x[2]*(x[3:(p+2)]%*%beta.tilde)+log(x[p+3])
                            if(exp(t) == Inf){
                              return(t)
                            }
                            else{
                              return(log(1+exp(t)))
                            }
                          }
                          )
      loglik.function <- 2*sum(loglik.vector)/(n*(n-1))

      BIC[b] <- 2*(loglik.function)+(df.seq[b])*log(n)/n
    }

    print("Evaluate BIC completed.")

    lambda.min <- lambda.seq[which.min(BIC)]

    theta.est.BIC <- as.numeric(coef(fit.adalasso, s = lambda.min))
    gamma.tilde_est <- theta.est.BIC[2]
    beta.est <- theta.est.BIC[-c(1, 2)]/gamma.tilde_est
    res_est <- c(gamma.tilde_est, beta.est)
    print("Estimates obtained")

    if(is.null(se.method) == TRUE){
      res_se <- rep(".", p+1)
    }else{
      print("Start to evaluate standard error.")
      res_pert <- NFMiss_se(y, z, u, gamma_tilde_hat = NULL, beta_hat = NULL,
                            se.method,
                            M, seed_num, regularize_TF,
                            data.pair, newcov.pair,
                            obs.if, w_hat = NULL)
      res_se <- apply(res_pert, 1, function(x){sd(x, na.rm = TRUE)})
      res_lb <- apply(res_pert, 1, function(x){quantile(x, CI.alpha/2, na.rm = TRUE)})
      res_ub <- apply(res_pert, 1, function(x){quantile(x, 1-CI.alpha/2, na.rm = TRUE)})
      print("Standard error obtained")
    }
  }


  if(is.null(CI.alpha) == FALSE){
    if(CI.alpha>0 & CI.alpha<1){
      res_out <- as.data.frame(matrix(data = NA, ncol = 4, nrow = p+1))
      colnames(res_out) <- c("Estimates", "Standard Error",
                             "Lower Bound", "Upper Bound")
      res_out[, 1] <- res_est
      res_out[, 2] <- res_se
      res_out[, 3] <- res_lb
      res_out[, 4] <- res_ub
    }
    else{
      stop("CI.alpha is invalid.")
    }
  }else if(is.null(CI.alpha) == TRUE){
    res_out <- as.data.frame(matrix(data = NA, ncol = 2, nrow = p+1))
    colnames(res_out) <- c("Estimates", "Standard Error")
    res_out[, 1] <- res_est
    res_out[, 2] <- res_se
  }
  rownames(res_out) <- cov_names

  return(res_out)
}
