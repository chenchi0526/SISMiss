#' Obtain the estimates for standard error via perturbation.
#'
#' @import glmnet
#'
#' @param M Number of resampling in perturbation.
#' @param y Vector. The outcome variable.
#' @param seed_num Number of seed to control randomness of perturbation term generated.
#' @param regularize If the objective is to fit a regularized linear regression
#'  for variable selection, \code{regularize = TRUE}. By default, \code{regularize = FALSE}.
#' @param data_pair Matrix. The pairwise version of original data obtained from
#' function \code{pair_trans()$data_pair}.
#' @param newcov_pair Matrix. The pairwise transformed covariates obtained from
#' function \code{pair_trans()$newcov_pair}.
#' @param obs_id Vector. The indice of missingness for subjects obtained from
#' function \code{pair_trans()$obs_id}.
#' @return A matrix. \code{M} sets of estimates based on perturbated objective function.


perturbate <- function(M, y, seed_num, regularize = FALSE,
                       data_pair, newcov_pair, obs_id){
  set.seed(seed_num)
  N <- length(y)
  r <- as.numeric(!is.na(y))
  n <- sum(is.na(y) == FALSE)
  p <- ncol(data_pair)-3

  data.pair <- data_pair
  obs.if <- obs_id
  newcov.pair <- newcov_pair

  v.exp <- matrix(data = NA, nrow = N, ncol = M)
  for(m in 1:M){
    v.exp[, m]=rexp(N, 1)
  }

  ## matrix: output of estimates from perturbation
  ## Results from one perturbation is stored in a column
  res_perturb <- matrix(data = NA, ncol = M, nrow = p+1)

  for(m in 1:M){
    print(paste0("m=", m))

    comb_v.exp <- rep(NA, N*(N-1)/2)

    k=0
    for(j in 2:N){
      for(i in 1:(j-1)){
        k <- k+1
        comb_v.exp[k] <- v.exp[i, m] * v.exp[j, m]
      }
    }
    log.fit <- glm(newcov.pair[obs.if, 1] ~ 0+newcov.pair[obs.if, 2:(p+2)],
                   weights = comb_v.exp[obs.if],
                   offset = newcov.pair[obs.if, p+3],
                   family = binomial(link = "logit"))
    theta.init <- as.numeric(log.fit$coefficients)
    gamma_tilde_init_est <- theta.init[1]
    beta_init_est <- theta.init[-1]/gamma_tilde_init_est

    ## Unregularized
    if(regularize == FALSE){
      res_perturb[, m] <- c(gamma_tilde_init_est, beta_init_est)
    }else if(regularize == TRUE){
      ## Regularized
      #theta.init <- as.numeric(log.fit$coefficients)
      p.fac <- c(0, abs(theta.init[-1])^(-1))
      ##ALASSO for BIC
      fit.adalasso <- glmnet(x = newcov.pair[obs.if, 2:(p+2)],
                             y = newcov.pair[obs.if, 1],
                             intercept = FALSE,
                             weights = comb_v.exp[obs.if],
                             offset = newcov.pair[obs.if, p+3],
                             family = 'binomial',
                             penalty.factor = p.fac)
      lambda.seq <- fit.adalasso$lambda
      df.seq <- fit.adalasso$df
      lambda.length <- length(lambda.seq)

      ##BIC(n)
      BIC <- c()

      data.pair_with_v <- cbind(data.pair, comb_v.exp)
      for(m_tun in 1:lambda.length){
        #print(m_tun)

        fit.theta <- as.numeric(coef(fit.adalasso, s = lambda.seq[m_tun]))[-1]
        gamma.tilde <- fit.theta[1]
        beta.tilde <- fit.theta[-1]

        ##BIC
        loglik.vector <- apply(data.pair_with_v[obs.if, ],
                               1,
                               function(x){
                                 t <- -gamma.tilde*x[2]*x[1]+x[2]*(x[3:(p+2)]%*%beta.tilde)+log(x[p+3])
                                 if(exp(t) == Inf){
                                   return(t*x[p+4])
                                 }
                                 else{
                                   return(log(1+exp(t))*x[p+4])
                                 }
                               }
        )

        loglik.function <- 2*sum(loglik.vector)/(n*(n-1))

        BIC[m_tun] <- 2*(loglik.function)+(df.seq[m_tun])*log(n)/n
      }
      lambda.min <- lambda.seq[which.min(BIC)]
      theta.est.BIC <- as.numeric(coef(fit.adalasso, s = lambda.min))
      gamma_tilde.est <- theta.est.BIC[2]
      beta.est <- theta.est.BIC[-c(1, 2)]/theta.est.BIC[2]

      res_perturb[, m] <- c(gamma_tilde.est, beta.est)
    }
    else{stop("Specify correct regularization option.")}
  }

  return(res_perturb)
}
