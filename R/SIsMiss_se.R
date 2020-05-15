#' Estimate standard error using asymptotic formula for unregularized regression.
#'
#' @param y Vector. The outcome variable
#' @param z Vector. The shadow variable (fully observed).
#' @param u Matrix. The covariate matrix \eqn{U}.
#' @param gamma_tilde_hat Initial estimates of scaled \eqn{\gamma}.
#' @param beta_hat Initial estimates of \eqn{\beta}.
#' @param se.method Charactor.
#' The method for estimating standard error of parameters.
#' Three options are available.
#' \describe{
#' \item{\code{"asymp"}} Obtain estimates of standard error via asymptotic
#' distribution. This option is applicable to \code{regularize=FALSE} only.
#' \item{\code{"perturb"}} Obtain estimates of standard error via perturbation method.
#' This option is applicable to both \code{regularize=FALSE} and \code{regularize=TRUE}.
#' \item{\code{NULL}} Standard error is not to be estimated.
#' }
#' @param M Number of resampling in perturbation.
#' @param seed_num Number of seed to control randomness of perturbation term generated.
#' @param regularize If the objective is to fit a regularized linear regression
#'  for variable selection, \code{regularize = TRUE}. By default,
#'  \code{regularize = FALSE}.
#' @param data_pair Matrix. The pairwise version of original data obtained from
#' function \code{pair_trans()$data_pair}.
#' @param newcov_pair Matrix. The pairwise transformed covariates obtained from
#' function \code{pair_trans()$newcov_pair}.
#' @param obs_id Vector. The indice of missingness for subjects obtained from
#' function \code{pair_trans}.
#' @param w_hat Matrix. The estimates of \eqn{W} matrix in the offset term
#' obtained by \code{pair_trans}.
#' @return A vector. The standard error of estimators in unregularized regression.

NFMiss_se <- function(y, z, u, gamma_tilde_hat, beta_hat,
                      se.method,
                      M = 500, seed_num, regularize,
                      data_pair, newcov_pair, obs_id, w_hat){
  regularize_TF <- regularize
  p <- ncol(u)

  z_u_model <- z_u_fit(z, u)
  eta_est <- z_u_model$eta_est
  z.s_est <- z_u_model$z.sigma_est

  g_t_h <- gamma_tilde_hat
  b_h <- beta_hat

  data.pair <- data_pair
  newcov.pair <- newcov_pair
  obs.if <- obs_id
  w.hat <- w_hat

  if(se.method == "asymp"){
    G_est <- G_matrix(z = z, u = u, z.sigma_est = z.s_est)
    AB_est <- AB_matrix(y, z, u,
                        gamma_tilde_hat = g_t_h, beta_hat = b_h,
                        z.sigma = z.s_est, w_hat = w.hat)
    A_est <- AB_est$A_est
    B_est <- AB_est$B_est

    se_out <- se_asymp(y, z, u,
                       gamma_tilde_hat = g_t_h, beta_hat = b_h,
                       eta_est, z.sigma = z.s_est,
                       A_est, B_est, G_est)
    return(se_out)
  }
  else if(se.method == "perturb"){
    res_pert <- perturbate(M, y, seed_num, regularize = regularize_TF,
                           data_pair = data.pair, newcov_pair = newcov.pair,
                           obs_id = obs.if)
    return(res_pert)
  }
}
