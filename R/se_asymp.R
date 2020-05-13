#' Generate the covariance matrix of estimators in unregularized regression.
#' @param y Vector. The outcome variable
#' @param z Vector. The shadow variable (fully observed).
#' @param u Matrix. The covariate matrix \eqn{u}.
#' @param gamma_tilde_hat Initial estimates of scaled \eqn{\gamma}.
#' @param beta_hat Initial estimates of \eqn{\beta}.
#' @param z.sigma The estimates for \eqn{\sigma} in model \eqn{Z} given\eqn{U}.
#' @param A_est Matrix. The estimates of \eqn{A} matrix obtained by \code{AB_matrix}.
#' @param B_est Matrix. The estimates of \eqn{B} matrix obtained by \code{AB_matrix}.
#' @param G_est Matrix. The estimates of \eqn{G} matrix obtained by \code{G_matrix}
#' @return A vector. The standard error of estimators in unregularized regression.
#'
se_asymp <- function(y, z, u, gamma_tilde_hat, beta_hat,
                     eta_est, z.sigma, A_est, B_est, G_est){
  N <- length(y)
  p <- ncol(u)
  r <- as.numeric(!is.na(y))

  gamma_tilde.init <- gamma_tilde_hat
  beta.init <- beta_hat * gamma_tilde_hat
  eta.est <- eta_est

  A <- A_est
  B <- B_est
  G <- G_est

  # Design matrix C
  C_matrix <- diag(c(1, rep(1/gamma_tilde_hat, p)))

  sigma.2.est <- matrix(data = 0, nrow = p+1, ncol = p+1)
  for(i in 1:N){
    sigma.2.vec.core <- rep(0, p+1)
    for(j in setdiff(1:N, i)){
      z.ij <- z[i] - z[j]
      u.ij <- as.numeric(u[i, ] - u[j, ])
      u.i.design <- c(1, as.numeric(u[i, ]))
      u.j.design <- c(1, as.numeric(u[j, ]))
      u.ij.design <- c(1, u.ij)
      y.ij <- y[i] - y[j]
      zeta.ij <- gamma_tilde.init*z.ij*(as.numeric(u.ij%*%beta.init) - y.ij)
      w.ij <- exp(-z.ij*as.numeric(u.ij%*%eta.est[-1])/z.sigma^2)
      exp.ij <- exp(zeta.ij + log(w.ij))

      M.i <- (z[i] - as.numeric(t(eta.est)%*%u.i.design))*u.i.design/z.sigma^2
      M.j <- (z[j] - as.numeric(t(eta.est)%*%u.j.design))*u.j.design/z.sigma^2
      M.ij <- (M.i + M.j)/2
      if(r[i] == 1 & r[j] == 1){
        N.ij.vector <- c(z.ij*(as.numeric(u.ij%*%beta.init) - y.ij), gamma_tilde.init*z.ij*u.ij)
        N.ij <- -(exp.ij/(1 + exp.ij))*N.ij.vector
      }else{
        N.ij <- rep(0, p+1)
      }

      sigma.2.vec.core <- sigma.2.vec.core + (B%*%solve(G)%*%M.ij - N.ij)/(N - 1)
    }
    sigma.2.est <- sigma.2.est + sigma.2.vec.core%*%t(sigma.2.vec.core)*(4/(N - 1))
  }

  sigma.est <- solve(A)%*%sigma.2.est%*%solve(A)/N
  se_out <- sqrt(diag(C_matrix %*% sigma.est %*% C_matrix))

  #names(se_out) <- name_out

  return(se_out)
}
