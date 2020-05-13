#' Generate covariates to fit the logistic regression,
#' which is the form of reparameterized objective function.
#'
#' @description \code{pair_trans} is used to generate transformed covariates in pairwise form.
#'
#' @param y Vector. The response variable subject to missing data.
#' @param z Vector. The shadow variable (fully observed).
#' @param u Matrix. The covariate matrix excluding shadow variable \eqn{z} (fully observed).
#' @param eta_est The estimates for coefficients in model \eqn{Z} given\eqn{U}.
#' @param z.sigma_est The estimates for \eqn{\sigma} in model \eqn{Z} given\eqn{U}.
#'
#' @return A list:
#' \describe{
#' \item{data_pair}{The matrix of pairwised \eqn{y}, \eqn{z}, \eqn{u}.}
#' \item{newcov_pair}{The matrix of pairwised transformed variables to be fit in the objective function.}
#' \item{obs_id}{The indice of missingness for subjects in \code{data}.}
#' \item{w_hat}{The estimates of \eqn{W} matrix in the offset term.}
#' }
#'
#'

pair_trans <- function(y, z, u, eta_est, z.sigma_est){
  if(sum(is.na(z)) > 0 | sum(!complete.cases(u)) > 0){
    stop("Remove missing value in covariates.")
  }else{
    N <- length(y)
    p <- ncol(u)
    r <- as.numeric(!is.na(y))
    n <- sum(r)

    # generate a matrix to store transformed covariates
    data.pair <- matrix(data = NA, nrow = N*(N-1)/2, ncol = p+3)
    newcov.pair <- matrix(data = NA, nrow = N*(N-1)/2, ncol = p+3)
    w.hat <- matrix(data = NA, nrow = N, ncol = N)
    obs.if <- c()
    k <- 0
    for(j in 2:N){
      for(i in 1:(j-1)){
        k <- k+1
        if(r[i] == 0 | r[j] == 0){
          w.hat[i,j] <- NA
          obs.if[k] <- FALSE
        }
        else{
          obs.if[k] <- TRUE
          y.ij <- y[i] - y[j]
          u.ij <- as.numeric(u[i,]-u[j,])
          u.ij_eta <- as.numeric(t(u.ij)%*%eta_est[-1])
          z.ij <- z[i] - z[j]
          # Offset w
          w.hat[i,j] <- exp(-z.ij*u.ij_eta/z.sigma_est^2)
          # Paired data
          # y.ij, z.ij, u.ij, w.ij
          data.pair[k,1] <- y.ij
          data.pair[k,2] <- z.ij
          data.pair[k,3:(p+2)] <- u.ij
          data.pair[k,p+3] <- w.hat[i,j]

          ##New covariates
          #response
          newcov.pair[k,1] <- as.numeric(z.ij>0)
          #cov for gamma.tilde
          newcov.pair[k,2] <- abs(z.ij)*y.ij
          #cov for beta.tilde
          newcov.pair[k,3:(p+2)] <- -abs(z.ij)*u.ij
          #offset term
          newcov.pair[k,p+3] <- -sign(z.ij)*log(w.hat[i,j])
        }
      }
    }
    return(list(data_pair = data.pair,
                newcov_pair = newcov.pair,
                obs_id = obs.if,
                w_hat = w.hat))
  }
}




