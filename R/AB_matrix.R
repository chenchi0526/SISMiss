#' Generate the matrix A and B for variance estimation in low dimension scenario.
#' @param y Vector. The outcome variable
#' @param z Vector. The shadow variable (fully observed).
#' @param u Matrix. The covariate matrix \eqn{u}.
#' @param gamma_tilde_hat Initial estimates of scaled \eqn{\gamma}.
#' @param beta_hat Initial estimates of \eqn{\beta}.
#' @param z.sigma The estimates for \eqn{\sigma} in model \eqn{Z} given\eqn{U}.
#' @param w_hat The estimates of \eqn{W} matrix in the offset term
#' obtained by \code{pair_trans}.
#' @return A list.
#' \describe{
#' \item{A_est}{The estimates for matrix \eqn{A}.}
#' \item{B_est}{The estimates for matrix \eqn{B}.}
#' }

AB_matrix <- function(y, z, u, gamma_tilde_hat, beta_hat, z.sigma, w_hat){
  N <- length(y)
  p <- ncol(u)
  r <- as.numeric(!is.na(y))

  w <- w_hat

  #Need to scale gamma here
  gamma_tilde.init <- gamma_tilde_hat
  beta.init <- beta_hat * gamma_tilde_hat
  A <- matrix(data = 0, ncol = p+1, nrow = p+1)
  B <- matrix(data = 0, ncol = p+1, nrow = p+1)
  for(j in 2:N){
    for(i in 1:(j-1)){
      if(r[i] == 1 & r[j] == 1){
        z.ij <- z[i]-z[j]
        u.ij <- as.numeric(u[i,] - u[j,])
        y.ij <- y[i] - y[j]
        zeta.ij <- as.numeric(gamma_tilde.init*z.ij*(u.ij%*%beta.init - y.ij))
        exp.ij <- as.numeric(exp(zeta.ij + log(w[i,j])))

        ##Components in A
        A.part.1.matrix <- matrix(data = NA, nrow = p+1, ncol = p+1)
        A.part.1.matrix[1, 1] <- 0
        A.part.1.matrix[2:(p+1), 2:(p+1)] <- 0
        A.part.1.matrix.top <- z.ij*u.ij
        A.part.1.matrix[1, 2:(p+1)] <- A.part.1.matrix.top
        A.part.1.matrix[2:(p+1), 1] <- A.part.1.matrix.top
        A.part.1 <- (exp.ij/(1 + exp.ij))*A.part.1.matrix
        A.part.2.matrix <- matrix(data = NA, nrow = p+1, ncol = p+1)
        A.part.2.matrix[1, 1] <- (z.ij*(u.ij%*%beta.init - y.ij))^2
        A.part.2.matrix.top <- as.numeric(gamma_tilde.init*(z.ij^2)*(u.ij%*%beta.init - y.ij))*u.ij
        A.part.2.matrix[1,2:(p+1)] <- A.part.2.matrix.top
        A.part.2.matrix[2:(p+1),1] <- A.part.2.matrix.top
        A.part.2.matrix.core <- ((gamma_tilde.init*z.ij)^2)*(u.ij%*%t(u.ij))
        A.part.2.matrix[2:(p+1), 2:(p+1)] <- A.part.2.matrix.core
        A.part.2 <- (exp.ij/(1 + exp.ij)^2)*A.part.2.matrix
        A.delta <- -2*(A.part.1 + A.part.2)/(N*(N - 1))

        ##Components in B
        theta.partial <- c(z.ij*(as.numeric(t(beta.init)%*%u.ij) - y.ij), gamma_tilde.init*z.ij*u.ij)
        eta.partial <- -(z.ij/z.sigma^2)*c(1,u.ij)
        B.delta <- -(2/(N*(N - 1)))*(exp.ij/((1 + exp.ij)^2))*theta.partial%*%t(eta.partial)
      }
      else{
        A.delta <- 0
        B.delta <- 0
      }
      A <- A + A.delta
      B <- B + B.delta
    }
  }
  return(list(A_est = A, B_est = B))
}
