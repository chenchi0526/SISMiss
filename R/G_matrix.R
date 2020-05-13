#' Generate the G matrix for variance estimation in low dimension scenario.
#' @param z Vector. The fully observed observed shadow variable.
#' @param u Matrix. The covariate matrix excluding shadow variable \eqn{z} (fully observed).
#' @param z.sigma_est The estimated standard error of \eqn{z}.
#' @return \eqn{G} matrix for variance estimation in low dimension scenario.

G_matrix <- function(z, u, z.sigma_est){
  N <- nrow(u)
  p <- ncol(u)
  G <- matrix(data = 0, nrow = p+1, ncol = p+1)

  for(i in 1:N){
    u.i.d <- as.numeric(c(1, u[i, ]))
    G <- G - (1/N)*(1/z.sigma_est^2)*u.i.d%*%t(u.i.d)
  }
  return(G)
}

