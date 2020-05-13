#' Fit the model of \eqn{Z} given\eqn{U}.
#'
#' @param z Vector. The shadow variable (fully observed).
#' @param u Matrix. The covariate matrix excluding shadow variable \eqn{z} (fully observed).
#'
#' @return A list.
#' \describe{
#' \item{eta.est}{The estimates for coefficients in model \eqn{Z} given\eqn{U}.}
#' \item{z.sigma}{The estimates for \eqn{\sigma} in model \eqn{Z} given\eqn{U}.}
#' }
z_u_fit <- function(z, u){
  if(sum(is.na(z)) > 0 | sum(!complete.cases(u)) > 0){
    stop("Remove missing value in covariates.")
  }else{
    z.fit <- glm(as.matrix(z) ~ as.matrix(u),
                 family = gaussian(link = "identity"))
    eta.est <- as.numeric(z.fit$coefficients)
    z.sigma <- sigma(z.fit)
    return(list(eta_est = eta.est,
                z.sigma_est = z.sigma))
  }
}
