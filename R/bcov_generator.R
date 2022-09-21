
#' Generate a matrix of baseline covariates
#'
#' @param num_bvar Number of baseline variables
#' @param middle Mean of multivariate normal, which generates covariates
#' @param stime Length of follow-up time
#' @param n Number of patients
#'
#' @return A data.frame with dimension [n*stime]x[num_bvar]
#' @export
#'
#' @examples
bcov_generator <- function(
    num_bvar,
# TODO (remove this arg?:)    diags,
    middle,
    stime,
    n) {

  covariance <- matrix(data = stats::runif(num_bvar^2, 0.1, 0.5), nrow = num_bvar)
  diag(covariance) <- 0.4
  covariance <- t(covariance) %*% covariance # force positive definite
  bcov <- as.data.frame((MASS::mvrnorm(n, rep.int(middle,num_bvar), covariance))) # generate baseline covariates for each patient
  bcov[,1] <- ifelse(bcov[,1] > middle, 1, 0) # cast as a binary, representing some kind of predisposition to death AND switching
  bcov <- bcov[rep(seq_len(nrow(bcov)), each = stime), ] # repeat each row of bcov stime times
}


