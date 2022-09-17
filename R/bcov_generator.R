# generate bcov (baseline covariates)

#' Generate a matrix of baseline covariates
#'
#' @param num_bvar
#' @param diags
#' @param middle
#' @param stime
#' @param n
#'
#' @return
#' @export
#'
#' @examples
gen_bcov <- function(
    num_bvar,
    diags,
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


