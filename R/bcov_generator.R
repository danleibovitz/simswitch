
#' Generate a matrix of baseline covariates
#'
#' @param num_bvar Number of baseline variables
#' @param middle Mean of multivariate normal, which generates covariates
#' @param stime Length of follow-up time
#' @param n Number of patients
#' @param diags Covariance between all baseline covariate generation distributions
#'
#' @return A data.frame with dimension [n*stime]x[num_bvar]
#' @export
#'
#' @examples
bcov_generator <- function(
    num_bvar,
    diags = 0.4,
    middle,
    stime,
    n) {

  # defend against incorrect argument classes
  if(class(num_bvar) != "numeric") stop("'num_bvar' must be of class 'numeric'")
  if(class(diags) != "numeric") stop("'diags' must be of class 'numeric'")
  if(class(middle) != "numeric") stop("'middle' must be of class 'numeric'")
  if(class(stime) != "numeric") stop("'stime' must be of class 'numeric'")
  if(class(n) != "numeric") stop("'n' must be of class 'numeric'")

  # defend against incorrect dimension
  if(length(num_bvar) != 1) stop("'num_bvar' must be of length 1")
  if(length(diags) != 1) stop("'diags' must be of length 1")
  if(length(middle) != 1) stop("'middle' must be of length 1")
  if(length(stime) != 1) stop("'stime' must be of length 1")
  if(length(n) != 1) stop("'n' must be of length 1")

  # defend against incorrect argument range
  if(num_bvar <1) stop("'num_bvar' must be a positive integer")
  if(diags <= 0 | diags >= 1) stop("'diags' must be between 0 and 1")
  if(stime < 1) stop("'stime' must be a positive integer")
  if(n < 1) stop("'n' must be a positive integer")


  covariance <- matrix(data = stats::runif(num_bvar^2, 0.1, 0.5), nrow = num_bvar)
  diag(covariance) <- diags
  covariance <- t(covariance) %*% covariance # force positive definite
  bcov <- as.data.frame((MASS::mvrnorm(n, rep.int(middle,num_bvar), covariance))) # generate baseline covariates for each patient
  bcov[,1] <- ifelse(bcov[,1] > middle, 1, 0) # cast as a binary, representing some kind of predisposition to death AND switching
  bcov <- bcov[rep(seq_len(nrow(bcov)), each = stime), ] # repeat each row of bcov stime times
}


