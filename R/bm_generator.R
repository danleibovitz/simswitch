

#' Generate matrix of coefficients for survival time generation
#'
#' @param stime Follow-up time for coefficient matrix
#' @param num_bvar Number of baseline coefficients
#' @param num_tvar Number of time-varying coefficients
#' @param bcov_names Names of baseline variables
#' @param tcov_names Names of time-varying variables
#' @param treat_beta Treatment coefficient
#' @param violate Tag for violation of estimation assumptions
#' @param n Number of patients
#'
#' @return A matrix of dimension [n*stime]x[num_tvar + num_bvar + 2]
#' @export
#'
#' @examples
bm_generator <- function(
    stime,
    num_bvar,
    num_tvar,
    bcov_names,
    tcov_names,
    treat_beta,
    violate,
    n) {
  beta.mat <- as.data.frame(matrix(nrow = stime, ncol = num_bvar + num_tvar + 2))
  names(beta.mat) <- c("time", "treat", bcov_names, tcov_names)
  beta.mat$time <- 1:stime
  beta.mat$treat <- treat_beta
  if("All" %in% violate | "RPSFTM" %in% violate){
    beta.mat$treat <- (0.08*sqrt(1:stime)) - 0.9 # make treatment effect time-dependent for violate == RPSFTM
  }
  for(j in 1:(num_bvar + num_tvar)){
    beta.mat[, j + 2] <- runif(1,0, 0.1) # Give non-treatment covariates exclusively small, accelerating impacts on survival
  }
  beta.mat$ids <- NA # This ONLY exists so that simsurv() doesn't complain that the betas dataframe has no 'ids' column, which is unnecessary anyhow

  # replicate beta.mat N times, for dumb reasons...
  beta.mat <- do.call("rbind", replicate(n, beta.mat, simplify = FALSE))
  beta.mat$ids <- rep(1:n, each=stime)
}
