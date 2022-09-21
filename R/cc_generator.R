

#' Generate matrix of coefficients for generation of time-varying covariates
#'
#' @param numb Number of baseline covariates
#' @param numt Number of time-varying covariates
#'
#' @return A list of 2 data.frames named 'baseline' and 'varying', representing baseline and time varying coefficients
#' for the generation of survival times via a cox-like hazard model.
#'
#' @examples
.cc_generate <- function(
    numb,
    numt) {
  covar_coef <- list(baseline = matrix(sample(1:(num_bvar*num_tvar), num_bvar*num_tvar), ncol = num_tvar),
                     varying = matrix(sample(1:(num_tvar*(num_tvar+1)), num_tvar*(num_tvar+1)), ncol = num_tvar)) # add a coefficient row for the treatment effect
  covar_coef$baseline <- LICORS::normalize(covar_coef$baseline, byrow = FALSE)
  covar_coef$varying <- LICORS::normalize(covar_coef$varying, byrow = FALSE)
  covar_coef$varying[1,] <- -covar_coef$varying[1,] # make the treatment effect protective
  covar_coef$varying <- covar_coef$varying/100
  covar_coef$varying[2:(num_tvar+1),1] <- covar_coef$varying[2:(num_tvar+1),1]*150 # bump up coefficients predictive of secondary baseline
  for(j in 1:num_tvar){
    covar_coef$varying[(j+1),j] <- 1
  }
  }
