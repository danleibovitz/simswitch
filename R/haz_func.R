# define haz() function

haz_func <- function(t, x, betas, b_haz, ncov = 0, ...) {

  # TODO change c(5:(5+ncov)) and all similar structures to calling the actual names in names_bcov, etc.
  # time <- ifelse(t > max(x[["time"]]), data.table::last(x[["time"]]), sapply( t, function(i) min(x[["time"]][x[["time"]] >= i]) ) ) # for vectorized t, get next largest values of time in discrete time sequence of x
  covariates <- x[x[["time"]] == t, c(5:(5+ncov))] # get covariates value, and all covariate values up to ncov. right now, 5 is the index of "treat"
  if(ncov != 0){ # if there are more than 1 covariates...
    exp( rowSums(covariates*betas[betas$time == t, c(2:(2+ncov))]) )*b_haz[t] # calculate exp() of linear dot product. right now, 2 is the index of "treat" in betas matrix
  }
  else{
    exp( covariates*betas[betas$time == t, c(2:(2+ncov))] )*b_haz[t]
    #exp( covariates*rep.int(1,length(t)) )*b_haz[time]
  }
}
