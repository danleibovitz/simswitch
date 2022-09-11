

#' full (unadjusted) dataframe generator
#'
#' @param n Number of patients
#' @param stime Number of follow-up times
#' @param prop_trt Proportion of \code{n} that is assigned to treatment
#' @param id_trt Patient ids of patients assigned to treatment
#' @param num_bvar
#' @param bcov
#' @param num_tvar
#'
#' @return a dataframe
#' @export
#'
#' @examples
fd_generator <- function(
    n,
    stime,
    prop_trt,
    id_trt,
    num_bvar,
    bcov,
    num_tvar) {

  # Defend against incorrect argument types
  if(class(n) != "integer") stop("'n' must be of class 'integer'")
  if(class(stime) != "integer") stop("'stime' must be of class 'integer'")
  if(class(prop_trt) != "numeric") stop("'prop_trt' must be of class 'numeric'")
  if(class(id_trt) != "integer") stop("'id_trt' must be of class 'integer'")
  if(class(num_bvar) != "integer") stop("'num_bvar' must be of class 'integer'")
  if(class(bcov) != "data.frame") stop("'bcov' must be of class 'data.frame")
  if(class(num_tvar) != "integer") stop("'num_tvar' must be of class 'integer'")

  # Defend against incorrect argument dimensions
  if(length(n) != 1) stop("'n' must be of length 1")
  if(length(stime) != 1) stop("stime'' must be of length 1")
  if(length(prop_trt) != 1) stop("'prop_trt' must be of length 1")
  if(length(num_bvar) != 1) stop("'num_bvar' must be of length 1")
  if(dim(bcov)[1] != stime*n | dim(bcov)[2] != num_bvar) stop("'bcov' must be of dimension [stime*n]x[num_bvar]")
  if(length(num_tvar) != stime) stop("'num_tvar' must be of length 1")

  ids <- rep(1:n, each=stime) # specify participant ids
  time <- rep(1:stime, n)
  id_trt <- sample(unique(ids), prop_trt*length(unique(ids))) # select participants for treatment
  trt <- rep(NA, length(ids)) # create treatment covariate
  trt[ids %in% id_trt] <- 1 # set treatment covariate to 1 for experimental participants
  trt[is.na(trt)] <- 0
  id_con <- subset(unique(ids), !(unique(ids) %in% id_trt)) # get control group ids

  # if(dep_cov == TRUE){ # generate covariates dependently
  xdat <- data.frame(ids = ids, arm = NA, switch = 0, time = time, treat = trt)
  xdat$arm <- ifelse(xdat$ids %in% id_trt, 1, 0)

  # set random covariates, number equal to num_bvar
  if(missing(bcov)){ # if a baseline covariate matrix is not specified, generate a random one
    bcov <- gen_bcov(num_bvar = num_bvar, diags = 0.4, middle = 1, stime = stime, n = n)
  } else{
    if(dim(bcov)[1] != (n*stime) | dim(bcov)[2] != num_bvar) stop("a pre-specified baseline covariate matrix must have dimension [stime*n]x[num_bvar]")
  }
  bcov_names <- paste0(rep("b", num_bvar), seq.int(1:num_bvar)) # rename bcov columns as b1, b2, etc
  names(bcov) <- bcov_names

  # set empty tcov columns of width num_tvar
  tcov <- as.data.frame(matrix(nrow = dim(xdat[1]), ncol = num_tvar))
  if(num_tvar > 1){
    tcov_names <- c("M", paste0(rep("v", num_tvar-1), seq.int(1:(num_tvar-1)))) #  rename tcov columns as M, v1, v2, etc. First name is always "M"
  }else{
    tcov_names <- c("M")
  }
  names(tcov) <- tcov_names

  fulldat <- cbind(xdat, bcov, tcov) # merge covariate and treatment data
  fulldat$Mtime <- 0 # add variablre for time since beginning of M
  return(fulldat) # set a control group, where there is no time-dependent confounding

}
