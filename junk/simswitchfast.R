# main function

## TODO
# re-format styling
# move features out that could be written as own functions. these can be internal, not exported.
# make K-M curve for each dataset, and make K-M curves

# load libraries
library(simsurv)
library(MASS)
  # library(coxed)
library(data.table)
library(survival)
library(ggplot2)
library(survminer)
library(rpsftm)
  # library(LICORS)
library(ipcwswitch)
library(boot)
library(tidyr)
library(purrr)

weihaz <- function(x, shape, scale){
  (shape/scale)*(x/scale)^(shape-1)
}

# Discrete-time survival time generator
# TODO I think... the survive() function is mixing up the order of the dataframe.
survive <- function(x, hazard, betas, ncov, stime, idvar, ids,
                    b_haz){
  n <- length(ids)
  store <- ids # create repository of patients with unobserved events
  df <- data.frame(ids = ids, eventtime = rep.int(0, n), status = rep.int(0, n)) # create response dataframe
  i <- 1 # set iterator

  # for all time windows: probability of failure in that window is hazard produced by haz() function
  # TODO replace for-loop with apply() df$eventtime <- apply(2:stime, function(x) ifelse(df$status[df$time == x-1] == 1))
  while(any(df$eventtime == 0) & i <= stime){ # TODO should this be a while loop, and stop when all events are observed?

    # randomly assign event to remaining patients
    # TODO very important -- call to hazard() must have appropriately indexed x and betas values
    df$eventtime[df$ids %in% store] <- ifelse(rbinom(n = length(store), size = 1, prob = hazard(t=i, x=x[x$ids %in% store,], betas = betas[betas$ids %in% store,], b_haz = b_haz, ncov = ncov)) == 1, i,0) # assign current time to for patients who observe in this window

    # update store
    store <- df$ids[df$eventtime == 0]
    i <- i + 1
  }

  # update statuses
  df$status[df$eventtime != 0] <- 1
  # apply administrative censoring to any remaining patients, and leave remaining statuses as 0
  df$eventtime[df$eventtime == 0] <- stime

  return(df)
}

#' Simulate TS data and apply treatment-effect-estimating methods
#'
#' @param add_tvar
#' @param b_allowance
#' @param b_haz
#' @param b_mag
#' @param b_scale
#' @param b_shape
#' @param bcov
#' @param beta.mat
#' @param bootrep
#' @param cens_flag
#' @param covar_coef
#' @param dep_func
#' @param haz
#' @param hide_tvar
#' @param ipcw_robust
#' @param m_allowance
#' @param m_inflation
#' @param m_fidelity
#' @param m_hard
#' @param m_haz
#' @param m_mag
#' @param m_scale
#' @param m_shape
#' @param n
#' @param num_bvar
#' @param num_tvar
#' @param prop_cens
#' @param prop_cens_allowance
#' @param prop_cont_event
#' @param prop_switch
#' @param prop_trt
#' @param prop_trt_event
#' @param recens
#' @param rerun_lim
#' @param s_allowance
#' @param s_haz
#' @param s_mag
#' @param s_scale
#' @param s_shape
#' @param stime
#' @param switch_coef
#' @param t_allowance
#' @param t_mag
#' @param treat_beta
#' @param tse_dist
#' @param verbose
#' @param violate
#'
#' @return
#' @export
#'
#' @examples
simswitch <- function(add_tvar, b_allowance, b_haz, b_mag, b_scale, b_shape, bcov, beta.mat, bootrep,
                      cens_flag, covar_coef, dep_func, haz, hide_tvar, ipcw_robust,
                      m_allowance, m_inflation, m_fidelity, m_hard, m_haz, m_mag, m_scale, m_shape,
                      n, num_bvar, num_tvar, prop_cens, prop_cens_allowance, prop_cont_event, prop_switch,
                      prop_trt, prop_trt_event, recens, rerun_lim, s_allowance, s_haz, s_mag, s_scale,
                      s_shape, seed, stime, switch_coef, t_allowance, t_mag, treat_beta, tse_dist, unfix, verbose, violate){


  # TODO
  # add option for switching from experimental to control

  if(missing(verbose)){
    verbose <- 2
  }
  if(verbose > 1){
    print("Setting parameters...")
  }
  # set ipcw_robust. If TRUE, we skip bootstrapping for IPCW. Bootstrapped IPCW may not be symmetrical, and may be more accurate, but takes hella long
  if(missing(ipcw_robust)){
    ipcw_robust <- TRUE
  }

  if(missing(tse_dist)){ # alternatives are weibull, lognormal, etc.
    tse_dist <- "loglogistic"
  }
  # set bootrep param
  if(missing(bootrep)){
    bootrep <- 1000
  }
  # set assumption violation flag
  # TODO implement automatic assumption violation. For RPSFTM, non-constant treatment effect. For TSE, switching after secondary baseline. For IPCW, unmeasured confounding
  if(missing(violate)){
    violate <- "None"
  } else{
    if(!( "All" %in% violate | "None" %in% violate | all(violate %in% c("RPSFTM", "TSE", "IPCW")) )){
      stop("Violate must be set to All, None, or a subset of RPSFTM, TSE and IPCW")
    }
  }

  # set seed for seed option
  if(!missing(seed)){set.seed(seed)}

  # set random censoring flag. One of no censoring, random censoring, non-random censoring
  if(missing(cens_flag)){
    cens_flag <- "Random"
  }else{
    if(!(cens_flag %in% c("None", "Random", "Nonrandom"))) stop("cens_flag must be one of None, Random, or Nonrandom")
  }
  if(missing(prop_cens) & cens_flag != "None"){
    prop_cens <- 0.1 # default censoring proportion is 0.3 across both groups. This is pre-administrative censoring.
  }
  if(missing(prop_cens_allowance)){
    prop_cens_allowance <- 0.1
  }
  if(prop_cens < 0 | prop_cens >= 1) stop("prop_cens must be in range [0,1)")

  if(missing(recens)){
    recens <- TRUE # set flag for recensoring of RPSFTM and Two-Stage Model
  }
  # set dependent covariates flag
  # if(missing(dep_cov)){
  #   dep_cov <- TRUE # keep default FALSE while testing dependency generation
  #   }
  if(missing(unfix)){
    unfix <- as.character(c())
  }
  if(!is.character(unfix) | any(! unfix %in% c("B", "M", "S", "T"))){
    stop("unfix must be of type character, and only contain \"B\", \"M\", \"S\", or \"T\" for
         Baseline, Metastatic disease, Switch, or Treatment")
  }
  # set survival time, i.e., administrative censoring time
  if(missing(stime)){
    stime <- 100
  }
  if(missing(num_bvar)){
    num_bvar <- 3
  }
  num_bvar <- as.integer(num_bvar)
  if(missing(num_tvar)){
    num_tvar <- 3
  }
  num_tvar <- as.integer(num_tvar)
  if(num_tvar < 1) stop("Must have at least 1 time-varying covariate")
  # add option for hiding relevent covariates from switching modeling
  if(missing(hide_tvar)){
    hide_tvar <- 0
  }
  if("All" %in% violate | "IPCW" %in% violate){
    hide_tvar <- 2 # exclude 2 tvar
  }
  if(hide_tvar > (num_tvar-1) ) stop("can't hide more time-varying covariates than there are time-varying covariates, and time-varying cov M must always be kept.")
  # add option for adding irrelevant covariates
  if(missing(add_tvar)){
    add_tvar <- 0
  }

  # set number of participants
  if(missing(n)){
    n <- 400
  }
  # set giant while loop flag
  rerun <- TRUE
  # set while loop limit
  if(missing(rerun_lim)){
    rerun_lim <- 200
  }
  # set proportion of participants on experimental treatment
  if(missing(prop_trt)){
    prop_trt <- 0.5
  }
  if(missing(prop_switch)){
    prop_switch <- 0.5
  }# set proportion of control participants to switch
  if(missing(switch_coef)){
    switch_coef <- c(runif(num_bvar, 0.1, 0.3), runif(num_tvar, 0.3, 1)) # default of switch_coef log hazard ratios. baseline switch coefs are smaller.
  }
  if(length(switch_coef) != sum(num_bvar, num_tvar)) stop("the switching hazard coefficients must be of the same length as all covariates")


  # Create patient id, randomization, and discrete time structure
  ids <- rep(1:n, each=stime) # specify participant ids
  time <- rep(1:stime, n)
  id_trt <- sample(unique(ids), prop_trt*length(unique(ids))) # select participants for treatment
  trt <- rep(NA, length(ids)) # create treatment covariate
  trt[ids %in% id_trt] <- 1 # set treatment covariate to 1 for experimental participants
  trt[is.na(trt)] <- 0
  id_con <- subset(unique(ids), !(unique(ids) %in% id_trt)) # get control group ids
  if(verbose > 1){
    print("Building data frames...")
  }
  # if(dep_cov == TRUE){ # generate covariates dependently
  xdat <- data.frame(ids = ids, arm = NA, switch = 0, time = time, treat = trt)
  xdat$arm <- ifelse(xdat$ids %in% id_trt, 1, 0)

  # TODO disease progression can be a OS * beta(a, b) and somehow dependent on baseline covars
  # set random covariates, number equal to num_bvar
  if(missing(bcov)){ # if a baseline covariate matrix is not specified, generate a random one
    bcov <- as.data.frame(abs(MASS::mvrnorm(n, rep.int(0,num_bvar), diag(0.4, num_bvar, num_bvar)))) # generate baseline covariates for each patient
    bcov <- bcov[rep(seq_len(nrow(bcov)), each = stime), ] # repeat each row of bcov stime times
    } else{
      if(dim(bcov)[1] != (n*stime) | dim(bcov)[2] != num_bvar) stop("a pre-specified baseline covariate matrix must have length equal to (number of patients)*(stime) and width equal to num_bvar")
  }
  bcov_names <- paste0(rep("b", num_bvar), seq.int(1:num_bvar)) # rename bcov columns as b1, b2, etc
  names(bcov) <- bcov_names

  # set empty tcov columns of width num_tvar
  tcov <- as.data.frame(matrix(nrow = dim(xdat[1]), ncol = num_tvar))
  tcov_names <- c("M", paste0(rep("v", num_tvar-1), seq.int(1:(num_tvar-1)))) #  rename tcov columns as M, v1, v2, etc. First name is always "M"
  names(tcov) <- tcov_names

  fulldat <- cbind(xdat, bcov, tcov) # merge covariate and treatment data
  fulldat$Mtime <- 0 # add variablre for time since beginning of M
  fulldat_cont <- fulldat # set a control group, where there is no time-dependent confounding

  # set a default covar_coef list. If user defines their own dep_func, it must call covar_coef argument even if it doesnt use it
  if(missing(covar_coef)){
    covar_coef <- list(baseline = matrix(sample(1:(num_bvar*num_tvar), num_bvar*num_tvar), ncol = num_tvar),
                       varying = matrix(sample(1:(num_tvar*(num_tvar+1)), num_tvar*(num_tvar+1)), ncol = num_tvar)) # add a coefficient row for the treatment effect
    covar_coef$baseline <- LICORS::normalize(covar_coef$baseline, byrow = FALSE)
    # sweep(covar_coef$baseline, 2, colSums(covar_coef$baseline), `/`)
    covar_coef$varying <- LICORS::normalize(covar_coef$varying, byrow = FALSE)
    # sweep(covar_coef$varying, 2, colSums(covar_coef$varying), `/`)
    covar_coef$varying[1,] <- -covar_coef$varying[1,] # make the treatment effect protective
    # TODO try to make more orderly...
    covar_coef$varying <- covar_coef$varying/200
    for(j in 1:num_tvar){
      covar_coef$varying[(j+1),j] <- 1
    }
  }

  # set dep_func (dependent function for generating covariates). must always return matrix of length n, width num_tvar + 1. at least first index of num_tvar is binary
  if(missing(dep_func)){
    dep_func <- function(dat, window, base_var, time_var, covar_coef, m_haz){ # should return a binary variable for M, and continuous for all other tvar
      if(window == 1){# if window is 1, only a function of baseline.
        retval <- (as.matrix(dat[dat$time == window, names(dat) %in% base_var]) %*% covar_coef$baseline) + MASS::mvrnorm(n = sum(dat$time == window), rep.int(0,num_tvar), diag(0.05,nrow = num_tvar))
        retval[,1] <- 0 # Set initial M values to 0
        retval <- cbind(retval, rep.int(0, length(retval[,1]))) # set initial Mtimes
        return(retval)
      }else{# if window is not 1, a function of baseline, previous tvar and previous treat
        retval <- (as.matrix(dat[dat$time == (window-1), names(dat) %in% c("treat", time_var)]) %*% covar_coef$varying) + MASS::mvrnorm(n = sum(dat$time==(window-1)), rep.int(0, num_tvar), diag(0.05, nrow = num_tvar))
        # retval[,1] <- plogis(4)
        # change column of metastatic disease to binary, based upon hazard function
        retval[,1] <- ifelse(dat$M[dat$time == i-1] == 1, 1, 0) # if the previous window M is 1, continue M
        retval[retval[,1] == 0, 1] <- (rbinom( n = length(retval[, 1]), size = 1,
                                               prob = 1 - exp(-exp(log(m_haz[i]) +
                                                                     as.matrix(dat[dat$time == i-1, names(dat) %in% c("treat", tcov_names)]) %*% covar_coef$varying[,1])) ))[retval[,1] == 0] # randomly assign the treat variable with probability
        retval <- cbind(retval, rep.int(0, length(retval[,1])))
        retval[,num_tvar + 1] <- ifelse(retval[,1] == 0, 0,
                                        ifelse(dat$M[dat$time == i-1] == 0, 1,
                                               dat$Mtime[dat$time == i-1]+1))
        # retval$M <- beta # replace M (metastatic disease) with a
        # fulldat$M
        return(retval)
      }
    }
  }

  # set baseline switch hazard function. lam is set as lambda of exp distribution with expected value of prop_switch by stime
  if(missing(s_haz)){
    s_haz <- weihaz(1:stime, 2, 0.7*stime)
  }
  if(!missing(s_shape) & !missing(s_scale)){
    s_haz <- weihaz(1:stime, s_shape, s_scale)
  }

  ## this section is a bit tedious. We have to set up parameters to iteratively search for the correct switching proportion
  switch_iter <- 0 # how many times have we searched?
  if(missing(s_allowance)){
    s_allowance <- 0.1
  }# how far from the proportion of switching requested is acceptable?
  s_direc <- 0 # was the last attempt too low or too high?
  prev_s_direc <- 0 # set holder variable for the previous s_direc
  if(missing(s_mag)){
    s_mag <- 0.5
  }# by what factor should we adjust the baseline hazard? must be:
  if(s_mag >= 1 | s_mag <= 0) stop("s_mag must be between 0 and 1, exclusive")

  if(missing(m_allowance)){
    m_allowance <- 0.1
  }
  m_direc <- 0
  prev_m_direc <- 0
  if(missing(m_mag)){
    m_mag <- 2
  }
  if(missing(m_inflation)){
    m_inflation <- 2
  }
  if(missing(m_fidelity)){
    m_fidelity <- 0.2
  }# represents the proportion of stime away from first M that switch can occur
  if(missing(m_hard)){
    m_hard <- TRUE
  }# represents weather switch can happen only after M, or if we don't care
  if(violate == "All" | "TSE" %in% violate){
    m_hard <- FALSE
  }
  if(!(m_hard %in% c(TRUE, FALSE))) stop("m_hard must be Boolean")
  if(!missing(m_shape) & !missing(m_scale)){
    m_haz <- weihaz(1:stime, m_shape, m_scale)
  }
  if(missing(m_haz)){ # TODO if m_haz IS specified, we need to provide some other way of adjusting m_haz with respect to current_m_prop
    m_shape <- 2.8
    m_scale <- 0.7*stime
    m_haz <- weihaz(1:stime, m_shape, m_scale)
  }


  if(missing(b_haz)){
    b_haz <- weihaz(1:stime, 1, stime) # default is an exponential dist with lambda = stime
  }
  # if(missing(lambdas)){
  #   lambdas <- 1/stime # default scale parameter of exponential distribution is 1/stime
  # }
  if(!missing(b_shape) & !missing(b_scale)){
    b_haz <- weihaz(1:stime, b_shape, b_scale)
  }
  if(missing(b_allowance)){
    b_allowance <- 0.1
  }
  if(missing(b_mag)){
    b_mag <- 0.5
  }
  b_direc <- 0
  prev_b_direc <- 0

  # TODO here set observed event parameters and threshhold
  if(missing(prop_trt_event)){ # if proportion of trt pts is out of window defined by b_allowance, we adjust b_haz by b_mag
    prop_trt_event <- 0.25
  }
  if(missing(prop_cont_event)){ # if proportion of control pts is out of window defined by t_allowance, we adjust treatment coefficients by t_mag
    prop_cont_event <- min(1, 1.75*prop_trt_event)
  }
  if(missing(t_allowance)){
    t_allowance <- 0.1
  }
  if(missing(t_mag)){
    t_mag <- 0.5
  }
  t_direc <- 0
  prev_t_direc <- 0

  if(verbose > 1){
    print("Setting survival coefficients...")
  }
  if(missing(treat_beta)){
    treat_beta <- -1 # default treatment coef is -1
  }
  if(missing(beta.mat)){
    beta.mat <- as.data.frame(matrix(nrow = stime, ncol = num_bvar + num_tvar + 2))
    names(beta.mat) <- c("time", "treat", bcov_names, tcov_names)
    beta.mat$time <- 1:stime
    beta.mat$treat <- treat_beta
    if(violate == "All" | "RPSFTM" %in% violate){
      beta.mat$treat <- exp((0.35*log(1:stime)) - 1.5) # make treatment effect time-dependent for violate == RPSFTM
    }
    for(j in 1:(num_bvar + num_tvar)){
      beta.mat[, j + 2] <- runif(1,0.1, 0.2) # Give non-treatment covariates exclusively small, accelerating impacts on survival
    }
    beta.mat$ids <- NA # This ONLY exists so that simsurv() doesn't complain that the betas dataframe has no 'ids' column, which is unnecessary anyhow

    # replicate beta.mat N times, for dumb reasons...
    beta.mat <- do.call("rbind", replicate(n, beta.mat, simplify = FALSE))
    beta.mat$ids <- rep(1:n, each=stime)
  }

  if(verbose > 1){
    print("Setting baseline hazard function...")
  }
  if(missing(haz)){
    haz <- function(t, x, betas, b_haz, ncov = 0, ...) {

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
  }

  if(verbose > 1){
    print("Filling data frames")
  }
  # TODO giant while loop should begin here. At this point, we have all parameters set, and a fulldat dataframe with no time
  # varying covariates, no switching and no secondary baseline
  while(rerun){ # iteratively update switching hazard function until we get the right proportion. The first term in this check is the number of pts whose final treatment indicator

    for(i in 1:stime){
      # set covars
      fulldat[fulldat$time == i, names(fulldat) %in% c(tcov_names, "Mtime")] <-
        dep_func(dat = fulldat, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz)
      # set treatment indicator
      if(i != 1){ # if its not the first time window
        fulldat$treat[fulldat$time == i] <- ifelse(fulldat$treat[fulldat$time == i-1] == 1, 1, 0) # if the previous window treat is 1, continue treatment
      }
      # if treatment has not yet begun, probability of begining in the next window is a hazard function
      if(m_hard){
        fulldat$treat[fulldat$time == i & fulldat$treat == 0 & fulldat$Mtime > 0 & fulldat$Mtime <= ceiling(m_fidelity*stime)] <-
          rbinom( n = length(fulldat$treat[fulldat$time == i & fulldat$treat == 0 & fulldat$Mtime > 0 & fulldat$Mtime <= ceiling(m_fidelity*stime)]), size = 1,
                  prob = 1 - exp(-exp(log(s_haz[i]) +
                                        as.matrix(fulldat[fulldat$time == i & fulldat$treat == 0, names(fulldat) %in% c(bcov_names, tcov_names)]) %*% switch_coef))) # randomly assign the treat variable with probability
      }else{
        fulldat$treat[fulldat$time == i & fulldat$treat == 0] <- rbinom( n = length(fulldat$treat[fulldat$time == i & fulldat$treat == 0]), size = 1,
                                                                         prob = 1 - exp(-exp(log(s_haz[i]) +
                                                                                               as.matrix(fulldat[fulldat$time == i & fulldat$treat == 0, names(fulldat) %in% c(bcov_names, tcov_names)]) %*% switch_coef))) # randomly assign the treat variable with probability
      }

    }

    # Generate fulldat_cont (the control dataset) by blocking switching
    fulldat_cont$treat <- fulldat_cont$arm
    for(i in 1:stime){
      # set covars
      fulldat_cont[fulldat_cont$time == i, names(fulldat_cont) %in% c(tcov_names, "Mtime")] <-
        dep_func(dat = fulldat_cont, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz)
    }

    #  else{ # generate covariates independently. much faster.
    #     id_switch <- sample(id_con, prop_switch*length(id_con)) # select participants for switch
    #     # collect into dataframe
    #     xdat <- data.frame(ids = ids, arm = NA, switch = NA, time = time, treat = trt)
    #     xdat$arm <- ifelse(xdat$ids %in% id_trt, 1, 0)
    #     xdat$switch <- ifelse(xdat$ids %in% id_switch, 1, 0)
    #     # set random time stime to trt == 1
    #     for(i in id_switch){
    #       start <- sample(1:stime, 1)
    #       xdat$treat[xdat$ids == i & xdat$time >= start] <- 1
    #     }
    #
    #     # set random covariates, number equal to num_bvar
    #     tcov <- mvrnorm(n*stime, rep.int(0,num_bvar), diag(1, num_bvar, num_bvar))
    #
    #     fulldat <- cbind(xdat, tcov) # merge covariate and treatment data
    #     # TODO rename covars. bi for baseline covars i and vi for tvar covars i
    # }


    # TODO Standardize all covariates ?
    ############
    # fulldat[fulldat$arm == 0, names(fulldat) %in% tcov_names[tcov_names != "M"]] <- scale(fulldat[fulldat$arm == 0, names(fulldat) %in% tcov_names[tcov_names != "M"]])
    # fulldat[fulldat$arm == 1,names(fulldat) %in% tcov_names[tcov_names != "M"]] <- scale(fulldat[fulldat$arm == 1,names(fulldat) %in% tcov_names[tcov_names != "M"]])
    # fulldat_cont[fulldat_cont$arm == 0, names(fulldat_cont) %in% tcov_names[tcov_names != "M"]] <- scale(fulldat_cont[fulldat_cont$arm == 0, names(fulldat_cont) %in% tcov_names[tcov_names != "M"]])
    # fulldat_cont[fulldat_cont$arm == 1, names(fulldat_cont) %in% tcov_names[tcov_names != "M"]] <- scale(fulldat_cont[fulldat_cont$arm == 1, names(fulldat_cont) %in% tcov_names[tcov_names != "M"]])
    ############


    ######################
    # Are the survival times of the confounded and unconfounded datasets probabalistically equivalent? visual inspection
    # fulldat$dat <- 0
    # fulldat_cont$dat <- 1
    # testdat <- rbind(fulldat[fulldat$arm == 1, ], fulldat_cont[fulldat_cont$arm ==1, ])
    # testdat$ids <- beta.mat$ids
    # testsdat <- simsurv( x = testdat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), maxt = 100, idvar = "ids", ids = unique(fulldat$ids),
    #                  b_haz = rep(0.01, stime))
    # testdat <- merge(testdat, testsdat, by = "ids")
    # survminer::ggsurvplot(
    #   fit = survfit(survival::Surv(eventtime, status) ~ dat, data = testdat[testdat$time == 1,]),
    #   conf.int = TRUE)
    #########################

    if(verbose > 1){
      print("Generating confounded survival times...")
    }
    # if(b_haz == "exponential"){
    #   #
    #   # TODO simsurv currently wont accept beta.mat as a dataframe. Here we are recasting test.mat as the first row of beta.mat
    #   # TODO this shouldnt be necessary, and we need to be able to pass beta.mat as a full dataframe.
    #   # TODO Currently, when we pass beta.mat as a dataframe, eventtime values are EXTREMELY small.
    #   test.mat <- as.numeric(beta.mat[1, !names(beta.mat) %in% c("time", "ids")])
    #   names(test.mat) <- names(beta.mat)[!names(beta.mat) %in% c("time", "ids")]
    #   #
    #   sdat <- simsurv( x = fulldat, dist = "exponential", lambdas = lambdas, betas = test.mat, maxt = stime)
    #   names(sdat) <- c("ids", "eventtime", "status") # Strangely, simsurv returns matrices with different column names depending on the use of the 'dist' argument. We correct that here.
    # }else{
      sdat <- survive( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
                     b_haz = b_haz)
    # }

    if(verbose > 1){
      print("Generating un-confounded survival times...")
    }
    # if(b_haz == "exponential"){
    #   #
    #   # TODO simsurv currently wont accept beta.mat as a dataframe. Here we are recasting test.mat as the first row of beta.mat
    #   # TODO this shouldnt be necessary, and we need to be able to pass beta.mat as a full dataframe.
    #   # TODO Currently, when we pass beta.mat as a dataframe, eventtime values are EXTREMELY small.
    #   test.mat <- as.numeric(beta.mat[1, !names(beta.mat) %in% c("time", "ids")])
    #   names(test.mat) <- names(beta.mat)[!names(beta.mat) %in% c("time", "ids")]
    #   #
    #   sdat_cont <- simsurv( x = fulldat_cont[fulldat_cont$arm == 0,], dist = "exponential", lambdas = lambdas, betas = test.mat, maxt = stime)
    #   names(sdat_cont) <- c("ids", "eventtime", "status") # Strangely, simsurv returns matrices with different column names depending on the use of the 'dist' argument. We correct that here.
    # }else{
      sdat_cont <- survive( x = fulldat_cont[fulldat_cont$arm == 0,], hazard = haz, betas = beta.mat[fulldat_cont$arm == 0,], ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat_cont$ids[fulldat_cont$arm==0]),
                          b_haz = b_haz)

    # }

    if(verbose > 1){
      print("Adding censoring...")
    }
    # merge datasets:
    if(switch_iter == 0){ # if its the first iteration
      fulldat <- merge(fulldat, sdat, by = "ids")
      fulldat_cont <- merge(fulldat_cont, sdat_cont, by = "ids", all = TRUE)
      # replace experimental group in unconfounded dataset with experimental group in confounded dat. Theoretically, these are generated identically.
      fulldat_cont[fulldat_cont$arm == 1, ] <- fulldat[fulldat$arm == 1, ]
    }else{
      fulldat <- fulldat[, !names(fulldat) %in% c("eventtime", "status")] # first remove old sdat
      fulldat <- merge(fulldat, sdat, by = "ids")

      fulldat_cont <- fulldat_cont[, !names(fulldat_cont) %in% c("eventtime", "status")] # first remove old sdat
      fulldat_cont <- merge(fulldat_cont, sdat_cont, by = "ids", all = TRUE)
      # replace experimental group in unconfounded dataset with experimental group in confounded dat. Theoretically, these are generated identically.
      fulldat_cont[fulldat_cont$arm == 1, ] <- fulldat[fulldat$arm == 1, names(fulldat) %in% names(fulldat_cont)]
    }


    # censor at prespecified, non-administrative censoring time. Censoring will have been generated independently of covars, or dependent on covars
    if(cens_flag == "Random"){
      fulldat$cens <- stime # set default censoring time
      fulldat_cont$cens <- stime # set default censoring time

      cens_ids <- sample(unique(fulldat$ids), ceiling(prop_cens*n)) # sample prop_cens of ids for censoring
      rand_cens <- rbeta(n = length(cens_ids), shape1 = 2, shape2 = 1.5) # censoring times are drawn from a beta dist
      # rand_cens <- ceiling(rand_cens)
      rand_cens <- rep(rand_cens, each=stime) # spread to length of dataset
      fulldat$cens[fulldat$ids %in% cens_ids] <- ceiling(rand_cens*fulldat$eventtime[fulldat$ids %in% cens_ids])
      fulldat_cont$cens[fulldat_cont$ids %in% cens_ids] <- ceiling(rand_cens*fulldat_cont$eventtime[fulldat_cont$ids %in% cens_ids])

      # censor observations. Take minimum of eventtime and cens, and change status to 0 if eventtime == cens
      fulldat$eventtime <- pmin(fulldat$eventtime, fulldat$cens)
      fulldat$status <- ifelse(fulldat$eventtime == fulldat$cens, 0, fulldat$status)
      fulldat_cont$eventtime <- pmin(fulldat_cont$eventtime, fulldat_cont$cens)
      fulldat_cont$status <- ifelse(fulldat_cont$eventtime == fulldat_cont$cens, 0, fulldat_cont$status)


    } else if(cens_flag == "Nonrandom"){
      # TODO nonrandom censoring function!
    }


    # Does the patient (control and experimental) observe secondary baseline (M)?
    # TODO i think pts who observe NO M at all are being given NA here
    fulldat$secbase_observed <- as.numeric(rep(sapply(unique(fulldat$ids), function(x)
      ifelse(sum(fulldat$M[fulldat$ids == x]) > 0 & fulldat$time[fulldat$Mtime == 1 & fulldat$ids == x] <= fulldat$eventtime[fulldat$Mtime == 1 & fulldat$ids == x], # id has at least some m, and time at first m is less than event time
             1,
             0)), each = stime))

    # Does the patient (only control) observe switch _before eventtime_?
    fulldat$switch_status <- rep(sapply(unique(fulldat$ids), function(x)
      ifelse(sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime]) > 0 & fulldat$arm[fulldat$ids == x][1] == 0,
             1,
             0)), each = stime)

    switch_iter <- switch_iter + 1 # iterate search indicator
    if(switch_iter > rerun_lim) stop("Your covariate model failed to converge. Try different hazards and/or coefficients")

    rerun <- FALSE # prempt a rerun, unless the following conditions are met

    # if not enough M occurences, adjust M hazard
    # TODO for now, M occurence is only being adjusted upward. do we want to give it a within-window type adjustment?
    current_m_prop <- length(fulldat$ids[fulldat$arm == 0 & fulldat$Mtime == 1 & fulldat$time <= fulldat$eventtime])
    if(! "M" %in% unfix){
      if(current_m_prop < min(n*prop_trt, prop_switch*prop_trt*n*m_inflation)){ # current_m_prop is proportion of control arm that observes M. if it's less than the min of (the whole control arm, the switchers*m_inflation ), increase m_haz and rerun
        m_haz <- m_haz*m_mag # change scale of m related haz, and rerun
        rerun <- TRUE
      }
    }
    print(paste("proportion M: ", current_m_prop))
    print(m_haz[1:7])

    # adjust switching proportion
    current_s_prop <- sum(fulldat$switch_status[fulldat$time == 1 & fulldat$arm == 0])/length(id_con) # get current switch proportion
    if(! "S" %in% unfix){
      if(abs(current_s_prop - prop_switch) > s_allowance){ # If we're outside the window
        prev_s_direc <- s_direc # save most recent s_direc value
        s_direc <- ifelse(current_s_prop > prop_switch, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
        if((prev_s_direc == -1 & s_direc == 1) | (prev_s_direc == 1 & s_direc == -1)){ # if we overshot the allowance window, and must backtrack
          s_mag <- s_mag/2 # halve switch magnitude if weve crossed over the window
        }
        s_haz <- s_haz + s_direc*s_mag*s_haz # modify the switch hazard by a factor of s_mag in the direction of s_direc
        rerun <- TRUE
      }
    }
    print(paste("proportion of switch: ",current_s_prop))
    print(paste("s_mag: ", s_mag))
    print(paste("s_direc: ", s_direc))
    print(s_haz[1:7])

    # adjust control group event occurence
    current_b_prop <- sum(fulldat$status[fulldat$time == 1 & fulldat$arm == 0])/length(id_con) # get current switch proportion
    if(! "B" %in% unfix){
      if(abs(current_b_prop - prop_cont_event) > b_allowance){ # If we're outside the window
        # if(b_haz == "exponential"){
        #   prev_b_direc <- b_direc # save most recent b_direc value
        #   b_direc <- ifelse(current_b_prop > prop_cont_event, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
        #   if((prev_b_direc == -1 & b_direc == 1) | (prev_b_direc == 1 & b_direc == -1)){ # if we overshot the allowance window, and must backtrack
        #     b_mag <- b_mag/2 # halve switch magnitude if weve crossed over the window
        #   }
        #   lambdas <- lambdas - b_direc*b_mag*lambdas # modify the switch hazard by a factor of s_mag in the direction of s_direc
        #   if(lambdas <= 0) stop("The scale parameter for the simsurv distribution is 0 or negative; it must be positive.")
        #   rerun <- TRUE
        # }else{
          prev_b_direc <- b_direc # save most recent b_direc value
          b_direc <- ifelse(current_b_prop > prop_cont_event, -1, 1) # identify the s_direc, i.e., in what direction the baseline hazard must be adjusted
          if((prev_b_direc == -1 & b_direc == 1) | (prev_b_direc == 1 & b_direc == -1)){ # if we overshot the allowance window, and must backtrack
            b_mag <- b_mag/2 # halve switch magnitude if weve crossed over the window
          }
          b_haz <- b_haz + b_direc*b_mag*b_haz # modify the switch hazard by a factor of s_mag in the direction of s_direc
          rerun <- TRUE
        # }
      }
    }

    # adjust trt group event occurence
    current_t_prop <- sum(fulldat$status[fulldat$time == 1 & fulldat$arm == 1])/length(id_trt) # get current switch proportion
    if(! "T" %in% unfix){
      if(abs(current_t_prop - prop_trt_event) > t_allowance){ # If we're outside the window
        prev_t_direc <- t_direc # save most recent b_direc value
        t_direc <- ifelse(current_t_prop > prop_trt_event, 1, -1) # identify the t_direc, i.e., in what direction the treatment effect must be adjusted
        if((prev_t_direc == -1 & t_direc == 1) | (prev_t_direc == 1 & t_direc == -1)){ # if we overshot the allowance window, and must backtrack
          t_mag <- t_mag/2 # halve switch magnitude if weve crossed over the window
        }
        beta.mat$treat <- beta.mat$treat + t_direc*t_mag*beta.mat$treat # modify the switch hazard by a factor of s_mag in the direction of s_direc
        rerun <- TRUE
      }
    }

  }

  if(verbose > 1){
    print("Performing naive method estimates...")
  }
  # TODO make giant loop. catch here if number of observed secbase/M, switch, or  is too low, and adjust M-hazard, switch-hazard, or b_haz (baseline hazard)
  #######
  # if(observed M is too low or high){
  #   adjust M_haz
  # }
  # if(observed switch is too low or high){
  #   adjust s_haz
  # }
  # if(observed death is too low or high){
  #   adjust b_haz
  # }
  #######

  ## Old plots
  # # plot KM curves for dependent and control datasets
  # plot(survfit(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == stime,]), mark.time = TRUE)
  # plot(survfit(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == stime,]), mark.time = TRUE)
  #
  #
  # # trtxdat <- xdat[xdat$time == 1, ]
  # # test <- cbind(trtxdat, sdat)
  # # plot(survfit(survival::Surv(eventtime, status) ~ arm, data = test), mark.time = TRUE)
  # survminer::ggsurvplot(
  #   fit = survfit(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,]),
  #    conf.int = TRUE)
  # survminer::ggsurvplot(
  #   fit = survfit(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,]),
  #   conf.int = TRUE)





  # Run control model (i.e., get do-treat, unconfounded estimates)
  unbiased_est <- function( data, indices){ # function to pass to boot::boot()
    d <- data[indices,] # allows boot to select sample
    fit <-   survival::coxph(survival::Surv(eventtime, status) ~ arm, data = d)
    return(exp(coef(fit)))
  }
  unbiased <- boot::boot(data = fulldat_cont[fulldat_cont$time == 1,], statistic = unbiased_est, R = bootrep)
  # Old, single rep model: survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,])

  unbiased_plot <- survminer::ggsurvplot(
    fit = survminer::surv_fit(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,]),
    xlab = "Time",
    ylab = "OS",
    title = "KM Plots for Unconfounded Analysis", conf.int = TRUE) %++%
    geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


  # Run Naive models:

  ## ITT ##
  # HR estimate:
  itt_est <- function( data, indices){ # function to pass to boot::boot()
    # select samples, run a survival model, and extract exp() of coefficient to return estimated hazard ratio
    return(exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = data[indices,]))))
  }
  itt <- boot::boot(data = fulldat[fulldat$time == 1,], statistic = itt_est, R = bootrep)
  # Old, single rep model: survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,])

  itt_plot <- survminer::ggsurvplot(
    fit = survminer::surv_fit(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,]),
    xlab = "Time",
    ylab = "OS",
    title = "KM Plots for ITT", conf.int = TRUE) %++%
    geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


  ## Per Protocol ##
  # Per Protocol (censor at switch with no correction)
  PPtime <- sapply(unique(fulldat$ids), function(x) stime - sum(fulldat$treat[fulldat$ids == x])) # get total time off treatment. since there is no back-switching, tihs is Per Protocol time for control patients
  PPtime <- rep(PPtime, each=stime) # expand PP time
  fulldat$PPtime <- ifelse(fulldat$arm == 1 | fulldat$eventtime <= PPtime, fulldat$eventtime, PPtime) # if patient is in experimental group or already sees event before PPtime, use eventtime. Else censor at switch
  fulldat$PPdeath <- ifelse(fulldat$PPtime < fulldat$eventtime, 0, fulldat$status) # within control group, if death had been recorded but is now censored, censor event indicator

  # HR estimate
  pp_est <- function(data, indices){ # function to pass to boot::boot()
    # select samples, run a survival model, and extract exp() of coefficient to return estimated hazard ratio
    return(exp(coef(survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = data[indices,] ))))
  }
  pp <- boot::boot(data = fulldat, statistic = pp_est, R = bootrep)
  # Old, single rep model: survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat)

  # KM curves:
  pp_plot <- survminer::ggsurvplot(
    fit = survminer::surv_fit(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat[fulldat$time == 1, ]),
    xlab = "Time",
    ylab = "OS",
    title = "KM Plots for PP", conf.int = TRUE) %++%
    geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)

  # Run adjustment models:
  if(verbose > 1){
    print("Performing complex method estimates: IPCW...")
  }

  ## IPCW ##

  fulldat$starttime <- fulldat$time - 1 # add starttime
  fulldat$switch_time <- rep(sapply(unique(fulldat$ids), function(x) ifelse( # get first
    fulldat$arm[fulldat$ids == x][1] == 0 & fulldat$switch_status[fulldat$ids == x][1] == 1, # if id is a control pt. and observes switch,
    fulldat$time[fulldat$ids == x & fulldat$treat == 1][1], # get the first treatment time
    0)), each = stime)
  fulldat$switch <- ifelse(fulldat$time == fulldat$switch_time, 1, 0)
  # cut fulldat at censoring or event time. Recast several fulldat variables
  fulldat_cut <- fulldat[fulldat$time <= fulldat$eventtime,] # cut data past eventtime
  cens_ids <- unique(fulldat_cut$ids[fulldat_cut$switch == 1]) # if you get censored, you are added to cens_ids
  fulldat_cut$status[fulldat_cut$ids %in% cens_ids] <- 0 # if you get censored, your status becomes 0
  fulldat_cut <- fulldat_cut[!(fulldat_cut$switch_status == 1 & fulldat_cut$switch_time < fulldat_cut$time),] # cut data from control arm switcher past the point of switch
  # TODO eventstatus is, for
  fulldat_cut$eventstatus <- ifelse(fulldat_cut$time == floor(fulldat_cut$eventtime) & fulldat_cut$status == 1, 1, 0) # new eventstatus is 1 only if event observed before switch.
  fulldat_cut$time <- as.numeric(fulldat_cut$time)
  fulldat_cut$arm <- as.factor(fulldat_cut$arm)

  ipform <- formula(paste("survival::Surv(starttime, time, eventstatus) ~ arm + cluster(ids) +", paste(c(bcov_names, tcov_names), collapse = "+"), collapse = " "))

  if(!ipcw_robust){
    t <- c()
    for(i in 1:10){ # TODO Cap ipcw bootrep
      bootids <- sample(unique(fulldat_cut$ids), n, replace = TRUE)
      # fit ipcw model
      if("All" %in% violate | "IPCW" %in% violate){
        ipdat <- ipcw(data = fulldat_cut[fulldat_cut$ids %in% bootids,], id = "ids", tstart = starttime,
                      tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names[1:(length(tcov_names) - hide_tvar)],
                      type = "kaplan-meier", trunc = 0.05)
      }else{
        ipdat <- ipcw(data = fulldat_cut[fulldat_cut$ids %in% bootids,], id = "ids", tstart = starttime,
                    tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names,
                    type = "kaplan-meier", trunc = 0.05)
      }
      # HR estimate
      msm <- survival::coxph(survival::Surv(starttime, time, eventstatus) ~ arm , data = ipdat, weights = ipdat$weights.trunc)
      t[i] <- exp(coef(msm))
    }
    # on full data
    if("All" %in% violate | "IPCW" %in% violate){
      ipdat <- ipcw(data = fulldat_cut, id = "ids", tstart = starttime,
                    tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names[1:(length(tcov_names) - hide_tvar)],
                    type = "kaplan-meier", trunc = 0.05)
    }else{
      ipdat <- ipcw(data = fulldat_cut, id = "ids", tstart = starttime,
                  tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names,
                  type = "kaplan-meier", trunc = 0.05)
    }
    # HR estimate
    msm <- survival::coxph(survival::Surv(starttime, time, eventstatus) ~ arm , data = ipdat, weights = ipdat$weights.trunc)
    t0 <- exp(coef(msm))
  } else{ # use robust SE, and skip bootstrapping
    # on full data
    if("All" %in% violate | "IPCW" %in% violate){
      ipdat <- ipcw(data = fulldat_cut, id = "ids", tstart = starttime,
                    tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names[1:(length(tcov_names) - hide_tvar)],
                    type = "kaplan-meier", trunc = 0.05)
    }else{
      ipdat <- suppressWarnings(ipcw(data = fulldat_cut, id = "ids", tstart = starttime,
                  tstop = time, cens = switch, arm = "arm", bas.cov = bcov_names, conf = tcov_names,
                  type = "kaplan-meier", trunc = 0.05))
    }
    # HR estimate
    msm <- survival::coxph(survival::Surv(starttime, time, eventstatus) ~ arm, data = ipdat, weights = ipdat$weights.trunc)
    robust_se <- sqrt(diag(msm$var))
    t <- exp(c(coef(msm), coef(msm) - robust_se, coef(msm) + robust_se))
    t0 <- exp(coef(msm))

  }

  # collect HR estimates
  msm_hr <- list(t = t, t0 = t0)

  # make weighted KM estimates
  # TODO This doesnt work... at all
  ipcw_plot <- survminer::ggsurvplot(
    fit = survminer::surv_fit(survival::Surv(starttime, time, eventstatus) ~ arm, data = ipdat, weights = ipdat$weights.trunc),
    xlab = "Time",
    ylab = "OS",
    title = "KM Plots for IPCW", conf.int = TRUE) %++%
    geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)

  if(verbose > 1){
    print("Performing complex method estimates: RPSFTM...")
  }
  ## RPSFTM ##

  # get proportion of treatment per patient
  # TODO should this function actually return treatment UP UNTIL eventtime, i.e., including continuous time until eventtime?
  rx <- sapply(unique(fulldat$ids), function(x) sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime])/
                 fulldat$eventtime[fulldat$ids==x & fulldat$time==1])
  rpsft_dat <- cbind(fulldat[fulldat$time == 1, ], rx) # create condensed RPSFTM dataset
  rpsft_dat$cens <- stime # add an administrative censoring time, which is simply stime
  rpsft_dat$arm <- as.factor(rpsft_dat$arm)

  # single rep mod: Build either recensored or unrecensored model
  if(recens == TRUE){
    # build rpsftm model
    mr <- suppressWarnings(rpsftm(survival::Surv(eventtime, status) ~ rand(arm, rx), data = rpsft_dat, low_psi = -2, hi_psi = 2, censor_time = cens))
  } else{
    mr <- suppressWarnings(rpsftm(survival::Surv(eventtime, status) ~ rand(arm, rx), data = rpsft_dat, low_psi = -2, hi_psi = 2))
  }

  # set counterfactuals with rpsftm model object:
  rpsft_dat$counterfact <- rpsft_dat$eventtime # set default counterfactual
  rpsft_dat$counterfact[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 1] # get rpsftm counterfactual times
  rpsft_dat$cf_status <- rpsft_dat$status
  rpsft_dat$cf_status[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 2] # get rpsftm counterfactual death/censoring

  # HR estimate
  rpsft_est <- function(data, indices){ # function to pass to boot::boot()
    d <- data[indices,] # allows boot to select sample
    # Build either recensored or unrecensored model
    if(recens == TRUE){
      # build rpsftm model
      mr <- rpsftm(survival::Surv(eventtime, status) ~ rand(arm, rx), data = d, low_psi = -3, hi_psi = 3, censor_time = cens)
    } else{
      mr <- rpsftm(survival::Surv(eventtime, status) ~ rand(arm, rx), data = d, low_psi = -3, hi_psi = 3)
    }

    # set counterfactuals with rpsftm model object:
    d$counterfact <- d$eventtime # set default counterfactual
    d$counterfact[d$arm == 0] <- mr$Sstar[d$arm == 0, 1] # get rpsftm counterfactual times
    d$cf_status <- d$status
    d$cf_status[d$arm == 0] <- mr$Sstar[d$arm == 0, 2] # get rpsftm counterfactual death/censoring
    fit <- survival::coxph(survival::Surv(counterfact, cf_status) ~ arm, data = d)
    return(exp(coef(fit)))
  }
  rpsft <- suppressWarnings(boot::boot(data = rpsft_dat, statistic = rpsft_est, R = bootrep )) # TODO reduce bootreps of rpsftm
  # Old, single rep model: survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat)

  # make KM estimates, censored
  rpsft_plot <- survminer::ggsurvplot(
    fit = survminer::surv_fit(survival::Surv(counterfact, cf_status) ~ arm, data = rpsft_dat),
    xlab = "Time",
    ylab = "OS",
    title = "KM Plots for RPSFTM", conf.int = TRUE) %++%
    geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)

  if(verbose > 1){
    print("Performing complex method estimates: TSE...")
  }
  ## TSE ##

  tsdat <- fulldat[fulldat$time == 1,] # Holder for wide format
  tscontrol <- fulldat[fulldat$secbase_observed == 1 & fulldat$arm == 0 & fulldat$Mtime == 1,] # take subset of pts who observe M and who are in arm == 0
  #tscontrol <- tscontrol[tscontrol$time < tscontrol$eventtime,] # subset again, removing pts who observe M (secondary baseline) after eventtime
  tscontrol$TSsurv <- tscontrol$eventtime - tscontrol$time
  tscontrol <- tscontrol[tscontrol$TSsurv > 0,] # exclude last-day-switchers

  TSEform <- formula(paste("survival::Surv(TSsurv, status) ~ switch_status +", paste(c(bcov_names, tcov_names), collapse = "+"), collapse = " "))


  # TODO are we bootstrapping the wrong dataset? should it not be tsdat that gets bootstrapped, not tscontrol?
  tse_est <- function(data, indices, tsdat){ # function to pass to boot::boot()
    d <- data[indices,] # allows boot to select sample

    # fit AF model
    mod <- survreg(formula = TSEform, dist = tse_dist, data = d)
    AF <- exp(coef(mod))[names(exp(coef(mod))) == "switch_status"] # get acceleration factor
    d$counterfact <- ifelse(d$switch_status == 0,
                            d$TSsurv + d$time, # observed survival
                            (d$TSsurv / AF) + d$time # counterfactual survival
    )

    #
    df_boot <- tsdat[!tsdat$ids %in% data$ids,] # take subset of tsdat which is NOT in tscontrol
    df_boot$counterfact <- df_boot$eventtime # reset counterfactuals in df_boot
    df_boot <- rbindlist(list(df_boot, d), fill = TRUE) # fold d back into tsdat
    # for(i in unique(d$ids)){
    #   print(paste("tsdat:", length(tsdat$counterfact[tsdat$ids == i])))
    #   print(paste("d: ", length(d$counterfact[d$ids == i])))
    #   tsdat$counterfact[tsdat$ids == i] <- d$counterfact[d$ids == i]
    # }
    fit <- survival::coxph(survival::Surv(counterfact, status) ~ arm, data = df_boot)
    return(exp(coef(fit)))
  }
  tse_wrapper <- possibly(tse_est, otherwise = NA)
  tse <- boot::boot(data = tscontrol, statistic = tse_wrapper, R = bootrep, tsdat = tsdat)

  #########
  # fit AF model
  mod <- survreg(formula = TSEform, dist = tse_dist, data = tscontrol)

  AF <- exp(coef(mod))[names(exp(coef(mod))) == "switch_status"] # get acceleration factor
  tscontrol$counterfact <- ifelse(tscontrol$switch_status == 0,
                                  tscontrol$TSsurv + tscontrol$time, # observed survival
                                  (tscontrol$TSsurv / AF) + tscontrol$time # counterfactual survival
  )

  tsdat$counterfact <- tsdat$eventtime # reset counterfactuals
  for(i in unique(tscontrol$ids)){
    tsdat$counterfact[tsdat$ids == i] <- tscontrol$counterfact[tscontrol$ids == i]
  }
  #########

  # make KM estimates, censored
  tse_plot <- survminer::ggsurvplot(
    fit = survminer::surv_fit(survival::Surv(counterfact, status) ~ arm, data = tsdat),
    xlab = "Time",
    ylab = "OS",
    title = "KM Plots for TSE", conf.int = TRUE) %++%
    geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


  # TODO Recensor!

  # build boxplot dataset
  df <- data.frame(Unbiased = rep(NA, bootrep), ITT = NA, PP = NA, IPCW = NA, RPSFTM = NA, TSE = NA)
  df$Unbiased[1:length(unbiased$t)] <- unbiased$t
  df$ITT[1:length(itt$t)] <- itt$t
  df$PP[1:length(pp$t)] <- pp$t
  df$IPCW[1:length(msm_hr$t)] <-  msm_hr$t
  df$RPSFTM[1:length(rpsft$t)] <- rpsft$t
  df$TSE[1:length(tse$t)] <- tse$t
  df <- df %>% pivot_longer(names(df), names_to = "Method", values_to = "est") # df wide to long, for ggplot
  compar_plot <- ggplot(df, aes(Method, est, color = Method)) + geom_boxplot() +
    scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.title.x=element_blank(),
                                                               axis.text.x=element_blank(),
                                                               axis.ticks.x=element_blank()) +
    ylab("HR Estimate") +
    ggtitle("HR Estimates Across Methods") +
    geom_hline(yintercept=unbiased$t0, linetype="dashed", size=0.1, alpha = 0.5)

  # TODO produce point-comparison estimate of bias.
  # bias_plot <- ggplot(df, aes())

  secondary_baseline_observed <- sum(fulldat$secbase_observed[fulldat$time == 1], na.rm = TRUE) / n
  switch_observed <- sum(fulldat$switch_status[fulldat$time == 1], na.rm = TRUE) / length(id_con)

  if(verbose > 1){
    print("Done!")
  }


  return(list(unbiased = unbiased, unbiased_plot = unbiased_plot, itt = itt,
              itt_plot = itt_plot, pp = pp, pp_plot = pp_plot,
              ipcw = msm_hr, ipcw_plot = ipcw_plot, rpsft = rpsft,
              rpsft_plot = rpsft_plot, tse = tse, tse_plot = tse_plot, compar_plot = compar_plot,
              params = list(bootrep = bootrep, violate = violate, cens_flag = cens_flag, prop_cens = prop_cens,
                            recens = recens, stime = stime, num_bvar = num_bvar, num_tvar = num_tvar,
                            hide_tvar = hide_tvar, add_tvar = add_tvar, n = n, rerun_lim = rerun_lim,
                            prop_trt = prop_trt, prop_switch = prop_switch, prop_cont_event = prop_cont_event, switch_coef = switch_coef, bcov = bcov,
                            covar_coef = covar_coef, s_haz = s_haz, s_allowance = s_allowance,
                            m_hard = m_hard, m_haz = m_haz, m_inflation = m_inflation, m_fidelity = m_fidelity, beta.mat = beta.mat, haz = haz, b_haz = b_haz,
                            secondary_baseline_observed = secondary_baseline_observed,
                            switch_observed = switch_observed,
                            iterations = switch_iter, proportion_events_control = current_b_prop,
                            proportion_events_experimental = current_t_prop,
                            proportion_switchers = current_s_prop,
                            confounded_data = fulldat,
                            unconfounded_data = fulldat_cont,
                            compar_df = df)))

}



has_no_multidim_mat <- function(l){
  if(!is.list(l)) stop("Something went wrong; l isn't a list")
  # get elements which are matrices
  matrices <- l[which(lapply(l,is.matrix) == TRUE)]
    if(length(matrices) < 1) return(TRUE) # cut short if there are no matrices
  # get boolean for whether matrix is multidimensional
  if(any(lapply(matrices, function(x) dim(x)[2]) > 1)) return(FALSE)
  else(return(TRUE))
}

has_more_than_one_multidim_mat <- function(l){
  if(!is.list(l)) stop("Something went wrong; l isn't a list")
  # get elements which are matrices
  matrices <- l[which(lapply(l,is.matrix) == TRUE)]
  if(length(matrices) < 1) return(FALSE) # cut short if there are no matrices
  # get booleans for whether matrix is multidimensional
  if(sum(lapply(matrices, function(x) dim(x)[2]) > 1) > 1) return(TRUE)
  else(return(FALSE))
}

has_one_multidim_mat <- function(l){
  if(!is.list(l)) stop("Something went wrong; l isn't a list")
  # get elements which are matrices
  matrices <- l[which(lapply(l,is.matrix) == TRUE)]
  if(length(matrices) < 1) return(FALSE) # cut short if there are no matrices
  # get booleans for whether matrix is multidimensional
  if(sum(lapply(matrices, function(x) dim(x)[2]) > 1) == 1) return(TRUE)
  else(return(FALSE))
}

which_multidim_mat <- function(l){
  if(!is.list(l)) stop("Something went wrong; l isn't a list")
  # get elements which are matrices
  matrices <- l[which(lapply(l,is.matrix) == TRUE)]
  # get booleans for whether matrix is multidimensional
  return(which(lapply(matrices, function(x) dim(x)[2]) > 1))
}


#' Title
#'
#' @param add_tvar
#' @param b_allowance
#' @param b_haz
#' @param b_mag
#' @param b_scale
#' @param b_shape
#' @param bcov
#' @param beta.mat
#' @param bootrep
#' @param cens_flag
#' @param covar_coef
#' @param dep_func
#' @param haz
#' @param hide_tvar
#' @param ipcw_robust
#' @param m_allowance
#' @param m_inflation
#' @param m_fidelity
#' @param m_hard
#' @param m_haz
#' @param m_mag
#' @param m_scale
#' @param m_shape
#' @param n
#' @param num_bvar
#' @param num_tvar
#' @param prop_cens
#' @param prop_cens_allowance
#' @param prop_cont_event
#' @param prop_switch
#' @param prop_trt
#' @param prop_trt_event
#' @param recens
#' @param rerun_lim
#' @param s_allowance
#' @param s_haz
#' @param s_mag
#' @param s_scale
#' @param s_shape
#' @param stime
#' @param switch_coef
#' @param t_allowance
#' @param t_mag
#' @param treat_beta
#' @param tse_dist
#' @param unfix
#' @param verbose
#' @param violate
#'
#' @return
#' @export
#'
#' @examples
simrange <- function(add_tvar, b_allowance, b_haz, b_mag, b_scale, b_shape, bcov, beta.mat, bootrep,
                     cens_flag, covar_coef, dep_func, haz, hide_tvar, ipcw_robust,
                     m_allowance, m_inflation, m_fidelity, m_hard, m_haz, m_mag, m_scale, m_shape,
                     n, num_bvar, num_tvar, prop_cens, prop_cens_allowance, prop_cont_event, prop_switch,
                     prop_trt, prop_trt_event, recens, rerun_lim, s_allowance, s_haz, s_mag, s_scale,
                     s_shape, stime, switch_coef, t_allowance, t_mag, treat_beta, tse_dist, unfix, verbose, violate){

  # variables that should accept a range value: add_tvar, b_haz, b_scale, b_shape, hide_tvar, m_inflation, m_fidelity, m_haz, m_scale, m_shape,
  # n, num_bvar, num_tvar, prop_cens, prop_cont_event, prop_switch,
  # prop_trt, prop_trt_event, s_haz, s_scale,
  # s_shape, stime, treat_beta

  ranges <- list(add_tvar, b_scale, b_shape, hide_tvar, m_inflation, m_fidelity, m_scale, m_shape,
                 n, num_bvar, num_tvar, prop_cens, prop_cont_event, prop_switch,
                 prop_trt, prop_trt_event, s_scale,
                 s_shape, stime, treat_beta)
  haz_ranges <- list(b_haz, m_haz, s_haz)

  # Ensure that one and only one range is specified
  # Ensure that at least one is selected
  if(all(lengths(ranges) <= 1) & has_no_multidim_mat(haz_ranges) ){
    stop("You must select at least one range variable. For help, see ?simrange")
  }
  # Ensure that no more than 1 is selected
  if(sum(lengths(ranges) > 1) > 1 | has_more_than_one_multidim_mat(haz_ranges) |
     (any(lengths(ranges) > 1) & has_one_multidim_mat(haz_ranges)) ){
    stop("You can only select one range variable. For help, see ?simrange")
  }

  # set marker for whether range variable is in ranges, or in haz_ranges
  if(has_one_multidim_mat(haz_ranges)) multidim <- TRUE

  # store range_var index
  if(multidim){
    haz_range_index <- which_multidim_mat(haz_ranges)
  }else{
    range_index <- which(lengths(ranges) > 1)
  }

  # store range variable in range_var
  if(multidim){
    range_var <- haz_ranges[[haz_range_index]] # TODO get single range variable
  } else{
    range_var <- ranges[[range_index]] # TODO get single range variable
  }

  # original variable becomes first index of range_var
  if(multidim){
    haz_ranges[[haz_range_index]] <- range_var[,1] # TODO get single range variable
  }else{
    ranges[[range_index]] <- range_var[1] # TODO get single range variable
  }


  # TODO can we somehow use the ... argument here to pass outer arguments to simswitch() call?
  first_sim <- simswitch(
    add_tvar = ,
    b_allowance = 0.1,
    bcov = ,
    beta.mat = ,
    b_haz = ,
    bootrep = ,
    covar_coef = ,
    # dep_func = sequence[[i-1]]$params$dep_func,
    # haz = sequence[[i-1]]$params$haz(),
    hide_tvar = ,
    # ipcw_robust = sequence[[i-1]]$params$ipcw_robust,
    # m_allowance = sequence[[i-1]]$params$m_allowance,
    m_fidelity = ,
    m_hard = ,
    m_haz = ,
    m_inflation = ,
    n = ,
    num_bvar = ,
    num_tvar = ,
    prop_cens = ,
    # prop_cens_allowance = sequence[[i-1]]$params$prop_cens_allowance,
    prop_cont_event = ,
    prop_switch = ,
    prop_trt = ,
    prop_trt_event = ,
    recens = ,
    rerun_lim = ,
    s_allowance = 0.025,
    s_haz = ,
    stime = ,
    switch_coef = ,
    t_allowance = 0.1,
    # tse_dist = sequence[[i-1]]$params$tse_dist,
    unfix = c("B", "M", "T")
  )

  # replace default values with previous simulation values
  ranges[["n"]] = first_sim$params$n
  ranges[["num_tvar"]] = first_sim$params$num_tvar
  ranges[["num_bvar"]] = first_sim$params$num_bvar
  ranges[["prop_trt"]] = first_sim$params$prop_trt
  ranges[["prop_switch"]] = first_sim$params$prop_switch
  ranges[["s_allowance"]] = 0.025
  ranges[["prop_cens"]] = first_sim$params$prop_cens
  ranges[["prop_cont_event"]] = first_sim$params$prop_cont_event
  ranges[["b_allowance"]] = 0.1
  ranges[["prop_trt_event"]] = first_sim$params$prop_trt
  ranges[["t_allowance"]] = 0.1
  ranges[["bootrep"]] = first_sim$params$bootrep
  ranges[["recens"]] = first_sim$params$recens
  ranges[["b_haz"]] = first_sim$params$b_haz
  ranges[["m_inflation"]] = first_sim$params$m_inflation
  ranges[["m_fidelity"]] = first_sim$params$m_fidelity
  ranges[["m_hard"]] = FALSE
  ranges[["unfix"]] = c("B", "M", "T")
  ranges[["bcov"]] = first_sim$params$bcov
  ranges[["beta.mat"]] = first_sim$params$beta.mat
  ranges[["covar_coef"]] = first_sim$params$covar_coef
  # ranges[["dep_func"]] = first_sim$params$dep_func
  # ranges[["haz"]] = first_sim$params$haz()
  # ranges[["ipcw_robust"]] = first_sim$params$ipcw_robust
  ranges[["m_haz"]] = first_sim$params$m_haz
  ranges[["s_haz"]] = first_sim$params$s_haz
  # ranges[["tse_dist"]] = first_sim$params$tse_dist
  ranges[["add_tvar"]] = first_sim$params$add_tvar
  ranges[["hide_tvar"]] = first_sim$params$hide_tvar
  # ranges[["m_allowance"]] = first_sim$params$m_allowance
  # ranges[["prop_cens_allowance"]] = first_sim$params$prop_cens_allowance
  ranges[["rerun_lim"]] = first_sim$params$rerun_lim
  ranges[["stime"]] = first_sim$params$stime
  ranges[["switch_coef"]] = first_sim$params$switch_coef

  # overwrite range variable
  ranges[[range_index]] <- range_var[2] # replace range variable with single value


  # fill out range sequence with sequential simulations
  if(multidim){
    # TODO fill in multivariate version

  }else{
    range_seq <- list() # set new list for simulation sequence
    range_seq[[1]] <- first_sim # set first simulation

    for(i in 2:length(range_var)){ # iterate through range_var values
      ranges[[range_index]] <- range_var[i] # replace range variable with single value
      range_seq[[i]] <- simswitch(
        n = range_seq[[i-1]]$params$n,
        num_tvar = range_seq[[i-1]]$params$num_tvar,
        num_bvar = range_seq[[i-1]]$params$num_bvar,
        prop_trt = range_seq[[i-1]]$params$prop_trt,
        prop_switch = switch_i[[i]],
        s_allowance = 0.025,
        prop_cens = range_seq[[i-1]]$params$prop_cens,
        prop_cont_event = range_seq[[i-1]]$params$prop_cont_event,
        b_allowance = 0.1,
        prop_trt_event = range_seq[[i-1]]$params$prop_trt,
        t_allowance = 0.1,
        bootrep = range_seq[[i-1]]$params$bootrep,
        recens = range_seq[[i-1]]$params$recens,
        b_haz = range_seq[[i-1]]$params$b_haz,
        m_inflation = range_seq[[i-1]]$params$m_inflation,
        m_fidelity = range_seq[[i-1]]$params$m_fidelity,
        m_hard = FALSE,
        unfix = c("B", "M", "T"),
        bcov = range_seq[[i-1]]$params$bcov,
        beta.mat = range_seq[[i-1]]$params$beta.mat,
        covar_coef = range_seq[[i-1]]$params$covar_coef,
        # dep_func = range_seq[[i-1]]$params$dep_func,
        # haz = range_seq[[i-1]]$params$haz(),
        # ipcw_robust = range_seq[[i-1]]$params$ipcw_robust,
        m_haz = range_seq[[i-1]]$params$m_haz,
        s_haz = range_seq[[i-1]]$params$s_haz,
        # tse_dist = range_seq[[i-1]]$params$tse_dist,
        add_tvar = range_seq[[i-1]]$params$add_tvar,
        hide_tvar = range_seq[[i-1]]$params$hide_tvar,
        # m_allowance = range_seq[[i-1]]$params$m_allowance,
        # prop_cens_allowance = range_seq[[i-1]]$params$prop_cens_allowance,
        rerun_lim = range_seq[[i-1]]$params$rerun_lim,
        stime = range_seq[[i-1]]$params$stime,
        switch_coef = range_seq[[i-1]]$params$switch_coef
      )
    }
    }


}


# TODO testrun the main function below

test <- simswitch(
                n = 400,
                stime = 100,
                num_tvar = 3,
                prop_trt = 0.5,
                prop_switch = 0.75,
                s_shape = 4,
                s_scale = 120,
                s_allowance = 0.05,
                prop_cens = 0,
                prop_cont_event = 0.7,
                b_allowance = 0.05,
                bootrep = 100,
                recens = FALSE,
                m_inflation = 2,
                m_fidelity = 0.15,
                m_hard = TRUE,
                m_shape = 3,
                m_scale = 140,
                unfix = c("T"),
                treat_beta = -0.43)
