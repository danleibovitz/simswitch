# main function

######
# TODO
# re-format styling
# move features out that could be written as own functions. these can be internal, not exported.
# make K-M curve for each dataset, and make K-M curves averages
# occasional error: Error in rpsft_dat$counterfact[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm ==  :
# replacement has length zero
######
#
# # load libraries
# library(MASS)
# # library(coxed)
# library(data.table)
# library(survival)
# library(ggplot2)
# library(survminer)
# library(rpsftm)
# # library(LICORS)
# library(ipcwswitch)
# library(boot)
# library(tidyr)
# library(purrr)

# recommended:
# library(parallel)

# set relevant classes

setOldClass("ggsurvplot")
setOldClass("gg")

setClassUnion("ggsurvplotOrNULL", c("ggsurvplot", "NULL"))

setClass("SwitchSimulation", representation(
  unbiased = "numeric",
  unbiased_plot = "ggsurvplot",
  itt = "numeric",
  itt_plot = "ggsurvplot",
  pp = "numeric",
  pp_plot = "ggsurvplot",
  ipcw = "numeric",
  ipcw_plot = "ggsurvplotOrNULL",
  rpsft = "numeric",
  rpsft_plot = "ggsurvplot",
  tse = "numeric",
  tse_plot = "ggsurvplot",
  compar_plot = "gg",
  bias_plot = "gg",
  params = "list"
  ))

#' Simulate TS data and apply treatment-effect-estimating methods
#'
#' @param add_tvar Select number of covariates to add which have no effect on overall survival, but which are included in IPCW, TSE and RPSFTM models.
#' @param b_allowance
#' @param b_haz Complete hazard vector
#' @param b_mag
#' @param b_scale
#' @param b_shape
#' @param bcov
#' @param beta.mat
#' @param cens_flag
#' @param covar_coef This is a new comment
#' @param dep_func
#' @param haz
#' @param hide_tvar
#' @param ipcw_exclude
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
simswitchloop <- function(add_tvar, b_allowance, b_haz, b_mag, b_scale, b_shape, bcov, beta.mat,
                      cens_flag, covar_coef, dep = dep_func, haz = haz_func, hide_tvar, ipcw_exclude,
                      m_allowance, m_inflation, m_fidelity, m_hard, m_haz, m_mag, m_scale, m_shape,
                      n, num_bvar, num_cores, num_tvar, prop_cens, prop_cens_allowance, prop_cont_event, prop_switch,
                      prop_trt, prop_trt_event, recens, reps, rerun_lim, s_allowance, s_haz, s_mag, s_scale,
                      s_shape, seed, stime, switch_coef, t_allowance, t_mag, treat_beta, tse_dist, unfix, verbose, violate){



  if(missing(num_cores)){ # set parmeters for parallelization
    para <- FALSE
  } else{
    # TODO How can we check here that the package 'parallel' is available? This is a requirement.
    para <- TRUE
    if(!(num_cores %% 1 == 0) | num_cores < 1){ # check that num_cores is a positive integer
      stop("The num_cores argument must be a positive integer.")
    }
    if(num_cores > detectCores()){ # make sure number of cores requested are available.
      num_cores <- detectCores()
      warning("Requested number of cores for parallelization exceeds available number Resetting to available number.")
    }
  }

  if(missing(verbose)){
    verbose <- 2
  }
  if(verbose > 1){
    print("Setting parameters...")
  }
  # set ipcw_exclude. If TRUE, we skip bootstrapping for IPCW. Bootstrapped IPCW may not be symmetrical, and may be more accurate, but takes hella long
  if(missing(ipcw_exclude)){
    ipcw_exclude <- TRUE
  }

  if(missing(tse_dist)){ # alternatives are weibull, lognormal, etc.
    tse_dist <- "loglogistic"
  }

  # set number of datasets to be repeated
  if(missing(reps)){
    reps <- 1000
  }
  # set assumption violation flag
  #  implement automatic assumption violation. For RPSFTM, non-constant treatment effect. For TSE, switching after secondary baseline. For IPCW, unmeasured confounding
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
    hide_tvar <- num_tvar - 1 # exclude all tvar except M
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
  if(missing(switch_coef)){ # baseline covariate attached to "predisposition" baseline is first, and slightly higher than other baselines.
    switch_coef <- c(runif(1, 1.5, 2), runif(num_bvar-1, 0.1, 0.3), runif(num_tvar, 0.2, 0.5)) # default of switch_coef log hazard ratios. baseline switch coefs are smaller.
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

  # set random covariates, number equal to num_bvar
  if(missing(bcov)){ # if a baseline covariate matrix is not specified, generate a random one
    bcov <- gen_bcov(num_bvar = num_bvar, diags = 0.4, middle = 1, stime = stime)
  } else{
    if(dim(bcov)[1] != (n*stime) | dim(bcov)[2] != num_bvar) stop("a pre-specified baseline covariate matrix must have length equal to (number of patients)*(stime) and width equal to num_bvar")
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
  fulldat_cont <- fulldat # set a control group, where there is no time-dependent confounding

  # set a default covar_coef list. If user defines their own dep_func, it must call covar_coef argument even if it doesnt use it
  if(missing(covar_coef)){
    covar_coef <- list(baseline = matrix(sample(1:(num_bvar*num_tvar), num_bvar*num_tvar), ncol = num_tvar),
                       varying = matrix(sample(1:(num_tvar*(num_tvar+1)), num_tvar*(num_tvar+1)), ncol = num_tvar)) # add a coefficient row for the treatment effect
    # covar_coef$baseline[1,] <- covar_coef$baseline[1,] + num_bvar*num_tvar*0.25 # beef up effect of baseline predisposition on first covars
    covar_coef$baseline <- LICORS::normalize(covar_coef$baseline, byrow = FALSE)
    # sweep(covar_coef$baseline, 2, colSums(covar_coef$baseline), `/`)
    covar_coef$varying <- LICORS::normalize(covar_coef$varying, byrow = FALSE)
    # sweep(covar_coef$varying, 2, colSums(covar_coef$varying), `/`)
    covar_coef$varying[1,] <- -covar_coef$varying[1,] # make the treatment effect protective
    #  try to make more orderly...
    covar_coef$varying <- covar_coef$varying/100
    covar_coef$varying[2:(num_tvar+1),1] <- covar_coef$varying[2:(num_tvar+1),1]*150 # bump up coefficients predictive of secondary baseline
    for(j in 1:num_tvar){
      covar_coef$varying[(j+1),j] <- 1
    }
  }

  # set dep_func (dependent function for generating covariates). must always return matrix of length n, width num_tvar + 1. at least first index of num_tvar is binary
  # if(missing(dep_func)){
  #   dep_func <- dep_func # if user does not provide a dependency function, use the default
  # }

  # set baseline switch hazard function. lam is set as lambda of exp distribution with expected value of prop_switch by stime
  if(missing(s_haz)){
    s_haz <- weihaz(1:stime, 2, 0.7*stime)
  }
  if(!missing(s_shape) & !missing(s_scale)){
    s_haz <- weihaz(1:stime, s_shape, s_scale)
  }

  ## this section is a bit tedious. We have to set up parameters to iteratively search for the correct switching proportion
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
  if("All" %in% violate | "TSE" %in% violate){
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

  if(verbose > 1){
    print("Setting baseline hazard function...")
  }
  # if(missing(haz)){
  #   haz <- haz # if user does not define a hazard function, use the default
  # }

  if(verbose > 1){
    print("Filling data frames")
  }

  # set empty estimate vectors
  unbiased <- c()
  itt <- c()
  pp <- c()
  msm_hr <- c()
  rpsft <- c()
  tse <- c()
  secondary_baseline_observed <- c()
  switch_observed <- c()
  switch_iter <- 0 # how many times have we searched on the inner loop?
  Minb0 <- c()
  Minb1 <- c()


  # giant while loop should begin here. At this point, we have all parameters set, and a fulldat dataframe with no time
  # varying covariates, no switching and no secondary baseline
  # TODO if para == TRUE, parallelize here. Can these be done all in parallel (including the 1st rep, where conditions are particular) or is it sequential?
  for(r in 1){

    if(r != 1) unfix <- c("B", "M", "S", "T") # if we have generated the first dataset, unfix all hazards

    while(rerun){ # iteratively update switching hazard function until we get the right proportion. The first term in this check is the number of pts whose final treatment indicator

      for(i in 1:stime){
        # set covars
        fulldat[fulldat$time == i, names(fulldat) %in% c(tcov_names, "Mtime")] <-
          dep(dat = fulldat, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz)
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
          dep(dat = fulldat_cont, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz)
      }


      if(verbose > 1){
        print("Generating confounded survival times...")
      }

      sdat <- survive( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
                       b_haz = b_haz)

      if(verbose > 1){
        print("Generating un-confounded survival times...")
      }

      sdat_cont <- survive( x = fulldat_cont[fulldat_cont$arm == 0,], hazard = haz, betas = beta.mat[fulldat_cont$arm == 0,], ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat_cont$ids[fulldat_cont$arm==0]),
                            b_haz = b_haz)


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
        ifelse(sum(fulldat$M[fulldat$ids == x]) > 0 && fulldat$time[fulldat$Mtime == 1 & fulldat$ids == x] <= fulldat$eventtime[fulldat$Mtime == 1 & fulldat$ids == x], # id has at least some m, and time at first m is less than event time
               1,
               0)), each = stime))

      # Does the patient (only control) observe switch _before eventtime_?
      fulldat$switch_status <- rep(sapply(unique(fulldat$ids), function(x)
        ifelse(sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime]) > 0 & fulldat$arm[fulldat$ids == x][1] == 0,
               1,
               0)), each = stime)

      if(r == 1){ # if we're in the first loop,
        switch_iter <- switch_iter + 1 # iterate search indicator
      }
      if(switch_iter > rerun_lim) stop("Your covariate model failed to converge. Try different hazards and/or coefficients")

      rerun <- FALSE # pre-empt a rerun, unless the following conditions are met

      # if not enough M occurences, adjust M hazard
      # TODO for now, M occurence is only being adjusted upward. do we want to give it a within-window type adjustment?
      current_m_prop <- length(fulldat$ids[fulldat$arm == 0 & fulldat$Mtime == 1 & fulldat$time <= fulldat$eventtime])
      if(! "M" %in% unfix){
        if(current_m_prop < min(n*prop_trt, prop_switch*prop_trt*n*m_inflation)){ # current_m_prop is proportion of control arm that observes M. if it's less than the min of (the whole control arm, the switchers*m_inflation ), increase m_haz and rerun
          m_haz <- m_haz*m_mag # change scale of m related haz, and rerun
          rerun <- TRUE
        }
      }


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

    rerun <- TRUE # reset inner loop rerun flag

    # if(verbose > 1){
    #   print("Performing naive method estimates...")
    # }


    # Run control model (i.e., get do-treat, unconfounded estimates)
    unbiased[r] <- exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,])))
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
    itt[r] <- exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,])))
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
    pp[r] <- exp(coef(survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat[fulldat$time == 1,] )))
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
    if(!ipcw_exclude){
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
      msm_hr[r] <- exp(coef(survival::coxph(survival::Surv(starttime, time, eventstatus) ~ arm , data = ipdat, weights = ipdat$weights.trunc)))
    }else{
      msm_hr[r] <- as.numeric(NA)
    }

    # make weighted KM estimates
    if(!ipcw_exclude){
      ipcw_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(starttime, time, eventstatus) ~ arm, data = ipdat, weights = ipdat$weights.trunc),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for IPCW", conf.int = TRUE) %++%
      geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)
    } else{
        ipcw_plot <- NULL
      }

    if(verbose > 1){
      print("Performing complex method estimates: RPSFTM...")
    }
    ## RPSFTM ##

    # get proportion of treatment per patient
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
    # TODO whole simuulator SOMETIMES bugs here:
    rpsft_dat$counterfact[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 1] # get rpsftm counterfactual times
    rpsft_dat$cf_status <- rpsft_dat$status
    rpsft_dat$cf_status[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 2] # get rpsftm counterfactual death/censoring
    rpsft[r] <- exp(coef(survival::coxph(survival::Surv(counterfact, cf_status) ~ arm, data = rpsft_dat)))

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

    # this or that ####
    # tse_est <- function(data, indices, tsdat){ # function to pass to boot::boot()
    #   d <- data[indices,] # allows boot to select sample
    #
    #   # fit AF model
    #   mod <- survreg(formula = TSEform, dist = tse_dist, data = d)
    #   AF <- exp(coef(mod))[names(exp(coef(mod))) == "switch_status"] # get acceleration factor
    #   d$counterfact <- ifelse(d$switch_status == 0,
    #                           d$TSsurv + d$time, # observed survival
    #                           (d$TSsurv / AF) + d$time # counterfactual survival
    #   )
    #
    #   #
    #   df_boot <- tsdat[!tsdat$ids %in% data$ids,] # take subset of tsdat which is NOT in tscontrol
    #   df_boot$counterfact <- df_boot$eventtime # reset counterfactuals in df_boot
    #   df_boot <- rbindlist(list(df_boot, d), fill = TRUE) # fold d back into tsdat
    #   # for(i in unique(d$ids)){
    #   #   print(paste("tsdat:", length(tsdat$counterfact[tsdat$ids == i])))
    #   #   print(paste("d: ", length(d$counterfact[d$ids == i])))
    #   #   tsdat$counterfact[tsdat$ids == i] <- d$counterfact[d$ids == i]
    #   # }
    #   fit <- survival::coxph(survival::Surv(counterfact, status) ~ arm, data = df_boot)
    #   return(exp(coef(fit)))
    # }
    # tse_wrapper <- possibly(tse_est, otherwise = NA)
    # tse <- boot::boot(data = tscontrol, statistic = tse_wrapper, R = bootrep, tsdat = tsdat)
    #########
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

    tse[r] <- exp(coef(survival::coxph(survival::Surv(counterfact, status) ~ arm, data = tsdat)))

    #########

    # make KM estimates, censored
    tse_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(counterfact, status) ~ arm, data = tsdat),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for TSE", conf.int = TRUE) %++%
      geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


    # TODO Recensor!



    secondary_baseline_observed[r] <- sum(fulldat$secbase_observed[fulldat$time == 1], na.rm = TRUE) / n
    switch_observed[r] <- sum(fulldat$switch_status[fulldat$time == 1], na.rm = TRUE) / length(id_con)
    Minb0[r] <- sum(fulldat$secbase_observed[fulldat$time == 1 & fulldat$arm == 0 & fulldat$b1==0])
    Minb1[r] <- sum(fulldat$secbase_observed[fulldat$time == 1 & fulldat$arm == 0 & fulldat$b1 == 1])

    if(verbose > 1){
      print(paste("repetition: ", r))
    }

  } # end giant loop. Build plots below.













  # giant while loop should begin here. At this point, we have all parameters set, and a fulldat dataframe with no time
  # varying covariates, no switching and no secondary baseline
  # TODO if para == TRUE, parallelize here. Can these be done all in parallel (including the 1st rep, where conditions are particular) or is it sequential?
  mclapply(X = 2:reps, mc.cores = 4, FUN = function(r){

    if(r != 1) unfix <- c("B", "M", "S", "T") # if we have generated the first dataset, unfix all hazards

    while(rerun){ # iteratively update switching hazard function until we get the right proportion. The first term in this check is the number of pts whose final treatment indicator

      for(i in 1:stime){
        # set covars
        fulldat[fulldat$time == i, names(fulldat) %in% c(tcov_names, "Mtime")] <-
          dep(dat = fulldat, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz)
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
          dep(dat = fulldat_cont, window = i, base_var=bcov_names, time_var=tcov_names, covar_coef = covar_coef, m_haz = m_haz)
      }


      if(verbose > 1){
        print("Generating confounded survival times...")
      }

      sdat <- survive( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
                       b_haz = b_haz)

      if(verbose > 1){
        print("Generating un-confounded survival times...")
      }

      sdat_cont <- survive( x = fulldat_cont[fulldat_cont$arm == 0,], hazard = haz, betas = beta.mat[fulldat_cont$arm == 0,], ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat_cont$ids[fulldat_cont$arm==0]),
                            b_haz = b_haz)


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
        ifelse(sum(fulldat$M[fulldat$ids == x]) > 0 && fulldat$time[fulldat$Mtime == 1 & fulldat$ids == x] <= fulldat$eventtime[fulldat$Mtime == 1 & fulldat$ids == x], # id has at least some m, and time at first m is less than event time
               1,
               0)), each = stime))

      # Does the patient (only control) observe switch _before eventtime_?
      fulldat$switch_status <- rep(sapply(unique(fulldat$ids), function(x)
        ifelse(sum(fulldat$treat[fulldat$ids == x & fulldat$time <= fulldat$eventtime]) > 0 & fulldat$arm[fulldat$ids == x][1] == 0,
               1,
               0)), each = stime)

      if(r == 1){ # if we're in the first loop,
        switch_iter <- switch_iter + 1 # iterate search indicator
      }
      if(switch_iter > rerun_lim) stop("Your covariate model failed to converge. Try different hazards and/or coefficients")

      rerun <- FALSE # pre-empt a rerun, unless the following conditions are met

      # if not enough M occurences, adjust M hazard
      # TODO for now, M occurence is only being adjusted upward. do we want to give it a within-window type adjustment?
      current_m_prop <- length(fulldat$ids[fulldat$arm == 0 & fulldat$Mtime == 1 & fulldat$time <= fulldat$eventtime])
      if(! "M" %in% unfix){
        if(current_m_prop < min(n*prop_trt, prop_switch*prop_trt*n*m_inflation)){ # current_m_prop is proportion of control arm that observes M. if it's less than the min of (the whole control arm, the switchers*m_inflation ), increase m_haz and rerun
          m_haz <- m_haz*m_mag # change scale of m related haz, and rerun
          rerun <- TRUE
        }
      }


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

    rerun <- TRUE # reset inner loop rerun flag

    # if(verbose > 1){
    #   print("Performing naive method estimates...")
    # }


    # Run control model (i.e., get do-treat, unconfounded estimates)
    unbiased[r] <- exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat_cont[fulldat_cont$time == 1,])))
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
    itt[r] <- exp(coef(survival::coxph(survival::Surv(eventtime, status) ~ arm, data = fulldat[fulldat$time == 1,])))
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
    pp[r] <- exp(coef(survival::coxph(survival::Surv(PPtime, PPdeath) ~ arm, data = fulldat[fulldat$time == 1,] )))
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
    if(!ipcw_exclude){
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
      msm_hr[r] <- exp(coef(survival::coxph(survival::Surv(starttime, time, eventstatus) ~ arm , data = ipdat, weights = ipdat$weights.trunc)))
    }else{
      msm_hr[r] <- as.numeric(NA)
    }

    # make weighted KM estimates
    if(!ipcw_exclude){
      ipcw_plot <- survminer::ggsurvplot(
        fit = survminer::surv_fit(survival::Surv(starttime, time, eventstatus) ~ arm, data = ipdat, weights = ipdat$weights.trunc),
        xlab = "Time",
        ylab = "OS",
        title = "KM Plots for IPCW", conf.int = TRUE) %++%
        geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)
    } else{
      ipcw_plot <- NULL
    }

    if(verbose > 1){
      print("Performing complex method estimates: RPSFTM...")
    }
    ## RPSFTM ##

    # get proportion of treatment per patient
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
    # TODO whole simuulator SOMETIMES bugs here:
    rpsft_dat$counterfact[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 1] # get rpsftm counterfactual times
    rpsft_dat$cf_status <- rpsft_dat$status
    rpsft_dat$cf_status[rpsft_dat$arm == 0] <- mr$Sstar[rpsft_dat$arm == 0, 2] # get rpsftm counterfactual death/censoring
    rpsft[r] <- exp(coef(survival::coxph(survival::Surv(counterfact, cf_status) ~ arm, data = rpsft_dat)))

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

    #########
    # tse_est <- function(data, indices, tsdat){ # function to pass to boot::boot()
    #   d <- data[indices,] # allows boot to select sample
    #
    #   # fit AF model
    #   mod <- survreg(formula = TSEform, dist = tse_dist, data = d)
    #   AF <- exp(coef(mod))[names(exp(coef(mod))) == "switch_status"] # get acceleration factor
    #   d$counterfact <- ifelse(d$switch_status == 0,
    #                           d$TSsurv + d$time, # observed survival
    #                           (d$TSsurv / AF) + d$time # counterfactual survival
    #   )
    #
    #   #
    #   df_boot <- tsdat[!tsdat$ids %in% data$ids,] # take subset of tsdat which is NOT in tscontrol
    #   df_boot$counterfact <- df_boot$eventtime # reset counterfactuals in df_boot
    #   df_boot <- rbindlist(list(df_boot, d), fill = TRUE) # fold d back into tsdat
    #   # for(i in unique(d$ids)){
    #   #   print(paste("tsdat:", length(tsdat$counterfact[tsdat$ids == i])))
    #   #   print(paste("d: ", length(d$counterfact[d$ids == i])))
    #   #   tsdat$counterfact[tsdat$ids == i] <- d$counterfact[d$ids == i]
    #   # }
    #   fit <- survival::coxph(survival::Surv(counterfact, status) ~ arm, data = df_boot)
    #   return(exp(coef(fit)))
    # }
    # tse_wrapper <- possibly(tse_est, otherwise = NA)
    # tse <- boot::boot(data = tscontrol, statistic = tse_wrapper, R = bootrep, tsdat = tsdat)
    #########
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

    tse[r] <- exp(coef(survival::coxph(survival::Surv(counterfact, status) ~ arm, data = tsdat)))

    #########

    # make KM estimates, censored
    tse_plot <- survminer::ggsurvplot(
      fit = survminer::surv_fit(survival::Surv(counterfact, status) ~ arm, data = tsdat),
      xlab = "Time",
      ylab = "OS",
      title = "KM Plots for TSE", conf.int = TRUE) %++%
      geom_hline(yintercept=0.5, linetype="dashed", size=0.1, alpha = 0.5)


    # TODO Recensor!



    secondary_baseline_observed[r] <- sum(fulldat$secbase_observed[fulldat$time == 1], na.rm = TRUE) / n
    switch_observed[r] <- sum(fulldat$switch_status[fulldat$time == 1], na.rm = TRUE) / length(id_con)
    Minb0[r] <- sum(fulldat$secbase_observed[fulldat$time == 1 & fulldat$arm == 0 & fulldat$b1==0])
    Minb1[r] <- sum(fulldat$secbase_observed[fulldat$time == 1 & fulldat$arm == 0 & fulldat$b1 == 1])

    if(verbose > 1){
      system(paste("echo 'now processing: repetition",r,"'"))
      print(paste("repetition: ", r))
    }


     }# end giant loop. Build plots below.
  )







  # build boxplot dataset
  df <- data.frame(Unbiased = rep(NA, reps), ITT = NA, PP = NA, IPCW = NA, RPSFTM = NA, TSE = NA)
  df$Unbiased[1:length(unbiased)] <- unbiased
  df$ITT[1:length(itt)] <- itt
  df$PP[1:length(pp)] <- pp
  df$IPCW[1:length(msm_hr)] <-  msm_hr # TODO msm_hr no longer exists as a list, its now a vector
  df$RPSFTM[1:length(rpsft)] <- rpsft
  df$TSE[1:length(tse)] <- tse
  df <- df %>% pivot_longer(names(df), names_to = "Method", values_to = "est") # df wide to long, for ggplot
  df$bias <- df$est - mean(df$est[df$Method=="Unbiased"])
  compar_plot <- ggplot(df, aes(Method, est, color = Method)) + geom_boxplot() +
    scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.title.x=element_blank(),
                                                               axis.text.x=element_blank(),
                                                               axis.ticks.x=element_blank()) +
    ylab("HR Estimate") +
    ggtitle("HR Estimates Across Methods") +
    geom_hline(yintercept=stats::median(unbiased), linetype="dashed", size=0.1, alpha = 0.5)


  # produce point-comparison estimate of bias.
  bias_plot <- ggplot(df[df$Method != "Unbiased",], aes(Method, bias, color = Method)) + geom_boxplot() +
    scale_color_brewer(palette = "Dark2") + theme_bw() + theme(axis.title.x=element_blank(),
                                                               axis.text.x=element_blank(),
                                                               axis.ticks.x=element_blank()) +
    ylab("HR Estimate Bias") +
    ggtitle("HR Estimate Bias Across Methods") +
    geom_hline(yintercept=0, linetype="dashed", size=0.1, alpha = 0.5)



  if(verbose > 1){
    print("Done!")
  }


  return(new("SwitchSimulation",
             unbiased = unbiased,
             unbiased_plot = unbiased_plot,
             itt = itt,
             itt_plot = itt_plot,
             pp = pp,
             pp_plot = pp_plot,
             ipcw = msm_hr,
             ipcw_plot = ipcw_plot,
             rpsft = rpsft,
             rpsft_plot = rpsft_plot,
             tse = tse,
             tse_plot = tse_plot,
             compar_plot = compar_plot,
             bias_plot = bias_plot,
             params = list(violate = violate, cens_flag = cens_flag, prop_cens = prop_cens,
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
                           compar_df = df,
                           Minb0 = Minb0,
                           Minb1 = Minb1)
             ))


}


