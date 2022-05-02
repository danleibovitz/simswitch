# # Set some variables:
#
# n <- 100
# k_covars <- 3
# start_means <- c(60, 80, 120) # length == k_covars
# end_means <- c(50, 90, 160) # length == k_covars
# start_spread <- c(5, 10, 10) # length == k_covars, all positive
# end_spread <- c(7, 10, 20) # length == k_covars, all positive
# rate <- 0.5 # bounded by 0, 1
# times <- 100 # positive, real
#
# # generate 1st time data
#
# df <- data.frame(matrix(nrow = n*times, ncol = k_covars))
# df[0:(n-1)*100 + 1, ] <- mvrnorm(n = 100, mu = start_means, Sigma = diag(start_spread))
#
#
# # construct a matrix for eigenvectors and eigenvalues
#
# evecs <- matrix(c(end_means, rnorm(k_covars*(k_covars - 1), 0, 10 )), ncol = k_covars)
# evals <- diag(c(1, rate, runif(k_covars-2, 0, rate)))
# transmat <- evecs %*% evals %*% solve(evecs)
# ptransmat <- transmat %^% 100
# eigen(ptransmat)
#
# test_end <- start_means %*% (transmat %^% 1000)
#
#

# override classunion to allow a NULL or Int class
setClassUnion("lmt", c("numeric", "NULL"))

# Define bridge object
Bridge <- setClass(Class = "Bridge",
                   slots = c(
                     data = "matrix",
                     time = "integer",
                     n = "integer",
                     uppr = "lmt",
                     lwr = "lmt",
                     trend = "numeric",
                     trt_trend = "numeric"
                   ))


# function for transforming brownian bridges into a logistic binary generator
binarize <- function(dat){
  # standardize all values
  dat <- (dat - min(dat)) / (max(dat) - min(dat))
  return(dat)
}


# a function for getting the data from a bridge object
get_bridge_data <- function(brij){

  # check class of bridge
  if(class(brij) != "Bridge") stop("treat() is a method for Bridge objects")

  return(brij@data)
}


# TODO add a treatment effect
# upper: an upper bound on the brownian bridge, about which any crossing values are just reflected in absolute value.
# lwr: a lower bound on the brownian bridge, about which any crossing values are just reflected in absolute value.
# trnsfrm: a transformation from R to a given Real subset, in lieu of setting an upper and/or lower bound
bb <- function(start, end, var, time, trt, uppr, lwr, trnsfrm){
  # optional arg handling
  if(missing(trt)){
    trt <- NULL
    } else{ # TODO there MUST be a trt value.
    # TODO trt must be... a vector of trt times?
  }
  if(length(time) != 1) stop("Length of TIME must be equal to 1") # check length of time
  if(time %% 1 == 0 & time > 0 ){
    time <- as.integer(time) # if time is a positive integer, cast it as integer
  } else{
    stop("TIME must be positive integer")
  }
  if(any(var <= 0)) stop("VAR values must all be positive.")
  if(missing(uppr)){uppr <- NULL}
  else{
    if(!is.null(uppr)){
      if(length(uppr) != 1) stop("Length of UPPR must be equal to 1") # check length of uppr
      if(any(uppr < start) | any(uppr < end)) stop("UPPR must be higher than all values in START and END")
      }
    }
  if(missing(lwr)){lwr <- NULL}
  else{
    if(!is.null(lwr)){
      if(length(lwr) != 1) stop("Length of LWR must be equal to 1") # check length of lwr
      if(any(lwr > start) | any(lwr > end)) stop("LWR must be lower than all values START and END")
    }

  }
  if(missing(trnsfrm)){
    trnsfrm <- NULL
  }else{
    if(!is.null(trnsfrm)){
      if(class(trnsfrm) != "function") stop("The class of trnsform must be function")
      }
    }
  # TODO START, END, VAR must all be same length. TIME, UPPR, LWR, must all be 1.
  if(length(start) != length(end)) stop("Lengths of START and END must be equal.")
  if(length(end) != length(var)) stop("Lengths of END and VAR must be equal.")
  # TRNSFRM must be function from R to subset of R
  n <- as.integer(length(start)) # define number of pts followed

  trend <- (end - start)/time # the slope of the drift, without treatment
  trt_trend <- trt/time # the slope of the drift, with treatment. This is an additive treatment effect

  # construct bridge
  bridge <- sapply(1:n, function(x) c(0, cumsum(rnorm(time-1, 0, var[x])))) # construct a standard brownian motion, with TIME steps and VAR variance. It must start at 0
  # Below: subtract t/T*X_T, to make this a [0,0] pinned brownian bridge. then, add start and subtract t/T*(end-start)
  # to make this a [start,end] pinned brownian bridge
  bridge <- bridge - sapply(1:n, function(x) seq(0,1,length.out=time))%*%diag(bridge[time,], nrow = dim(bridge)[2]) # create a brownian bridge pinned at [0,0]
  bridge <- t(t(bridge) + start) + sapply(1:n, function(x) seq(0,1,length.out=time))%*%diag(end-start, nrow = n) # move pins from [0,0] to [start, end]

  # apply treatment changes
  if(!is.null(trt)){}
  if(!is.null(uppr) & !is.null(lwr)){
    # TODO scale values within lwr-uppr range, but recommend the trnsfrm argument
  }else{
    if(!is.null(uppr)){
      bridge <- -abs(bridge - uppr) + uppr
    }
    if(!is.null(lwr)){
      bridge <- abs(bridge - lwr) + lwr
    }
  }
  if(!is.null(trnsfrm)){
    # apply transformation, given a transformation function
    bridge <- trnsfrm(bridge)
  }


  return(Bridge(data = bridge,
                time = time,
                n = n,
                uppr = uppr,
                lwr = lwr,
                trend = trend,
                trt_trend = trt_trend))
}



# a function for adjusting covariate values, according to a treatment vector and treatment effects
treat <- function(brij, trt_time, trt_vec){

  # check class of bridge
  if(class(brij) != "Bridge") stop("treat() is a method for Bridge objects")

  # check dimensions and values
  if(trt_time < 1 | trt_time > brij@time) stop("trt_time must be within range of follow-up time.")
  if(trt_time %% 1 != 0) stop("trt_time must be an integer")
  if(length(trt_vec) != brij@n) stop("Length of TRT_VEC must be equal to BRIJ$n, i.e., the number of patients in the Bridge object")
  if(any(trt_vec != 0 & trt_vec != 1)) stop("All values in TRT_VEC must be either 0 or 1") # if any of the trt_vec values are not 0 or 1, stop

  addmat <- brij@data # construct a matrix of equal size
  addmat[] <- 0 # fill all values with 0
  addmat[trt_time:brij@time, trt_vec == 1] <- 1 # replace all values in treated patients, at all times at and after trt_time, with 1s
  addmat <- t(t(addmat)*brij@trt_trend) # replace 1 values with trt_trend, specific to each patient
  brij@data <- brij@data + addmat # simply sum the matrices. This iterates each patient's treated value (and all subsequent values) by the patient-specific trt_trend
  # TODO need to adjust here based on brij uppr and lwr values

  return(brij)
}

test <- bb(start = c(10, 20, 30), end = c(34,28,15), var = c(0.2,0.3,1), trt = c(-30,54,-45), time = 100,
           trnsfrm = NULL, uppr = NULL, lwr = NULL)
testup <- treat(brij = test, trt_time = 50, trt_vec = c(0,0,0))
testup <- treat(brij = testup, trt_time = 60, trt_vec = c(1,1,1))
testup <- treat(brij = testup, trt_time = 70, trt_vec = c(1,1,0))

testm <- bb(c(10, 15), c(5, -40), c(0.5, 3), 200)
plot(1:100, test, type="l")

bb <- rbridge(frequency = 10000)
bbscale <- bb*50 + 100*seq(0,1,length.out = 10000)
plot(1:10000, bbscale, type="l")
var(bb*100)
var(bbscale)


plot(1:100, test@data[,1], type = "l")
lines(1:100, test@data[,2], type="l")
lines(1:100, test@data[,3], type="l")
lines(1:100, testup@data[,1], type="l", col="red")
lines(1:100, testup@data[,2], type="l", col="red")
lines(1:100, testup@data[,3], type="l", col="red")


num_tvar <- 3
starts <- c(-5, 0, 10)
var1s <- c(4,1,7)
gains <- c(8, 2, 20)
var2s <- c(2,5,7)

# create covariate bridges
gentimes <- function(num_tvar,
                     starts,
                     var1s ,
                     gains ,
                     var2s ){

}

