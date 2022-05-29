# define dep_func() function

#' Title
#'
#' @param dat
#' @param window
#' @param base_var
#' @param time_var
#' @param covar_coef
#' @param m_haz
#' @param num_t
#' @param tcov_n
#'
#' @return
#' @export
#'
#' @examples
dep_func <- function(dat, window, base_var, time_var, covar_coef, m_haz, num_t, tcov_n){ # should return a binary variable for M, and continuous for all other tvar
  if(window == 1){# if window is 1, only a function of baseline.
    retval <- (as.matrix(dat[dat$time == window, names(dat) %in% base_var]) %*% covar_coef$baseline) + MASS::mvrnorm(n = sum(dat$time == window), rep.int(0,num_t), diag(0.1,nrow = num_t))
    retval[,1] <- 0 # Set initial M values to 0
    retval <- cbind(retval, rep.int(0, length(retval[,1]))) # set initial Mtimes
    return(retval)
  }else{# if window is not 1, a function of baseline, previous tvar and previous treat
    retval <- (as.matrix(dat[dat$time == (window-1), names(dat) %in% c("treat", time_var)]) %*% covar_coef$varying) + MASS::mvrnorm(n = sum(dat$time==(window-1)), rep.int(0, num_t), diag(0.015, nrow = num_t)) # TODO sum(dat$time==(window-1)) can probably just be n
    # change column of metastatic disease to binary, based upon hazard function
    retval[,1] <- ifelse(dat$M[dat$time == window-1] == 1, 1, 0) # if the previous window M is 1, continue M
    retval[retval[,1] == 0, 1] <- (rbinom( n = length(retval[, 1]), size = 1,
                                           prob = 1 - exp(-exp(log(m_haz[window]) +
                                                                 as.matrix(dat[dat$time == window-1, names(dat) %in% c("treat", tcov_n)]) %*% covar_coef$varying[,1])) ))[retval[,1] == 0] # randomly assign the treat variable with probability
    retval <- cbind(retval, rep.int(0, length(retval[,1])))
    retval[,num_t + 1] <- ifelse(retval[,1] == 0, 0,  # add a column representing Mtime, and calculate Mtime.
                                    ifelse(dat$M[dat$time == window-1] == 0, 1,
                                           dat$Mtime[dat$time == window-1]+1))
    return(retval)
  }
}
