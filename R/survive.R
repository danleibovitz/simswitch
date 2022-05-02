# Discrete-time survival time generator
survive <- function(x, hazard, betas, ncov, stime, idvar, ids,
                    b_haz){
  n <- length(ids)
  store <- ids # create repository of patients with unobserved events
  df <- data.frame(ids = ids, eventtime = rep.int(0, n), status = rep.int(0, n)) # create response dataframe
  i <- 1 # set iterator

  # for all time windows: probability of failure in that window is hazard produced by haz() function
  while(any(df$eventtime == 0) & i <= stime){ # TODO should this be a while loop, and stop when all events are observed?

    # randomly assign event to remaining patients
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
