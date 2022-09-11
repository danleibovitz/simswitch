# Discrete-time survival time generator


#' Title
#'
#' @param x
#' @param hazard
#' @param betas
#' @param ncov
#' @param stime
#' @param idvar
#' @param ids
#' @param b_haz
#'
#' @return
#' @export
#'
#' @examples
survive <- function(
    x,
    hazard,
    betas,
    ncov,
    stime,
    idvar,
    ids,
    b_haz) {

  # sdat <- survive( x = fulldat, hazard = haz, betas = beta.mat, ncov = (num_bvar + num_tvar), stime = stime, idvar = "ids", ids = unique(fulldat$ids),
  #                 b_haz = b_haz)
  # Defend against incorrect argument types
  if(class(x) != "data.frame") stop("'x' must be of class 'data.frame'")
  if(class(hazard) != "function") stop("'hazard' must be of class 'function'")
  if(class(betas) != "data.frame") stop("'betas' must be of class 'data.frame'")
  if(class(ncov) != ) stop("'ncov' must be of class _")
  if(class(stime) != "integer") stop("'stime' must be of class 'integer'")
  if(class(idvar) != ) stop("'idvar' must be of class _")
  if(class(ids) != ) stop("'ids' must be of class _")
  if(class(b_haz) != "numeric") stop("'b_haz' must be of class 'numeric'")

  # Defend against incorrect argument dimensions
  if(dim(x) != ) stop("'x' must be of class _")
  if(dim(betas) != ) stop("'betas' must be of dim _")
  if(dim(ncov) != ) stop("'ncov' must be of dim _")
  if(length(stime) != 1) stop("'stime' must be of length 1")
  if(dim(idvar) != ) stop("'idvar' must be of dim _")
  if(dim(ids) != ) stop("'ids' must be of dim _")
  if(length(b_haz) != stime) stop("'b_haz' must be of length equal to \code{stime}")

  # Defend against hazard() function with incorrect arguments
  # TODO how do you check for the arguments a function takes?
  if(args(hazard) != c("t", "x", "betas", "b_haz", "ncov")) stop("'hazard' must accept arguments:
                                                                 't', 'x', 'betas', 'b_haz', 'ncov'")

  n <- length(ids)
  store <- ids # create repository of patients with unobserved events
  df <- data.frame(ids = ids, eventtime = rep.int(0, n), status = rep.int(0, n)) # create response dataframe
  i <- 1 # set iterator

  # for all time windows: probability of failure in that window is hazard produced by haz() function
  while(any(df$eventtime == 0) & i <= stime){ # TODO should this be a while loop, and stop when all events are observed?

    # randomly assign event to remaining patients
    df$eventtime[df$ids %in% store] <- ifelse(rbinom(n = length(store), size = 1, prob = hazard(
      t=i, x=x[x$ids %in% store,], betas = betas[betas$ids %in% store,], b_haz = b_haz, ncov = ncov
      )) == 1, i,0) # assign current time to for patients who observe in this window

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
