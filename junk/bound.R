# function for bounding a vector

bound <- function(x, minv, maxv){

  return(sapply(x, function(y) min(max(y,minv),maxv)))

}
