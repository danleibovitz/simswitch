#' Title
#'
#' @param x
#' @param shape
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
weihaz <- function(x, shape, scale){
  (shape/scale)*(x/scale)^(shape-1)
}
