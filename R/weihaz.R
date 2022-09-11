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
weihaz <- function(
    scale,
    shape,
    x) {

  # Defend against incorrect argument types
  if (class(shape) != "numeric" & class(shape) != "integer") stop("shape parameter must be numeric")
  if (class(scale) != "numeric" & class(scale) != "integer") stop("scale parameter must be numeric")
  if (class(x) != "numeric" & class(x) != "integer") stop("x parameter must be numeric")

  # Defend against incorrect argument dimensions


  # Defend against incorrect argument magnitudes
  if (shape <= 0) stop("shape parameter must be positive")
  if (scale <= 0) stop("scale parameter must be positive")
  if (range(x)[1] <= 0) stop("x parameter must all be positive")

  return((shape/scale)*(x/scale)^(shape-1))
}
