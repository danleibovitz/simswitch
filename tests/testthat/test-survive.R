#' @param x A data.frame resulting from a call to fd_generator()
#' @param hazard A hazard generating function. The default is haz_func(). The function must take arguments:
#' t, x, betas, b_haz, and ncov. For details, see the documentation of haz_func()
#' @param betas A data.frame with length equal to the length of 'x', and width equal to (ncov + 3)
#' @param ncov The number of covariates (the number of baseline covariates plus the number of time-varying covariates)
#' @param stime The number of follow-up times
#' @param idvar The name of the variable in 'x' which holds patient ids.
#' @param ids Patient ids.
#' @param b_haz A vector of "baseline hazards" for all patients, with length equal to 'stime'
#'

# Define test arguments
ex <- fd_generator()
h <- haz_func
bs <-
nc <- 4
st <- 10
idv <- "ids"
is <- unique(x$ids)
bh <- weihaz()

# test that incorrect argument types throw errors ####
test_that("Incorrect 'x' argument throws error", {
  expect_error(
    survive(x = "ex", hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'hazard' argument throws error", {
  expect_error(
    survive(x = ex, hazard = "h", betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'betas' argument throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = "bs", ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ncov' argument throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = "nc", stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'stime' argument throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = "st", idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'idvar' argument throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = "idv", ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ids' argument throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = "is", b_haz = bh)
  )
})
test_that("Incorrect 'b_haz' argument throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = "bh")
  )
})

# test incorrect argument dimensions throw errors ####
test_that("Incorrect 'ex' argument dimension throws error", {
  expect_error(
    survive(x = ex[1:3,], hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'betas' argument dimension throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs[1:3,], ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ncov' argument dimension throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = c(nc, nc), stime = st, idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'stime' argument dimension throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = c(st,st), idvar = idv, ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'idvar' argument dimension throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = c(idv, idv), ids = is, b_haz = bh)
  )
})
test_that("Incorrect 'ids' argument dimension throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is[1:3], b_haz = bh)
  )
})
test_that("Incorrect 'b_haz' argument dimension throws error", {
  expect_error(
    survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh[1:3])
  )
})

# test that results of correct type and dimension are returned ####
test_that("Incorrect 'b_haz' argument dimension throws error", {
  expect_equal(
    class(survive(x = ex, hazard = h, betas = bs, ncov = nc, stime = st, idvar = idv, ids = is, b_haz = bh)),
      "data.frame"
  )
})
