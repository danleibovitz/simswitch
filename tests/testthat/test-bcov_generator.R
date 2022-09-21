
# Define test arguments
nb <- 4
ds <- 0.5
md <- 3
st <- 10
en <- 8

# test that incorrect argument types throw errors ####
test_that("Incorrect 'b_haz' argument throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = en)
  )
})

# test incorrect argument dimensions throw errors ####
test_that("Incorrect 'b_haz' argument dimension throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = en)
  )
})

# test that results of correct type and dimension are returned ####
test_that("Correc function call returns correct class", {
  expect_equal(
    class(bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = en)),
    "data.frame"
  )
})
