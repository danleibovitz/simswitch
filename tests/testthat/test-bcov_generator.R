
# Define test arguments
nb <- 4
ds <- 0.5
md <- 3
st <- 10
en <- 8

# test that incorrect argument types throw errors ####
test_that("Incorrect 'nb' argument throws error", {
  expect_error(
    bcov_generator(num_bvar = "nb", diags = ds, middle = md, stime = st, n = en)
  )
})
test_that("Incorrect 'ds' argument throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = "ds", middle = md, stime = st, n = en)
  )
})
test_that("Incorrect 'md' argument throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = "md", stime = st, n = en)
  )
})
test_that("Incorrect 'st' argument throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = "st", n = en)
  )
})
test_that("Incorrect 'en' argument throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = "en")
  )
})

# test incorrect argument dimensions throw errors ####
test_that("Incorrect 'num_bvar' dimension throws error", {
  expect_error(
    bcov_generator(num_bvar = c(nb, nb), diags = ds, middle = md, stime = st, n = en)
  )
})
test_that("Incorrect 'diags' dimension throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = c(ds, ds), middle = md, stime = st, n = en)
  )
})
test_that("Incorrect 'middle' dimension throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = c(md, md), stime = st, n = en)
  )
})
test_that("Incorrect 'stime' dimension throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = c(st, st), n = en)
  )
})
test_that("Incorrect 'n' dimension throws error", {
  expect_error(
    bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = c(en, en))
  )
})


# test that results of correct type and dimension are returned ####
test_that("Correc function call returns correct class", {
  expect_equal(
    class(bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = en)),
    "data.frame"
  )
})
test_that("Correc function call returns correct dimension", {
  expect_equal(
    dim(bcov_generator(num_bvar = nb, diags = ds, middle = md, stime = st, n = en)),
    c(en*st, nb)
  )
})

