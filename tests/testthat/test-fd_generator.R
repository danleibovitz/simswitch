
# Define test arguments

# if(class(n) != "integer") stop("'n' must be of class 'integer'")
# if(class(stime) != "integer") stop("'stime' must be of class 'integer'")
# if(class(prop_trt) != "numeric") stop("'prop_trt' must be of class 'numeric'")
# if(class(id_trt) != "integer") stop("'id_trt' must be of class 'integer'")
# if(class(num_bvar) != "integer") stop("'num_bvar' must be of class 'integer'")
# if(class(bcov) != "data.frame") stop("'bcov' must be of class 'data.frame")
# if(class(num_tvar) != "integer") stop("'num_tvar' must be of class 'integer'")

en            <- as.integer(10)
number_b      <- as.integer(4)
number_t      <- as.integer(4)
proportion_t  <- 0.4
estime        <- as.integer(12)
id_treat      <- sample(x = 1:en, size = en*proportion_t, replace = FALSE)
becove        <- data.frame(matrix(nrow = en*estime, ncol = number_b))

# test that incorrect argument types throw errors ####
test_that("Incorrect 'bcov' argument throws error", {
  expect_error(
    fd_generator(bcov = "becove", id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'id_trt' argument throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = "id_treat", n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'n' argument throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = "en", num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'num_bvar' argument throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = "number_b",
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'num_tvar' argument throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = "number_t", prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'prop_trt' argument throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = "proportion_t", stime = estime)
  )
})
test_that("Incorrect 'stime' argument throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = "estime")
  )
})


# test incorrect argument dimensions throw errors ####
test_that("Incorrect 'bcov' dimension throws error", {
  expect_error(
    fd_generator(bcov = becove[1:2,], id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'n' dimension throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = c(en,en), num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'num_bvar' dimension throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = c(number_b, number_b),
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'num_tvar' dimension throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = c(number_t, number_t), prop_trt = proportion_t, stime = estime)
  )
})
test_that("Incorrect 'prop_trt' dimension throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = c(proportion_t, proportion_t), stime = estime)
  )
})
test_that("Incorrect 'stime' dimension throws error", {
  expect_error(
    fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = c(estime, estime))
  )
})

# test that results of correct type and dimension are returned ####
test_that("Correct call to fd_generator returns data.frame", {
  expect_equal(
    class(fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                 num_tvar = number_t, prop_trt = proportion_t, stime = estime)),
    "data.frame"
  )
})
test_that("Correct call to fd_generator returns data.frame with correct variable names", {
  expect_true(
    all(c("ids", "arm", "switch", "time", "treat", "Mtime") %in%
      names(fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                       num_tvar = number_t, prop_trt = proportion_t, stime = estime)))
  )
})
test_that("Correct call to fd_generator returns data.frame of correct dimension", {
  expect_equal(
    dim(fd_generator(bcov = becove, id_trt = id_treat, n = en, num_bvar = number_b,
                       num_tvar = number_t, prop_trt = proportion_t, stime = estime)),
    c(estime*en, number_t + number_b + 6)
  )
})
