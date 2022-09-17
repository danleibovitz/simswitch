
# Define test arguments

# test that incorrect argument types throw errors ####
test_that("Incorrect 'dat' argument throws error", {
  expect_error(
    .cc_generate()
  )
})

# test incorrect argument dimensions throw errors ####

# test that fulldat contains correct variables ####
# test that results of correct type and dimension are returned ####
