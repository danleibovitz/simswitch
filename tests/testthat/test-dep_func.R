# dat = fulldat | data.frame, 40,000 x 12, names(fulldat) ==  [1] "ids"    "arm"    "switch" "time"   "treat"  "b1" "b2" "b3" "M" "v1" "v2" "Mtime"
# window = i, 1
# base_var=bcov_names, "b1" "b2" "b3"
# time_var=tcov_names, "M"  "v1" "v2"
# covar_coef = covar_coef, list(baseline, varying). baseline 3x3 (colsums = c(1,1,1)), varying 4x3
# m_haz = m_haz, vector 1x100, all positive
# num_t = num_tvar, 3
# tcov_n = tcov_names, "M"  "v1" "v2"

# TODO store a 'fulldat' dataset for running tests on.

# test that incorrect argument types throw errors ####
test_that("multiplication works", {
  expect_error(dep_func())
})

# test incorrect argument dimensions throw errors ####

# test that results of correct type and dimension are returned ####


