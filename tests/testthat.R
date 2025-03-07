Sys.setenv(R_TESTS="")

library(testthat)
library(MOCHA)
test_check("MOCHA")
