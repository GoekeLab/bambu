Sys.setenv("R_TESTS" = "")

library(testthat)
library(bambu)

test_check("bambu")
