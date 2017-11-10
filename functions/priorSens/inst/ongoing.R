source("~/sw/R/util_functions.R")
my_sweave(file="test_routines_current", dir_file="~/projects/prior_sensitivity/priorSens/inst")

options(priorSens_test_computefinal = FALSE)

require(devtools)
require(testthat)
load_all("priorSens", reset = FALSE)
test_package("priorSens")


require(Rcpp)
require(compiler)
require(microbenchmark)
require(ggplot2)

sourceCpp("~/projects/prior_sensitivity/priorSens/src/test.cpp")