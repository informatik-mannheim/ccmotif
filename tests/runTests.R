# Tests invoked from devtools:check()
library(testthat)

# test_results = test_dir("tests/testthat", reporter = "summary")
#test_results = test_dir("tests/testthat", reporter = LocationReporter)
test_results = test_dir("testthat", reporter = "summary")
test_results