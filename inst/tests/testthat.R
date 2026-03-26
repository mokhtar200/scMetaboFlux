# This file is part of the standard testthat setup.

library(testthat)

# Determine package root
pkg_root <- find.package("scMetaboFlux")[1]
helpers_path <- file.path(pkg_root, "tests", "testthat", "helpers.R")

# Load helper functions
if (file.exists(helpers_path)) {
  source(helpers_path, local = TRUE)
}

# Load scMetaboFlux
library(scMetaboFlux)

# Set test options
testthat::test_dir(
  file.path(pkg_root, "tests", "testthat"),
  filter = "test-",
  stop_on_failure = FALSE,
  stop_on_warning = FALSE,
  reporter = testthat::default_reporter()
)
