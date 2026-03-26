#!/usr/bin/env Rscript
# ============================================================================
# TEST RUNNER FOR scMetaboFlux
# ============================================================================
# Run all tests for the scMetaboFlux package
# Usage: Rscript run_tests.R [--coverage] [--filter PATTERN]
# ============================================================================

suppressPackageStartupMessages({
  library(testthat)
  library(methods)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
run_coverage <- "--coverage" %in% args

# Find package root
pkg_root <- if (file.exists("DESCRIPTION")) "." else ".."

# Change to package root
setwd(pkg_root)

# Load helper functions
source(file.path("tests", "testthat", "helpers.R"))

# Source scMetaboFlux
devtools::load_all(".", quiet = TRUE)

# Set test reporter
reporter <- testthat::default_reporter()

# Progress reporter for CI
if (Sys.getenv("CI") == "true") {
  reporter <- testthat::MinimalReporter$new()
}

# Run tests with coverage if requested
if (run_coverage) {
  if (requireNamespace("covr", quietly = TRUE)) {
    test_results <- covr::package_coverage(
      path = ".",
      type = "all",
      quiet = FALSE
    )
    
    cat("\n=== Coverage Report ===\n")
    print(covr::report(test_results))
    
    # Fail if coverage is too low
    total_coverage <- sum(test_results$covered) / sum(test_results$count) * 100
    if (total_coverage < 50) {
      stop("Coverage too low: ", round(total_coverage, 1), "%")
    }
  } else {
    cat("covr package not available, skipping coverage\n")
    test_dir("tests/testthat", stop_on_failure = FALSE)
  }
} else {
  # Run all tests
  cat("Running scMetaboFlux tests...\n\n")
  
  test_dir(
    "tests/testthat",
    stop_on_failure = FALSE,
    stop_on_warning = FALSE,
    reporter = reporter
  )
}

cat("\n=== Test Run Complete ===\n")
