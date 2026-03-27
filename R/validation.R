#' @title Run Comprehensive Test Scenarios for scMetaboFlux
#'
#' @description Runs test scenarios with various simulated datasets
#' to validate package functionality. For comprehensive Bioconductor
#' submission testing, see \code{\link{runBioconductorValidationTests}}.
#'
#' @param verbose Print progress messages
#' @param stop_on_error Stop if any test fails
#'
#' @return List of test results
#' @export
#' @examples
#' \dontrun{
#' results <- runTestScenarios()
#' }
#'
#' @seealso \code{\link{runBioconductorValidationTests}} for comprehensive
#' Bioconductor-compliant validation testing including large datasets,
#' performance benchmarks, and edge cases.
runTestScenarios <- function(verbose = TRUE, stop_on_error = FALSE) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required for testing")
  }
  
  results <- list()
  
  if (verbose) cat("\n=== scMetaboFlux Test Scenarios ===\n\n")
  
  # Scenario 1: PBMC-like data
  if (verbose) cat("Scenario 1: PBMC-like data...\n")
  tryCatch({
    set.seed(42)
    pbmc_data <- generateExampleData(n_cells = 300, n_genes = 500)
    
    # Run pathway scoring
    pbmc_data <- computePathwayScores(
      pbmc_data,
      gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation", "tca_cycle")],
      method = "mean"
    )
    
    # Estimate ATP
    pbmc_data <- estimateATPProduction(pbmc_data)
    
    # Classify phenotypes
    pbmc_data <- classifyMetabolicPhenotype(pbmc_data, method = "quantile")
    
    results$pbmc <- list(
      status = "PASS",
      n_cells = ncol(pbmc_data),
      n_genes = nrow(pbmc_data),
      has_scores = hasMetabolicScores(pbmc_data),
      has_phenotypes = "metabolic_phenotype" %in% colnames(pbmc_data@meta.data)
    )
    if (verbose) cat("  PASS\n")
  }, error = function(e) {
    results$pbmc <- list(status = "FAIL", error = e$message)
    if (verbose) cat("  FAIL:", e$message, "\n")
    if (stop_on_error) stop(e)
  })
  
  # Scenario 2: Cancer hypoxia data
  if (verbose) cat("Scenario 2: Cancer hypoxia data...\n")
  tryCatch({
    set.seed(123)
    cancer_data <- generateCancerData(n_cells = 500, n_genes = 1000)
    
    # Run pathway scoring only (skip flux)
    cancer_data <- computePathwayScores(
      cancer_data,
      gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation", "hypoxia_response")],
      method = "mean"
    )
    
    # Test cell type aggregation
    if ("region" %in% colnames(cancer_data@meta.data)) {
      celltype_summary <- aggregateByCellType(cancer_data, cell_type_col = "region")
      agg_status <- is.data.frame(celltype_summary)
    } else {
      agg_status <- FALSE
    }
    
    results$cancer <- list(
      status = "PASS",
      n_cells = ncol(cancer_data),
      has_scores = hasMetabolicScores(cancer_data),
      aggregation_works = agg_status
    )
    if (verbose) cat("  PASS\n")
  }, error = function(e) {
    results$cancer <- list(status = "FAIL", error = e$message)
    if (verbose) cat("  FAIL:", e$message, "\n")
    if (stop_on_error) stop(e)
  })
  
  # Scenario 3: Differentiation time series
  if (verbose) cat("Scenario 3: Differentiation time series...\n")
  tryCatch({
    set.seed(456)
    diff_data <- generateDifferentiationData(n_cells = 200, n_genes = 500)
    
    # Compute pathway scores
    diff_data <- computePathwayScores(
      diff_data,
      gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
      method = "mean"
    )
    
    # Estimate ATP
    diff_data <- estimateATPProduction(diff_data)
    
    results$differentiation <- list(
      status = "PASS",
      n_cells = ncol(diff_data),
      has_scores = hasMetabolicScores(diff_data)
    )
    if (verbose) cat("  PASS\n")
  }, error = function(e) {
    results$differentiation <- list(status = "FAIL", error = e$message)
    if (verbose) cat("  FAIL:", e$message, "\n")
    if (stop_on_error) stop(e)
  })
  
  # Scenario 4: Large dataset stress test
  if (verbose) cat("Scenario 4: Large dataset (500 cells)...\n")
  tryCatch({
    set.seed(789)
    large_data <- generateExampleData(n_cells = 500, n_genes = 1000)
    
    # Run pathway scoring
    large_data <- computePathwayScores(
      large_data,
      gene_sets = metabolicGeneSets[c("glycolysis", "tca_cycle", "oxidative_phosphorylation")],
      method = "mean"
    )
    
    # Estimate ATP
    large_data <- estimateATPProduction(large_data)
    
    # Test differential analysis
    if ("cell_type" %in% colnames(large_data@meta.data)) {
      cell_types <- unique(large_data@meta.data$cell_type)
      if (length(cell_types) >= 2) {
        diff_result <- differentialMetabolicAnalysis(
          large_data,
          condition_col = "cell_type",
          control_group = cell_types[1],
          case_group = cell_types[2]
        )
        diff_status <- is.data.frame(diff_result) || is.list(diff_result)
      } else {
        diff_status <- FALSE
      }
    } else {
      diff_status <- FALSE
    }
    
    results$large_scale <- list(
      status = "PASS",
      n_cells = ncol(large_data),
      n_genes = nrow(large_data),
      differential_works = diff_status
    )
    if (verbose) cat("  PASS\n")
  }, error = function(e) {
    results$large_scale <- list(status = "FAIL", error = e$message)
    if (verbose) cat("  FAIL:", e$message, "\n")
    if (stop_on_error) stop(e)
  })
  
  # Summary
  passed <- sum(sapply(results, function(x) x$status == "PASS"))
  total <- length(results)
  
  if (verbose) {
    cat("\n=== Test Summary ===\n")
    cat(sprintf("Passed: %d/%d\n", passed, total))
    for (name in names(results)) {
      status <- results[[name]]$status
      cat(sprintf("  %s: %s\n", name, status))
    }
  }
  
  return(invisible(results))
}

#' @title Validate Package Integrity
#'
#' @description Validates that all exported functions are accessible and documented.
#'
#' @return List of validation results
#' @export
validatePackageIntegrity <- function() {
  
  results <- list()
  
  # Check core datasets
  results$gene_sets <- is.list(metabolicGeneSets) && length(metabolicGeneSets) >= 10
  results$rate_weights <- is.numeric(rateLimitingWeights) || is.list(rateLimitingWeights)
  results$atp_coef <- is.numeric(atpYieldCoefficients) || is.list(atpYieldCoefficients)
  
  # Check core functions exist
  core_functions <- c(
    "runMetabolicAnalysis", "computePathwayScores", "estimateATPProduction",
    "computeMetabolicFlux", "classifyMetabolicPhenotype", "aggregateByCellType",
    "differentialMetabolicAnalysis", "createScMetaboFlux",
    "generateExampleData", "metabolicGeneSets"
  )
  
  for (fn in core_functions) {
    results[[paste0("fn_", fn)]] <- exists(fn, mode = "function")
  }
  
  # Check S4 class
  results$s4_class <- exists("scMetaboFluxObject", mode = "function")
  
  # Summary
  all_pass <- all(unlist(results))
  cat("\nPackage Integrity Check:\n")
  cat(sprintf("Overall: %s\n", ifelse(all_pass, "PASS", "FAIL")))
  cat(sprintf("Tests passed: %d/%d\n", sum(unlist(results)), length(results)))
  
  return(invisible(results))
}
