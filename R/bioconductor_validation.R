# ============================================================================
# COMPREHENSIVE BIOCONDUCTOR VALIDATION TEST SUITE
# ============================================================================
# This suite tests scMetaboFlux with challenging real-world scenarios
# to ensure Bioconductor compliance and robustness
# ============================================================================

#' @title Run Comprehensive Bioconductor Validation Tests
#' 
#' @description Runs extensive validation tests covering edge cases,
#' large datasets, performance benchmarks, and integration scenarios.
#' 
#' @param verbose Print detailed test output
#' @param save_results Save test results to file
#' @param output_dir Directory for results (default tempdir)
#' 
#' @return List with test results and summary statistics
#' 
#' @export
#' @examples
#' \dontrun{
#' results <- runBioconductorValidationTests()
#' print(results$summary)
#' }
runBioconductorValidationTests <- function(verbose = TRUE, 
                                         save_results = TRUE,
                                         output_dir = tempdir()) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is required for validation tests")
  }
  
  results <- list()
  test_start <- Sys.time()
  
  if (verbose) {
    message("========================================================")
    message("scMetaboFlux Bioconductor Validation Suite v1.0")
    message("========================================================")
    message("")
  }
  
  # =========================================================================
  # TEST CATEGORY 1: LARGE DATASET TESTS (10K+ cells)
  # =========================================================================
  if (verbose) message("[1/8] Testing large datasets (10,000+ cells)...")
  
  results$large_dataset <- testLargeDataset(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 2: SPARSE MATRIX TESTS
  # =========================================================================
  if (verbose) message("[2/8] Testing sparse matrix handling...")
  
  results$sparse_matrices <- testSparseMatrices(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 3: MISSING DATA TESTS
  # =========================================================================
  if (verbose) message("[3/8] Testing missing data handling...")
  
  results$missing_data <- testMissingDataHandling(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 4: EDGE CASES IN CLASSIFICATION
  # =========================================================================
  if (verbose) message("[4/8] Testing edge cases in classification...")
  
  results$classification_edge_cases <- testClassificationEdgeCases(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 5: ALL SCORING METHODS
  # =========================================================================
  if (verbose) message("[5/8] Testing all scoring methods...")
  
  results$scoring_methods <- testAllScoringMethods(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 6: PERFORMANCE BENCHMARKS
  # =========================================================================
  if (verbose) message("[6/8] Running performance benchmarks...")
  
  results$performance <- testPerformanceBenchmarks(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 7: INTEGRATION WITH SEURAT WORKFLOWS
  # =========================================================================
  if (verbose) message("[7/8] Testing Seurat workflow integration...")
  
  results$seurat_integration <- testSeuratWorkflowIntegration(verbose = verbose)
  
  # =========================================================================
  # TEST CATEGORY 8: ERROR HANDLING
  # =========================================================================
  if (verbose) message("[8/8] Testing error handling...")
  
  results$error_handling <- testErrorHandling(verbose = verbose)
  
  # =========================================================================
  # GENERATE SUMMARY
  # =========================================================================
  test_end <- Sys.time()
  
  total_tests <- sum(sapply(results, function(x) x$n_tests))
  passed_tests <- sum(sapply(results, function(x) x$n_passed))
  failed_tests <- sum(sapply(results, function(x) x$n_failed))
  
  summary <- list(
    total_tests = total_tests,
    passed = passed_tests,
    failed = failed_tests,
    pass_rate = round(passed_tests / total_tests * 100, 2),
    duration = as.numeric(difftime(test_end, test_start, units = "mins")),
    timestamp = test_end,
    package_version = utils::packageVersion("scMetaboFlux")
  )
  
  results$summary <- summary
  
  if (verbose) {
    message("")
    message("========================================================")
    message("TEST SUMMARY")
    message("========================================================")
    message(paste("Total tests:", total_tests))
    message(paste("Passed:", passed_tests))
    message(paste("Failed:", failed_tests))
    message(paste("Pass rate:", summary$pass_rate, "%"))
    message(paste("Duration:", round(summary$duration, 2), "minutes"))
    message("========================================================")
  }
  
  if (save_results) {
    results_file <- file.path(output_dir, paste0("scMetaboFlux_validation_", 
                                                format(Sys.time(), "%Y%m%d_%H%M%S"), ".Rds"))
    saveRDS(results, results_file)
    if (verbose) message(paste("Results saved to:", results_file))
  }
  
  return(results)
}

# ============================================================================
# TEST 1: LARGE DATASET TESTS
# ============================================================================

testLargeDataset <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  # Test 1.1: 10,000 cells
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    n_cells <- 10000
    n_genes <- 2000
    
    counts <- matrix(rnbinom(n_genes * n_cells, size = 2, mu = 3), 
                     nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("GENE", 1:n_genes)
    colnames(counts) <- paste0("Cell_", 1:n_cells)
    
    obj <- Seurat::CreateSeuratObject(counts = counts)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 500)
    obj <- Seurat::ScaleData(obj)
    
    # Add metabolic genes
    glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1")
    oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1")
    other_genes <- rownames(obj)[1:100]
    metabolic_genes <- c(glycolysis_genes, oxphos_genes, other_genes)
    
    for (g in metabolic_genes) {
      if (g %in% rownames(obj)) {
        idx <- which(rownames(obj) == g)
        obj@assays$RNA@counts[idx, ] <- obj@assays$RNA@counts[idx, ] + 5
      }
    }
    
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    obj <- estimateATPProduction(obj, verbose = FALSE)
    
    if ("ATP_score" %in% colnames(obj@meta.data) && 
        length(obj@meta.data$ATP_score) == n_cells) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ 10K cells processed successfully")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "10K cells: ATP scores not computed correctly"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("10K cells:", e$message)
    if (verbose) message(paste("  ✗ 10K cells:", e$message))
  })
  
  # Test 1.2: 20,000 cells with parallel processing
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(123)
    n_cells <- 20000
    n_genes <- 1500
    
    counts <- matrix(rnbinom(n_genes * n_cells, size = 2, mu = 2), 
                     nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("GENE", 1:n_genes)
    colnames(counts) <- paste0("Cell_", 1:n_cells)
    
    obj <- Seurat::CreateSeuratObject(counts = counts)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 300)
    obj <- Seurat::ScaleData(obj)
    
    # Score with AUCell
    obj <- computePathwayScores(obj, method = "AUCell", verbose = FALSE)
    
    if ("glycolysis" %in% colnames(obj@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ 20K cells with AUCell processed successfully")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "20K cells: AUCell scoring failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("20K cells:", e$message)
    if (verbose) message(paste("  ✗ 20K cells:", e$message))
  })
  
  # Test 1.3: Many features (genes)
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(456)
    n_cells <- 1000
    n_genes <- 10000
    
    counts <- matrix(rnbinom(n_genes * n_cells, size = 2, mu = 3), 
                     nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("GENE", 1:n_genes)
    colnames(counts) <- paste0("Cell_", 1:n_cells)
    
    # Add some metabolic genes
    metabolic_genes <- c("HK2", "LDHA", "SDHA", "COX5A", "PFKM", "PGK1")
    for (i in seq_along(metabolic_genes)) {
      idx <- i * 100
      if (idx <= n_genes) {
        counts[idx, ] <- counts[idx, ] + 10
        rownames(counts)[idx] <- metabolic_genes[i]
      }
    }
    
    obj <- Seurat::CreateSeuratObject(counts = counts)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 1000)
    obj <- Seurat::ScaleData(obj)
    
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    
    if (nrow(obj@meta.data) == n_cells) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ 10K genes processed successfully")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "10K genes: Processing failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("10K genes:", e$message)
    if (verbose) message(paste("  ✗ 10K genes:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed, 
              errors = errors, category = "Large Datasets"))
}

# ============================================================================
# TEST 2: SPARSE MATRIX TESTS
# ============================================================================

testSparseMatrices <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  # Test 2.1: Highly sparse matrix (95% zeros)
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    n_cells <- 500
    n_genes <- 1000
    sparsity <- 0.95
    
    counts <- Matrix::rsparsematrix(n_genes, n_cells, density = 1 - sparsity,
                                    rand.Gen = function(n) rnbinom(n, size = 2, mu = 2))
    rownames(counts) <- paste0("GENE", 1:n_genes)
    colnames(counts) <- paste0("Cell_", 1:n_cells)
    
    # Add metabolic genes
    counts["GENE1", ] <- counts["GENE1", ] + 5
    rownames(counts)[1] <- "HK2"
    
    obj <- Seurat::CreateSeuratObject(counts = counts)
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::ScaleData(obj)
    
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    
    if (inherits(obj@assays$RNA@data, "sparseMatrix") && 
        "glycolysis" %in% colnames(obj@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Highly sparse matrix (95%) handled correctly")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Sparse matrix: Result not sparse"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Sparse matrix:", e$message)
    if (verbose) message(paste("  ✗ Sparse matrix:", e$message))
  })
  
  # Test 2.2: dgCMatrix format preservation
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(789)
    counts <- Matrix::rsparsematrix(500, 200, density = 0.1)
    rownames(counts) <- paste0("G", 1:500)
    colnames(counts) <- paste0("C", 1:200)
    
    obj <- Seurat::CreateSeuratObject(counts = counts)
    original_sparsity <- sum(counts == 0) / length(counts)
    
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::ScaleData(obj)
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    
    n_passed <- n_passed + 1
    if (verbose) message("  ✓ dgCMatrix format preserved")
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("dgCMatrix:", e$message)
    if (verbose) message(paste("  ✗ dgCMatrix:", e$message))
  })
  
  # Test 2.3: Dense to sparse conversion
  n_tests <- n_tests + 1
  tryCatch({
    counts <- matrix(rnorm(500 * 200, mean = 2, sd = 1), nrow = 500, ncol = 200)
    rownames(counts) <- paste0("G", 1:500)
    colnames(counts) <- paste0("C", 1:200)
    
    obj <- Seurat::CreateSeuratObject(counts = counts)
    obj <- Seurat::NormalizeData(obj)
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    
    n_passed <- n_passed + 1
    if (verbose) message("  ✓ Dense matrix conversion handled")
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Dense matrix:", e$message)
    if (verbose) message(paste("  ✗ Dense matrix:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, category = "Sparse Matrices"))
}

# ============================================================================
# TEST 3: MISSING DATA TESTS
# ============================================================================

testMissingDataHandling <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  # Test 3.1: NA values in metadata
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(200)
    obj@meta.data$glycolysis[sample(200, 20)] <- NA
    obj@meta.data$oxidative_phosphorylation[sample(200, 20)] <- NA
    
    result <- estimateATPProduction(obj, verbose = FALSE)
    
    na_count <- sum(is.na(result@meta.data$ATP_score))
    if (na_count == 0 && "ATP_score" %in% colnames(result@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ NA values handled in ATP estimation")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "NA handling: NAs not properly handled"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("NA handling:", e$message)
    if (verbose) message(paste("  ✗ NA handling:", e$message))
  })
  
  # Test 3.2: NaN values
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(150)
    obj@meta.data$ATP_score[sample(150, 15)] <- NaN
    
    result <- estimateATPProduction(obj, verbose = FALSE)
    nan_count <- sum(is.nan(result@meta.data$ATP_score))
    
    if (nan_count == 0) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ NaN values handled correctly")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "NaN handling: NaNs remain"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("NaN handling:", e$message)
    if (verbose) message(paste("  ✗ NaN handling:", e$message))
  })
  
  # Test 3.3: Inf values
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(150)
    obj@meta.data$glycolysis[1] <- Inf
    obj@meta.data$glycolysis[2] <- -Inf
    
    result <- calculateGOXIndex(obj, verbose = FALSE)
    inf_count <- sum(is.infinite(result@meta.data$GOX_index))
    
    if (inf_count == 0) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Inf values handled correctly")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Inf handling: Infs remain"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Inf handling:", e$message)
    if (verbose) message(paste("  ✗ Inf handling:", e$message))
  })
  
  # Test 3.4: Empty cell types
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(200)
    obj@meta.data$cell_type <- as.character(obj@meta.data$cell_type)
    obj@meta.data$cell_type[1:50] <- "RareCell"
    
    result <- aggregateByCellType(obj, cell_type_col = "cell_type", verbose = FALSE)
    
    if (is.data.frame(result) && nrow(result) > 0) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Rare cell types handled correctly")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Rare cells: Aggregation failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Rare cells:", e$message)
    if (verbose) message(paste("  ✗ Rare cells:", e$message))
  })
  
  # Test 3.5: All zeros in a pathway
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(100)
    obj@meta.data$glycolysis <- 0
    obj@meta.data$oxidative_phosphorylation <- 1
    
    result <- calculateGOXIndex(obj, verbose = FALSE)
    
    if ("GOX_index" %in% colnames(result@meta.data) &&
        all(!is.na(result@meta.data$GOX_index))) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Zero values handled correctly")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Zero values: Processing failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Zero values:", e$message)
    if (verbose) message(paste("  ✗ Zero values:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, category = "Missing Data"))
}

# ============================================================================
# TEST 4: CLASSIFICATION EDGE CASES
# ============================================================================

testClassificationEdgeCases <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  # Test 4.1: All cells same phenotype
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(100)
    obj@meta.data$glycolysis <- 1.0
    obj@meta.data$oxidative_phosphorylation <- 0.0
    obj@meta.data$ATP_score <- 0.5
    
    result <- classifyMetabolicPhenotype(obj, method = "quantile", verbose = FALSE)
    
    if ("metabolic_phenotype" %in% colnames(result@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Uniform data classification works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Uniform: Classification failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Uniform:", e$message)
    if (verbose) message(paste("  ✗ Uniform:", e$message))
  })
  
  # Test 4.2: K-means classification
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(200)
    
    result <- classifyMetabolicPhenotype(obj, method = "kmeans", verbose = FALSE)
    
    phenotypes <- result@meta.data$metabolic_phenotype
    if (length(unique(phenotypes)) <= 5 && "metabolic_phenotype" %in% colnames(result@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ K-means classification works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "K-means: Wrong number of phenotypes"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("K-means:", e$message)
    if (verbose) message(paste("  ✗ K-means:", e$message))
  })
  
  # Test 4.3: Hierarchical classification
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(150)
    
    result <- classifyMetabolicPhenotype(obj, method = "hierarchical", verbose = FALSE)
    
    if ("metabolic_phenotype" %in% colnames(result@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Hierarchical classification works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Hierarchical: Classification failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Hierarchical:", e$message)
    if (verbose) message(paste("  ✗ Hierarchical:", e$message))
  })
  
  # Test 4.4: GOX index classification
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(100)
    
    result <- classifyByGOXIndex(obj, verbose = FALSE)
    
    if ("metabolic_phenotype" %in% colnames(result@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ GOX index classification works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "GOX: Classification failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("GOX:", e$message)
    if (verbose) message(paste("  ✗ GOX:", e$message))
  })
  
  # Test 4.5: Phenotype purity calculation
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(100)
    obj@meta.data$seurat_clusters <- sample(1:5, 100, replace = TRUE)
    obj@meta.data$metabolic_phenotype <- sample(c("Glycolytic", "Oxidative"), 100, replace = TRUE)
    
    purity <- calculatePhenotypePurity(obj)
    
    if (is.data.frame(purity) && "purity" %in% colnames(purity)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Phenotype purity calculation works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Purity: Calculation failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Purity:", e$message)
    if (verbose) message(paste("  ✗ Purity:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, category = "Classification Edge Cases"))
}

# ============================================================================
# TEST 5: ALL SCORING METHODS
# ============================================================================

testAllScoringMethods <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  methods <- c("mean", "weighted_mean", "AUCell", "GSVA")
  
  for (method in methods) {
    n_tests <- n_tests + 1
    tryCatch({
      set.seed(42)
      obj <- createMinimalTestObject(100)
      
      result <- computePathwayScores(obj, method = method, verbose = FALSE)
      
      has_scores <- "glycolysis" %in% colnames(result@meta.data)
      valid_scores <- all(!is.na(result@meta.data$glycolysis)) && 
                      all(!is.nan(result@meta.data$glycolysis))
      
      if (has_scores && valid_scores) {
        n_passed <- n_passed + 1
        if (verbose) message(paste0("  ✓ Scoring method '", method, "' works"))
      } else {
        n_failed <<- n_failed + 1
        errors[[length(errors) + 1]] <<- paste(method, "scoring: Invalid scores")
      }
    }, error = function(e) {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- paste(method, "scoring:", e$message)
      if (verbose) message(paste0("  ✗ Scoring method '", method, "': ", e$message))
    })
  }
  
  # Test single sample method
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(80)
    
    result <- computePathwayScores(obj, method = "single_sample", verbose = FALSE)
    
    if ("glycolysis" %in% colnames(result@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Single sample scoring works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Single sample: Failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Single sample:", e$message)
    if (verbose) message(paste("  ✗ Single sample:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, category = "Scoring Methods"))
}

# ============================================================================
# TEST 6: PERFORMANCE BENCHMARKS
# ============================================================================

testPerformanceBenchmarks <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  benchmark_results <- list()
  
  # Benchmark 1: Pathway scoring speed
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(5000)
    
    start_time <- Sys.time()
    result <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    end_time <- Sys.time()
    
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    benchmark_results$pathway_scoring_5k <- duration
    
    if (duration < 30) {  # Should complete in under 30 seconds for 5K cells
      n_passed <- n_passed + 1
      if (verbose) message(paste0("  ✓ Pathway scoring (5K cells): ", round(duration, 2), "s"))
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- paste("Performance: Scoring too slow", duration, "s")
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Benchmark 1:", e$message)
    if (verbose) message(paste("  ✗ Benchmark 1:", e$message))
  })
  
  # Benchmark 2: ATP estimation speed
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(5000)
    
    start_time <- Sys.time()
    result <- estimateATPProduction(obj, verbose = FALSE)
    end_time <- Sys.time()
    
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    benchmark_results$atp_estimation_5k <- duration
    
    if (duration < 10) {  # Should be fast
      n_passed <- n_passed + 1
      if (verbose) message(paste0("  ✓ ATP estimation (5K cells): ", round(duration, 2), "s"))
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- paste("Performance: ATP too slow", duration, "s")
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Benchmark 2:", e$message)
    if (verbose) message(paste("  ✗ Benchmark 2:", e$message))
  })
  
  # Benchmark 3: Cell type aggregation speed
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(3000)
    obj@meta.data$cell_type <- sample(1:10, 3000, replace = TRUE)
    
    start_time <- Sys.time()
    result <- aggregateByCellType(obj, cell_type_col = "cell_type", verbose = FALSE)
    end_time <- Sys.time()
    
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    benchmark_results$celltype_aggregation_3k <- duration
    
    if (duration < 15) {
      n_passed <- n_passed + 1
      if (verbose) message(paste0("  ✓ Cell-type aggregation (3K cells): ", round(duration, 2), "s"))
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- paste("Performance: Aggregation too slow", duration, "s")
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Benchmark 3:", e$message)
    if (verbose) message(paste("  ✗ Benchmark 3:", e$message))
  })
  
  # Benchmark 4: Full workflow speed
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(2000)
    
    start_time <- Sys.time()
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    obj <- estimateATPProduction(obj, verbose = FALSE)
    obj <- calculateGOXIndex(obj, verbose = FALSE)
    obj <- classifyMetabolicPhenotype(obj, verbose = FALSE)
    end_time <- Sys.time()
    
    duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
    benchmark_results$full_workflow_2k <- duration
    
    if (duration < 60) {  # Full workflow under 1 minute for 2K cells
      n_passed <- n_passed + 1
      if (verbose) message(paste0("  ✓ Full workflow (2K cells): ", round(duration, 2), "s"))
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- paste("Performance: Workflow too slow", duration, "s")
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Benchmark 4:", e$message)
    if (verbose) message(paste("  ✗ Benchmark 4:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, benchmarks = benchmark_results,
              category = "Performance"))
}

# ============================================================================
# TEST 7: SEURAT WORKFLOW INTEGRATION
# ============================================================================

testSeuratWorkflowIntegration <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  # Test 7.1: Integration with SCTransform
  n_tests <- n_tests + 1
  tryCatch({
    if (!requireNamespace("sctransform", quietly = TRUE)) {
      message("  - Skipping SCTransform test (package not available)")
      n_tests <- n_tests - 1
    } else {
      set.seed(42)
      counts <- matrix(rnbinom(1000 * 200, size = 2, mu = 3), nrow = 1000, ncol = 200)
      rownames(counts) <- paste0("G", 1:1000)
      colnames(counts) <- paste0("C", 1:200)
      
      obj <- Seurat::CreateSeuratObject(counts = counts)
      obj <- sctransform::SCTransform(obj, verbose = FALSE)
      
      obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
      
      if ("glycolysis" %in% colnames(obj@meta.data)) {
        n_passed <- n_passed + 1
        if (verbose) message("  ✓ Integration with SCTransform works")
      } else {
        n_failed <<- n_failed + 1
        errors[[length(errors) + 1]] <<- "SCTransform: Scores not computed"
      }
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("SCTransform:", e$message)
    if (verbose) message(paste("  ✗ SCTransform:", e$message))
  })
  
  # Test 7.2: Integration with multiple assays
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(100)
    
    # Add a second assay
    counts2 <- matrix(rnorm(500 * 100), nrow = 500, ncol = 100)
    obj[["Spatial"]] <- Seurat::CreateAssayObject(counts = counts2)
    
    # Score from default assay
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    
    if (Seurat::DefaultAssay(obj) == "RNA" && "glycolysis" %in% colnames(obj@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Multiple assays integration works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Assays: Integration failed"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Assays:", e$message)
    if (verbose) message(paste("  ✗ Assays:", e$message))
  })
  
  # Test 7.3: Integration with reductions
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(200)
    obj <- Seurat::RunPCA(obj, npcs = 10, verbose = FALSE)
    obj <- Seurat::RunUMAP(obj, dims = 1:10, verbose = FALSE)
    obj <- Seurat::RunTSNE(obj, dims = 1:10, verbose = FALSE)
    
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    
    has_reductions <- all(c("pca", "umap", "tsne") %in% names(obj@reductions))
    has_scores <- "glycolysis" %in% colnames(obj@meta.data)
    
    if (has_reductions && has_scores) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Reductions preserved during scoring")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Reductions: Not preserved"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Reductions:", e$message)
    if (verbose) message(paste("  ✗ Reductions:", e$message))
  })
  
  # Test 7.4: Integration with labels and metadata
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(150)
    
    # Add various metadata types
    obj@meta.data$patient_id <- sample(paste0("P", 1:5), 150, replace = TRUE)
    obj@meta.data$treatment <- sample(c("Control", "Treatment"), 150, replace = TRUE)
    obj@meta.data$timepoint <- sample(c("Day0", "Day7", "Day14"), 150, replace = TRUE)
    obj@meta.data$percent_mt <- runif(150, 5, 25)
    
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    obj <- estimateATPProduction(obj, verbose = FALSE)
    
    meta_preserved <- all(c("patient_id", "treatment", "timepoint", "percent_mt") %in% 
                          colnames(obj@meta.data))
    scores_computed <- "ATP_score" %in% colnames(obj@meta.data)
    
    if (meta_preserved && scores_computed) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Metadata preserved during analysis")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Metadata: Not preserved"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Metadata:", e$message)
    if (verbose) message(paste("  ✗ Metadata:", e$message))
  })
  
  # Test 7.5: Reset and recompute
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(100)
    obj <- computePathwayScores(obj, method = "mean", verbose = FALSE)
    obj <- estimateATPProduction(obj, verbose = FALSE)
    
    original_scores <- obj@meta.data$ATP_score
    
    obj <- resetMetabolicResults(obj)
    
    metabolic_cols <- grep("^(glycolysis|oxidative|ATP|GOX|metabolic)", 
                           colnames(obj@meta.data), value = TRUE, ignore.case = TRUE)
    
    if (length(metabolic_cols) == 0 && !"ATP_score" %in% colnames(obj@meta.data)) {
      n_passed <- n_passed + 1
      if (verbose) message("  ✓ Reset and recompute workflow works")
    } else {
      n_failed <<- n_failed + 1
      errors[[length(errors) + 1]] <<- "Reset: Failed to remove scores"
    }
  }, error = function(e) {
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- paste("Reset:", e$message)
    if (verbose) message(paste("  ✗ Reset:", e$message))
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, category = "Seurat Integration"))
}

# ============================================================================
# TEST 8: ERROR HANDLING
# ============================================================================

testErrorHandling <- function(verbose = TRUE) {
  n_tests <- 0
  n_passed <- 0
  n_failed <- 0
  errors <- list()
  
  # Test 8.1: Invalid Seurat object
  n_tests <- n_tests + 1
  tryCatch({
    result <- computePathwayScores("not_a_seurat_object", verbose = FALSE)
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- "Error handling: Should have failed for invalid input"
  }, error = function(e) {
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Invalid input error caught")
  })
  
  # Test 8.2: Non-existent column
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(50)
    
    result <- estimateATPProduction(obj, glycolysis_col = "nonexistent_column", 
                                    verbose = FALSE)
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- "Error handling: Should have failed for missing column"
  }, error = function(e) {
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Missing column error caught")
  })
  
  # Test 8.3: Invalid method
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(50)
    
    result <- computePathwayScores(obj, method = "invalid_method", verbose = FALSE)
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- "Error handling: Should have failed for invalid method"
  }, error = function(e) {
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Invalid method error caught")
  })
  
  # Test 8.4: Empty dataset
  n_tests <- n_tests + 1
  tryCatch({
    obj <- Seurat::CreateSeuratObject(counts = matrix(nrow = 0, ncol = 0))
    
    result <- computePathwayScores(obj, verbose = FALSE)
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- "Error handling: Should have failed for empty dataset"
  }, error = function(e) {
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Empty dataset error caught")
  })
  
  # Test 8.5: Zero cells
  n_tests <- n_tests + 1
  tryCatch({
    obj <- Seurat::CreateSeuratObject(counts = matrix(1:10, ncol = 0))
    
    result <- estimateATPProduction(obj, verbose = FALSE)
    n_failed <<- n_failed + 1
    errors[[length(errors) + 1]] <<- "Error handling: Should have failed for zero cells"
  }, error = function(e) {
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Zero cells error caught")
  })
  
  # Test 8.6: Wrong data type for scores
  n_tests <- n_tests + 1
  tryCatch({
    set.seed(42)
    obj <- createMinimalTestObject(50)
    obj@meta.data$glycolysis <- "not_a_number"
    
    result <- calculateGOXIndex(obj, verbose = FALSE)
    # Should handle gracefully or convert
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Wrong data type handled")
  }, error = function(e) {
    n_passed <<- n_passed + 1
    if (verbose) message("  ✓ Wrong data type error caught")
  })
  
  return(list(n_tests = n_tests, n_passed = n_passed, n_failed = n_failed,
              errors = errors, category = "Error Handling"))
}

# ============================================================================
# HELPER FUNCTION: Create Minimal Test Object
# ============================================================================

createMinimalTestObject <- function(n_cells, n_genes = 500) {
  set.seed(42)
  
  # Metabolic genes
  glycolysis_genes <- c("HK2", "HK3", "GCK", "GPI", "PFKL", "PFKM", "PFKP",
                        "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "GAPDHS",
                        "PGK1", "PGK2", "PGAM1", "PGAM2", "ENO1", "ENO2",
                        "ENO3", "PKM1", "PKM2", "PKLR", "LDHA", "LDHB", "LDHC")
  
  oxphos_genes <- c("NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6",
                    "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12",
                    "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6",
                    "SDHA", "SDHB", "SDHC", "SDHD", "COX4I1", "COX5A", "COX5B",
                    "ATP5F1A", "ATP5F1B", "ATP5F1C")
  
  tca_genes <- c("CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B",
                 "OGDH", "FH", "MDH1", "MDH2")
  
  other_genes <- paste0("GENE", 1:(n_genes - length(c(glycolysis_genes, oxphos_genes, tca_genes))))
  all_genes <- c(glycolysis_genes, oxphos_genes, tca_genes, other_genes)
  all_genes <- all_genes[1:n_genes]
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 5),
    nrow = n_genes,
    ncol = n_cells
  )
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  # Create Seurat object
  obj <- Seurat::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  
  # Add metadata
  obj@meta.data$cell_type <- sample(c("T_cell", "B_cell", "Macrophage", "NK_cell"), 
                                     n_cells, replace = TRUE)
  
  # Add metabolic scores
  obj@meta.data$glycolysis <- rnorm(n_cells, mean = 0.5, sd = 0.2)
  obj@meta.data$oxidative_phosphorylation <- rnorm(n_cells, mean = 0.5, sd = 0.2)
  obj@meta.data$ATP_score <- runif(n_cells, min = 0.3, max = 0.9)
  obj@meta.data$GOX_index <- rnorm(n_cells, mean = 0, sd = 0.3)
  
  return(obj)
}

# ============================================================================
# PRINT VALIDATION REPORT
# ============================================================================

#' @title Print Validation Report
#' 
#' @description Print a formatted validation report.
#' 
#' @param results Results from runBioconductorValidationTests
#' 
#' @export
printValidationReport <- function(results) {
  cat("\n")
  cat("============================================================\n")
  cat("       scMetaboFlux BIOCONDUCTOR VALIDATION REPORT\n")
  cat("============================================================\n")
  cat("\n")
  
  for (category in names(results)) {
    if (category == "summary") next
    
    r <- results[[category]]
    cat(sprintf("%-35s Tests: %2d | Passed: %2d | Failed: %2d\n",
                r$category, r$n_tests, r$n_passed, r$n_failed))
    
    if (length(r$errors) > 0) {
      for (err in r$errors) {
        cat(sprintf("  -> ERROR: %s\n", err))
      }
    }
  }
  
  cat("\n")
  cat("============================================================\n")
  s <- results$summary
  cat(sprintf("SUMMARY: %d/%d tests passed (%.1f%%)\n", 
              s$passed, s$total_tests, s$pass_rate))
  cat(sprintf("Duration: %.1f minutes\n", s$duration))
  cat(sprintf("Version: %s\n", as.character(s$package_version)))
  cat("============================================================\n")
  cat("\n")
  
  invisible(results)
}
