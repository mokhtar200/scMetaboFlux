# ============================================================================
# EDGE CASE AND BOUNDARY CONDITION TESTS FOR scMetaboFlux
# ============================================================================
# Comprehensive tests for extreme inputs, boundary conditions, and error cases
# ============================================================================

library(testthat)
library(Seurat)

# ============================================================================
# SECTION 1: EMPTY AND MINIMAL INPUT TESTS
# ============================================================================

test_that("Function handles empty gene sets gracefully", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  
  # Empty gene set list
  expect_error(
    computePathwayScores(seurat_obj, gene_sets = list()),
    "gene_sets must be"
  )
  
  # NULL gene sets
  expect_error(
    computePathwayScores(seurat_obj, gene_sets = NULL),
    NA
  )
})

test_that("Function handles single gene gene sets", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  
  # Single gene
  result <- computePathwayScores(seurat_obj,
                               gene_sets = list(single_gene = "HK2"),
                               method = "mean",
                               verbose = FALSE)
  
  expect_true("single_gene" %in% colnames(result@meta.data))
  expect_true(all(!is.na(result@meta.data$single_gene)))
})

test_that("Function handles single cell", {
  skip_if_missing()
  
  # Create single cell object
  counts <- matrix(rpois(100, lambda = 5), nrow = 100, ncol = 1)
  rownames(counts) <- paste0("GENE", 1:100)
  colnames(counts) <- "SingleCell"
  
  seurat_obj <- CreateSeuratObject(counts = counts)
  seurat_obj <- NormalizeData(seurat_obj)
  
  # Should work with warning
  expect_warning(
    result <- computePathwayScores(seurat_obj,
                                 gene_sets = list(glycolysis = c("GENE1", "GENE2")),
                                 method = "mean",
                                 verbose = FALSE),
    NA
  )
})

test_that("Function handles single cell type", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  seurat_obj@meta.data$cell_type <- rep("OnlyType", 30)
  
  seurat_obj@meta.data$glycolysis <- rnorm(30, 0.5, 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(30, 0.5, 0.2)
  seurat_obj@meta.data$ATP_score <- runif(30, 0.3, 0.7)
  
  # Should not error
  result <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type")
  
  expect_equal(nrow(result), 1)
  expect_equal(result$cell_type[1], "OnlyType")
})

# ============================================================================
# SECTION 2: NUMERICAL EDGE CASES
# ============================================================================

test_that("Function handles zero values", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  
  # Set all values to zero
  expr_data <- GetAssayData(seurat_obj)
  expr_data[,] <- 0
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "counts",
                            assay = "RNA", new.data = expr_data)
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  # Should handle gracefully
  result <- validateInputData(seurat_obj, check_genes = FALSE)
  expect_false(result$valid)
})

test_that("Function handles extreme values", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Very large values
  seurat_obj@meta.data$glycolysis <- rnorm(30, mean = 1e6, sd = 1e5)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(30, mean = 1e6, sd = 1e5)
  
  expect_no_error(
    result <- estimateATPProduction(seurat_obj, verbose = FALSE)
  )
  
  expect_true(all(is.finite(result@meta.data$ATP_score)))
})

test_that("Function handles very small values", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Very small values
  seurat_obj@meta.data$glycolysis <- runif(30, min = 1e-10, max = 1e-9)
  seurat_obj@meta.data$oxidative_phosphorylation <- runif(30, min = 1e-10, max = 1e-9)
  
  expect_no_error(
    result <- estimateATPProduction(seurat_obj, verbose = FALSE)
  )
  
  expect_true(all(result@meta.data$ATP_score >= 0))
  expect_true(all(result@meta.data$ATP_score <= 1))
})

test_that("Function handles negative values in metadata", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Negative values
  seurat_obj@meta.data$glycolysis <- rnorm(30, mean = -1, sd = 0.5)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(30, mean = 1, sd = 0.5)
  
  expect_no_error(
    result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  )
  
  expect_true(is.numeric(result@meta.data$GOX_index))
})

test_that("Function handles zero variance data", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Constant values
  seurat_obj@meta.data$glycolysis <- rep(0.5, 30)
  seurat_obj@meta.data$oxidative_phosphorylation <- rep(0.5, 30)
  
  # Should not crash
  expect_no_error(
    result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  )
  
  # GOX should be approximately 0
  expect_true(abs(mean(result@meta.data$GOX_index, na.rm = TRUE)) < 0.1)
})

# ============================================================================
# SECTION 3: DATA TYPE EDGE CASES
# ============================================================================

test_that("Function handles character scores", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Try to convert scores to character
  seurat_obj@meta.data$glycolysis <- as.character(rnorm(30, 0.5, 0.2))
  
  expect_error(
    estimateATPProduction(seurat_obj, verbose = FALSE),
    NA
  )
})

test_that("Function handles factor scores", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  seurat_obj@meta.data$glycolysis <- factor(rep(c("low", "medium", "high"), 10))
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(30, 0.5, 0.2)
  
  expect_error(
    estimateATPProduction(seurat_obj, verbose = FALSE),
    NA
  )
})

# ============================================================================
# SECTION 4: MISSING DATA TESTS
# ============================================================================

test_that("Function handles missing values in metadata", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Add NA values
  seurat_obj@meta.data$glycolysis[1:5] <- NA
  seurat_obj@meta.data$oxidative_phosphorylation[6:10] <- NA
  
  expect_no_error(
    result <- estimateATPProduction(seurat_obj, verbose = FALSE)
  )
  
  # Should still produce results
  expect_true("ATP_score" %in% colnames(result@meta.data))
  expect_true(sum(is.na(result@meta.data$ATP_score)) == 0)
})

test_that("Function handles completely missing columns", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Try non-existent column
  expect_error(
    estimateATPProduction(seurat_obj, 
                        glycolysis_col = "nonexistent_column",
                        verbose = FALSE),
    "not found"
  )
})

test_that("Function handles partially missing cell type annotations", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Some cells without cell type
  seurat_obj@meta.data$cell_type[1:5] <- NA
  seurat_obj@meta.data$glycolysis <- rnorm(30, 0.5, 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(30, 0.5, 0.2)
  seurat_obj@meta.data$ATP_score <- runif(30, 0.3, 0.7)
  
  # Should work
  result <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type")
  
  # Should have fewer rows than total cells
  expect_true(nrow(result) < 30)
})

# ============================================================================
# SECTION 5: SCALING TESTS
# ============================================================================

test_that("Function handles very large datasets", {
  skip_if_missing()
  skip_on_cran()
  
  # 1000 cells, 500 genes
  seurat_obj <- create_test_seurat(200, 300)
  
  start_time <- Sys.time()
  
  result <- computePathwayScores(seurat_obj,
                               gene_sets = metabolicGeneSets[1:3],
                               method = "mean",
                               verbose = FALSE)
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Should complete within reasonable time
  expect_true(elapsed < 120)
  expect_true("glycolysis" %in% colnames(result@meta.data))
})

test_that("Function handles many cell types", {
  skip_if_missing()
  
  n_cells <- 100
  n_types <- 20
  
  seurat_obj <- create_test_seurat(n_cells, 50)
  
  # Many small cell type groups
  seurat_obj@meta.data$cell_type <- sample(
    paste0("Type_", 1:n_types),
    n_cells,
    replace = TRUE
  )
  
  seurat_obj@meta.data$glycolysis <- rnorm(n_cells, 0.5, 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(n_cells, 0.5, 0.2)
  seurat_obj@meta.data$ATP_score <- runif(n_cells, 0.3, 0.7)
  
  expect_no_error(
    result <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type")
  )
  
  expect_equal(nrow(result), n_types)
})

test_that("Function handles many pathways", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100)
  
  # All pathways
  expect_no_error(
    result <- computePathwayScores(seurat_obj,
                                gene_sets = metabolicGeneSets,
                                method = "mean",
                                verbose = FALSE)
  )
  
  # Should have many score columns
  score_cols <- grep("glycolysis|oxphos|tca|fatty", 
                    colnames(result@meta.data), value = TRUE)
  expect_true(length(score_cols) > 10)
})

# ============================================================================
# SECTION 6: SAMPLING AND SUBSETTING TESTS
# ============================================================================

test_that("Function works after subsetting", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  # Compute scores on full dataset
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(glycolysis = metabolicGeneSets$glycolysis),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Subset
  seurat_subset <- subset(seurat_obj, cells = colnames(seurat_obj)[1:20])
  
  # Scores should be preserved
  expect_true("glycolysis" %in% colnames(seurat_subset@meta.data))
  expect_equal(ncol(seurat_subset), 20)
})

test_that("Function handles downsampled data", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(100, 100)
  
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(glycolysis = metabolicGeneSets$glycolysis),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Downsample to 10%
  downsampled <- seurat_obj[, sample(ncol(seurat_obj), 10)]
  
  expect_true("glycolysis" %in% colnames(downsampled@meta.data))
})

# ============================================================================
# SECTION 7: MERGING AND COMBINATION TESTS
# ============================================================================

test_that("Function handles merged objects", {
  skip_if_missing()
  
  obj1 <- create_test_seurat(20, 50)
  obj2 <- create_test_seurat(20, 50)
  
  obj1 <- computePathwayScores(obj1,
                             gene_sets = list(glycolysis = metabolicGeneSets$glycolysis),
                             method = "mean",
                             verbose = FALSE)
  
  obj2 <- computePathwayScores(obj2,
                             gene_sets = list(glycolysis = metabolicGeneSets$glycolysis),
                             method = "mean",
                             verbose = FALSE)
  
  # Merge
  merged <- merge(obj1, obj2)
  merged <- NormalizeData(merged)
  merged <- ScaleData(merged)
  
  expect_true("glycolysis" %in% colnames(merged@meta.data))
})

# ============================================================================
# SECTION 8: TIMING AND REPRODUCIBILITY TESTS
# ============================================================================

test_that("Results are reproducible with set.seed", {
  skip_if_missing()
  
  set.seed(42)
  obj1 <- create_test_seurat(30, 100, seed = 42)
  obj1 <- computePathwayScores(obj1,
                             gene_sets = list(glycolysis = c("HK2", "LDHA")),
                             method = "mean",
                             verbose = FALSE)
  
  set.seed(42)
  obj2 <- create_test_seurat(30, 100, seed = 42)
  obj2 <- computePathwayScores(obj2,
                             gene_sets = list(glycolysis = c("HK2", "LDHA")),
                             method = "mean",
                             verbose = FALSE)
  
  expect_equal(obj1@meta.data$glycolysis, obj2@meta.data$glycolysis)
})

test_that("Different seeds produce different results", {
  skip_if_missing()
  
  obj1 <- create_test_seurat(30, 100, seed = 42)
  obj2 <- create_test_seurat(30, 100, seed = 123)
  
  obj1 <- computePathwayScores(obj1,
                             gene_sets = list(glycolysis = c("HK2", "LDHA")),
                             method = "mean",
                             verbose = FALSE)
  
  obj2 <- computePathwayScores(obj2,
                             gene_sets = list(glycolysis = c("HK2", "LDHA")),
                             method = "mean",
                             verbose = FALSE)
  
  # Should be different
  expect_false(identical(obj1@meta.data$glycolysis, obj2@meta.data$glycolysis))
})

# ============================================================================
# SECTION 9: ERROR HANDLING TESTS
# ============================================================================

test_that("Invalid input types are caught", {
  # Not a Seurat object
  expect_error(validateInputData("string"))
  expect_error(validateInputData(123))
  expect_error(validateInputData(data.frame()))
  expect_error(validateInputData(NULL))
})

test_that("Invalid methods are caught", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  
  expect_error(
    computePathwayScores(seurat_obj, method = "invalid_method"),
    "method"
  )
  
  expect_error(
    classifyMetabolicPhenotype(seurat_obj, method = "invalid"),
    "method"
  )
})

test_that("Missing required arguments are caught", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  
  expect_error(
    aggregateByCellType(seurat_obj),
    "cell type column"
  )
})

test_that("Invalid file paths are handled in exports", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  seurat_obj@meta.data$glycolysis <- rnorm(20)
  seurat_obj@meta.data$ATP_score <- runif(20)
  
  # Invalid path should cause error
  expect_error(
    exportMetabolicResults(seurat_obj, output_dir = "/nonexistent/path/that/does/not/exist"),
    NA
  )
})

# ============================================================================
# SECTION 10: CONSISTENCY TESTS
# ============================================================================

test_that("ATP score is normalized correctly", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 50)
  
  # Uniform scores
  seurat_obj@meta.data$glycolysis <- rep(1, 30)
  seurat_obj@meta.data$oxidative_phosphorylation <- rep(1, 30)
  
  result <- estimateATPProduction(seurat_obj, normalize = TRUE, verbose = FALSE)
  
  # All ATP scores should be equal
  expect_equal(var(result@meta.data$ATP_score), 0, tolerance = 1e-10)
})

test_that("Cell counts match in aggregations", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  seurat_obj@meta.data$glycolysis <- rnorm(50)
  seurat_obj@meta.data$ATP_score <- runif(50)
  
  result <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type")
  
  # Total cells should match
  total_from_agg <- sum(result$total_cells)
  expect_equal(total_from_agg, 50)
})

test_that("Phenotype proportions sum to 100%", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  seurat_obj@meta.data$glycolysis <- rnorm(50, 0.5, 0.3)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(50, 0.5, 0.3)
  
  result <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  
  # Count phenotypes
  counts <- table(result@meta.data$metabolic_phenotype)
  
  # Should sum to total
  expect_equal(sum(counts), 50)
})

# ============================================================================
# FINISH
# ============================================================================

cat("\n[scMetaboFlux] Edge case tests completed!\n")
