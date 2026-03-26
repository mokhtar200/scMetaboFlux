# ============================================================================
# INTEGRATION TESTS FOR scMetaboFlux COMPLETE WORKFLOWS
# ============================================================================
# Tests for complete analysis pipelines and real-world scenarios
# ============================================================================

library(testthat)
library(Seurat)

# Helper function
.create_integration_seurat <- function(n_cells = 100, seed = 42) {
  set.seed(seed)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1")
  tca_genes <- c("CS", "IDH1", "IDH2", "OGDH", "FH")
  other_genes <- paste0("GENE", 1:50)
  
  all_genes <- c(glycolytic_genes, oxphos_genes, tca_genes, other_genes)
  
  expr_matrix <- matrix(
    rnbinom(length(all_genes) * n_cells, size = 2, mu = 5),
    nrow = length(all_genes),
    ncol = n_cells
  )
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  
  seurat_obj@meta.data$cell_type <- sample(
    c("T_cell", "B_cell", "Macrophage"),
    n_cells, replace = TRUE
  )
  
  seurat_obj@meta.data$condition <- sample(
    c("Control", "Disease"),
    n_cells, replace = TRUE
  )
  
  seurat_obj@meta.data$patient_id <- sample(
    paste0("Patient_", 1:3),
    n_cells, replace = TRUE
  )
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = 10)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)
  
  return(seurat_obj)
}

# ============================================================================
# TEST 1: Basic Analysis Pipeline
# ============================================================================

test_that("Basic analysis pipeline runs end-to-end", {
  skip_if_missing()
  
  # Create test data
  seurat_obj <- .create_integration_seurat(50)
  
  # Step 1: Compute pathway scores
  seurat_obj <- computePathwayScores(
    seurat_obj,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean",
    verbose = FALSE
  )
  
  expect_true("glycolysis" %in% colnames(seurat_obj@meta.data))
  expect_true("oxidative_phosphorylation" %in% colnames(seurat_obj@meta.data))
  
  # Step 2: Estimate ATP
  seurat_obj <- estimateATPProduction(seurat_obj, verbose = FALSE)
  
  expect_true("ATP_score" %in% colnames(seurat_obj@meta.data))
  
  # Step 3: Calculate GOX index
  seurat_obj <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  
  expect_true("GOX_index" %in% colnames(seurat_obj@meta.data))
  
  # Step 4: Compute flux
  seurat_obj <- computeMetabolicFlux(
    seurat_obj,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    verbose = FALSE
  )
  
  expect_true("glycolysis_flux" %in% colnames(seurat_obj@meta.data))
  
  # Step 5: Classify phenotype
  seurat_obj <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  
  expect_true("metabolic_phenotype" %in% colnames(seurat_obj@meta.data))
  
  # Step 6: Aggregate by cell type
  stats <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type")
  
  expect_s3_class(stats, "data.frame")
  expect_true("cell_type" %in% colnames(stats))
})

# ============================================================================
# TEST 2: Complete Workflow with Wrapper Function
# ============================================================================

test_that("runMetabolicAnalysis wrapper works correctly", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(50)
  
  result <- runMetabolicAnalysis(
    seurat_obj,
    pathways_to_analyze = c("glycolysis", "oxidative_phosphorylation"),
    scoring_method = "mean",
    classify_phenotype = TRUE,
    cell_type_col = "cell_type",
    verbose = FALSE
  )
  
  # Check all expected columns exist
  expected_cols <- c(
    "glycolysis", "oxidative_phosphorylation",
    "ATP_score", "GOX_index",
    "metabolic_phenotype"
  )
  
  for (col in expected_cols) {
    expect_true(col %in% colnames(result@meta.data))
  }
  
  # Check misc slot has summary
  expect_true(!is.null(result@misc$metabolic_summary))
})

# ============================================================================
# TEST 3: Disease vs Control Comparison Workflow
# ============================================================================

test_that("Disease vs control comparison workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(80)
  
  # Add metabolic differences between conditions
  case_cells <- seurat_obj@meta.data$condition == "Disease"
  
  # Add simulated metabolic shift
  glycolysis_signal <- rnorm(80, mean = 0.5, sd = 0.2)
  glycolysis_signal[case_cells] <- glycolysis_signal[case_cells] + 0.5
  
  seurat_obj@meta.data$glycolysis <- glycolysis_signal
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(80, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$ATP_score <- runif(80, 0.3, 0.9)
  
  # Differential analysis
  diff_results <- differentialMetabolicAnalysis(
    seurat_obj,
    condition_col = "condition",
    control_group = "Control",
    case_group = "Disease",
    method = "wilcox.test"
  )
  
  expect_type(diff_results, "list")
  expect_true("results" %in% names(diff_results))
  
  gly_result <- diff_results$results[diff_results$results$metric == "glycolysis", ]
  
  # Disease should have higher glycolysis
  expect_true(gly_result$case_mean > gly_result$control_mean)
})

# ============================================================================
# TEST 4: Cell Type Comparison Workflow
# ============================================================================

test_that("Cell type comparison workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(60)
  
  # Make cell types have different metabolic profiles
  t_cells <- seurat_obj@meta.data$cell_type == "T_cell"
  b_cells <- seurat_obj@meta.data$cell_type == "B_cell"
  
  seurat_obj@meta.data$glycolysis <- rnorm(60, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$glycolysis[t_cells] <- rnorm(sum(t_cells), mean = 0.3, sd = 0.1)
  seurat_obj@meta.data$glycolysis[b_cells] <- rnorm(sum(b_cells), mean = 0.8, sd = 0.1)
  
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(60, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$ATP_score <- runif(60, 0.3, 0.9)
  
  # Compare cell types
  comparison <- compareCellTypesMetabolism(
    seurat_obj,
    cell_type_col = "cell_type",
    method = "anova"
  )
  
  expect_type(comparison, "list")
  
  # Aggregate by cell type
  agg_stats <- aggregateByCellType(
    seurat_obj,
    cell_type_col = "cell_type",
    score_cols = c("glycolysis", "ATP_score")
  )
  
  expect_equal(nrow(agg_stats), 3)  # 3 cell types
  
  # Get cell type rankings
  rankings <- rankCellTypesByMetabolism(
    seurat_obj,
    cell_type_col = "cell_type",
    rank_by = "ATP_score"
  )
  
  expect_true("rank" %in% colnames(rankings))
})

# ============================================================================
# TEST 5: Visualization Pipeline
# ============================================================================

test_that("All visualization functions work together", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(40)
  seurat_obj@meta.data$metabolic_phenotype <- sample(
    c("Glycolytic", "Oxidative", "Energetically Balanced"),
    40, replace = TRUE
  )
  
  # Test each visualization function
  expect_error(p1 <- plotMetabolicUMAP(seurat_obj, combine = FALSE), NA)
  expect_true(is_valid_ggplot(p1[[1]]))
  
  expect_error(p2 <- plotPhenotypeDistribution(seurat_obj), NA)
  expect_true(is_valid_ggplot(p2))
  
  expect_error(p3 <- plotATPDistribution(seurat_obj, type = "violin"), NA)
  expect_true(is_valid_ggplot(p3))
  
  expect_error(p4 <- plotGlycolysisVsOxphos(seurat_obj), NA)
  expect_true(is_valid_ggplot(p4))
  
  expect_error(p5 <- plotPathwayCorrelation(seurat_obj), NA)
  expect_true(is_valid_ggplot(p5))
  
  expect_error(p6 <- plotGOXIndex(seurat_obj), NA)
  expect_true(is_valid_ggplot(p6))
})

# ============================================================================
# TEST 6: Export and Import Workflow
# ============================================================================

test_that("Export and data handling workflow", {
  skip_if_missing()
  
  skip_on_cran()
  
  seurat_obj <- .create_integration_seurat(30)
  seurat_obj@meta.data$metabolic_phenotype <- sample(
    c("Glycolytic", "Oxidative"),
    30, replace = TRUE
  )
  
  temp_dir <- tempdir()
  
  # Export results
  exportMetabolicResults(
    seurat_obj,
    output_dir = temp_dir,
    include_embeddings = TRUE
  )
  
  # Check files were created
  expect_true(file.exists(file.path(temp_dir, "metabolic_scores.csv")))
  expect_true(file.exists(file.path(temp_dir, "phenotype_distribution.csv")))
  
  # Export cell type summary
  exportCellTypeSummary(
    seurat_obj,
    cell_type_col = "cell_type",
    output_dir = temp_dir
  )
  
  expect_true(file.exists(file.path(temp_dir, "celltype_metabolism_aggregated.csv")))
})

# ============================================================================
# TEST 7: Reset and Re-run Workflow
# ============================================================================

test_that("Reset and re-run workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(30)
  
  # First analysis
  seurat_obj <- computePathwayScores(
    seurat_obj,
    gene_sets = metabolicGeneSets["glycolysis"],
    method = "mean",
    verbose = FALSE
  )
  
  expect_true("glycolysis" %in% colnames(seurat_obj@meta.data))
  
  # Reset
  seurat_obj <- resetMetabolicResults(seurat_obj)
  
  expect_false("glycolysis" %in% colnames(seurat_obj@meta.data))
  
  # Re-run with different method
  seurat_obj <- computePathwayScores(
    seurat_obj,
    gene_sets = metabolicGeneSets["glycolysis"],
    method = "weighted_mean",
    verbose = FALSE
  )
  
  expect_true("glycolysis" %in% colnames(seurat_obj@meta.data))
})

# ============================================================================
# TEST 8: Object Validation Workflow
# ============================================================================

test_that("Object validation workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(20)
  
  # Validate input
  validation <- validateInputData(seurat_obj, check_genes = FALSE)
  expect_true(validation$valid)
  
  # Check scores existence
  expect_false(hasMetabolicScores(seurat_obj))
  
  # Add scores
  seurat_obj@meta.data$glycolysis <- rnorm(20)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(20)
  seurat_obj@meta.data$ATP_score <- runif(20)
  
  expect_true(hasMetabolicScores(seurat_obj))
  
  # Get parameters
  params <- getAnalysisParameters(seurat_obj)
  expect_true(is.null(params))  # No analysis run yet
})

# ============================================================================
# TEST 9: Multi-Method Comparison
# ============================================================================

test_that("Multi-method scoring comparison", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(30)
  
  # Compare methods
  methods <- c("mean", "weighted_mean")
  
  results <- lapply(methods, function(m) {
    obj <- computePathwayScores(
      seurat_obj,
      gene_sets = metabolicGeneSets["glycolysis"],
      method = m,
      verbose = FALSE
    )
    obj@meta.data$glycolysis
  })
  
  # Results should be correlated but not identical
  cor_val <- cor(results[[1]], results[[2]], use = "pairwise")
  expect_true(cor_val > 0.9)  # Highly correlated
  expect_false(identical(results[[1]], results[[2]]))  # But not identical
})

# ============================================================================
# TEST 10: Flux Analysis Workflow
# ============================================================================

test_that("Complete flux analysis workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(40)
  
  # Compute flux
  seurat_obj <- computeMetabolicFlux(
    seurat_obj,
    gene_sets = metabolicGeneSets[c("glycolysis", "tca_cycle", "oxidative_phosphorylation")],
    verbose = FALSE
  )
  
  expect_true("glycolysis_flux" %in% colnames(seurat_obj@meta.data))
  expect_true("tca_cycle_flux" %in% colnames(seurat_obj@meta.data))
  
  # Calculate coupling
  coupling <- calculateFluxCoupling(seurat_obj)
  expect_true(is.matrix(coupling))
  
  # Calculate contribution
  seurat_obj <- calculateFluxContribution(seurat_obj, verbose = FALSE)
  
  contrib_cols <- grep("^contribution_", colnames(seurat_obj@meta.data), value = TRUE)
  expect_true(length(contrib_cols) > 0)
  
  # Calculate entropy
  seurat_obj <- calculateFluxEntropy(seurat_obj, verbose = FALSE)
  expect_true("flux_entropy" %in% colnames(seurat_obj@meta.data))
  
  # Calculate total flux
  seurat_obj <- calculateTotalMetabolicFlux(seurat_obj, verbose = FALSE)
  expect_true("total_flux" %in% colnames(seurat_obj@meta.data))
})

# ============================================================================
# TEST 11: Phenotype Analysis Workflow
# ============================================================================

test_that("Complete phenotype analysis workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(50)
  
  # Classify phenotypes
  seurat_obj <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  
  expect_true("metabolic_phenotype" %in% colnames(seurat_obj@meta.data))
  
  # Get statistics
  phen_stats <- getPhenotypeStatistics(seurat_obj)
  expect_type(phen_stats, "list")
  expect_true("cell_counts" %in% names(phen_stats))
  
  # Get distribution
  dist <- getPhenotypeDistribution(seurat_obj)
  expect_s3_class(dist, "data.frame")
  
  # Calculate purity
  purity <- calculatePhenotypePurity(seurat_obj)
  expect_s3_class(purity, "data.frame")
  expect_true("purity" %in% colnames(purity))
  
  # Refine classification
  seurat_obj <- refinePhenotypeClassification(
    seurat_obj,
    k = 5,
    name = "metabolic_phenotype_refined"
  )
  expect_true("metabolic_phenotype_refined" %in% colnames(seurat_obj@meta.data))
})

# ============================================================================
# TEST 12: S4 Object Workflow
# ============================================================================

test_that("S4 object creation and manipulation workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(30)
  
  # Create scMetaboFlux object
  smf_obj <- createScMetaboFlux(seurat_obj, name = "test_workflow")
  
  expect_s4_class(smf_obj, "scMetaboFluxObject")
  expect_equal(smf_obj@name, "test_workflow")
  
  # Extract results
  scores <- extractResults(smf_obj, what = "scores")
  expect_s3_class(scores, "data.frame")
  
  # Print and summary
  expect_error(print(smf_obj), NA)
  expect_error(summary(smf_obj), NA)
})

# ============================================================================
# TEST 13: Large Dataset Workflow
# ============================================================================

test_that("Large dataset workflow (performance test)", {
  skip_if_missing()
  skip_on_cran()
  
  # Create larger dataset
  seurat_obj <- .create_integration_seurat(200)
  
  start_time <- Sys.time()
  
  result <- runMetabolicAnalysis(
    seurat_obj,
    pathways_to_analyze = c("glycolysis", "oxidative_phosphorylation", "tca_cycle"),
    scoring_method = "mean",
    classify_phenotype = TRUE,
    cell_type_col = "cell_type",
    verbose = FALSE
  )
  
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Should complete in reasonable time
  expect_true(elapsed < 120)
  
  # Results should be valid
  expect_true("ATP_score" %in% colnames(result@meta.data))
  expect_true("metabolic_phenotype" %in% colnames(result@meta.data))
})

# ============================================================================
# TEST 14: Z-Score Normalization Workflow
# ============================================================================

test_that("Z-score normalization workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(40)
  
  # Add scores
  seurat_obj@meta.data$glycolysis <- rnorm(40, mean = 5, sd = 2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(40, mean = 10, sd = 3)
  
  # Z-score
  seurat_obj <- zscorePathwayScores(seurat_obj)
  
  expect_true("zscore_glycolysis" %in% colnames(seurat_obj@meta.data))
  expect_true("zscore_oxidative_phosphorylation" %in% colnames(seurat_obj@meta.data))
  
  # Check z-score properties
  zscore_gly <- seurat_obj@meta.data$zscore_glycolysis
  expect_true(abs(mean(zscore_gly)) < 0.1)  # Mean ~ 0
  expect_true(abs(sd(zscore_gly) - 1) < 0.1)  # SD ~ 1
})

# ============================================================================
# TEST 15: Cross-Validation Simulation
# ============================================================================

test_that("Cross-validation simulation workflow", {
  skip_if_missing()
  
  seurat_obj <- .create_integration_seurat(60)
  
  # Create train/test split
  train_idx <- sample(1:60, 40, replace = FALSE)
  test_idx <- setdiff(1:60, train_idx)
  
  seurat_obj@meta.data$glycolysis <- rnorm(60, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(60, mean = 0.5, sd = 0.2)
  
  # Calculate ATP on full data
  seurat_obj <- estimateATPProduction(seurat_obj, verbose = FALSE)
  
  # Simulate flux
  seurat_obj <- computeMetabolicFlux(
    seurat_obj,
    gene_sets = list(glycolysis = c("HK2", "LDHA")),
    verbose = FALSE
  )
  
  # Bootstrap confidence intervals
  ci <- bootstrapATPConfidenceIntervals(
    seurat_obj,
    n_boot = 50
  )
  
  expect_s3_class(ci, "data.frame")
  expect_true("mean" %in% ci$metric)
  expect_true("CI_lower" %in% ci$metric)
  expect_true("CI_upper" %in% ci$metric)
})

# ============================================================================
# FINISH
# ============================================================================

cat("\n[scMetaboFlux] Integration tests completed!\n")
