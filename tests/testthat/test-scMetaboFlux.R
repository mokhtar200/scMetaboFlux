# ============================================================================
# COMPREHENSIVE TEST SUITE FOR scMetaboFlux
# ============================================================================
# Main test file - sources helpers and runs all tests
# ============================================================================

# Source helper functions from shared file
source(file.path("helpers.R"))

# ============================================================================
# SECTION 1: BASIC FUNCTIONALITY TESTS
# ============================================================================

test_that("metabolicGeneSets is properly defined", {
  expect_type(metabolicGeneSets, "list")
  expect_true(length(metabolicGeneSets) > 0)
  expect_true("glycolysis" %in% names(metabolicGeneSets))
  expect_true("oxidative_phosphorylation" %in% names(metabolicGeneSets))
  expect_true("tca_cycle" %in% names(metabolicGeneSets))
  expect_true("fatty_acid_oxidation" %in% names(metabolicGeneSets))
})

test_that("metabolicGeneSets contain valid genes", {
  expect_true(length(metabolicGeneSets$glycolysis) > 10)
  expect_true("HK2" %in% metabolicGeneSets$glycolysis)
  expect_true("LDHA" %in% metabolicGeneSets$glycolysis)
})

test_that("rateLimitingWeights is properly defined", {
  expect_type(rateLimitingWeights, "double")
  expect_true("HK2" %in% names(rateLimitingWeights))
  expect_true("SDHA" %in% names(rateLimitingWeights))
})

test_that("atpYieldCoefficients is properly defined", {
  expect_type(atpYieldCoefficients, "double")
  expect_equal(atpYieldCoefficients["glycolysis"], 2)
  expect_equal(atpYieldCoefficients["oxidative_phosphorylation"], 32)
})

# ============================================================================
# SECTION 2: GENE MAPPING TESTS
# ============================================================================

test_that("createMetabolicGeneSetCollection works", {
  gsc <- createMetabolicGeneSetCollection()
  expect_s4_class(gsc, "GeneSetCollection")
  expect_equal(length(gsc), length(metabolicGeneSets))
})

test_that("mapGenesToPathways works", {
  genes <- c("HK2", "LDHA", "SDHA", "GENE_NOT_METABOLIC")
  mapped <- mapGenesToPathways(genes)
  
  expect_type(mapped, "list")
  expect_true("mapped_genes" %in% names(mapped))
  expect_true("glycolysis" %in% names(mapped$mapped_genes))
  expect_true("oxidative_phosphorylation" %in% names(mapped$mapped_genes))
})

test_that("mapGenesToPathways coverage is calculated", {
  genes <- c("HK2", "LDHA", "SDHA", "BRCA1", "GAPDH")
  mapped <- mapGenesToPathways(genes)
  
  expect_true("coverage" %in% names(mapped))
  expect_true(mapped$coverage$total_mapped > 0)
})

test_that("getMetabolicGeneSetStats works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  stats <- getMetabolicGeneSetStats(seurat_obj, pathways_to_include = c("glycolysis"))
  
  expect_s3_class(stats, "data.frame")
  expect_true("pathway" %in% colnames(stats))
})

# ============================================================================
# SECTION 3: PATHWAY SCORING TESTS
# ============================================================================

test_that("computePathwayScores works with mean method", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  result <- computePathwayScores(
    seurat_obj,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean",
    verbose = FALSE
  )
  
  expect_true("glycolysis" %in% colnames(result@meta.data))
  expect_true("oxidative_phosphorylation" %in% colnames(result@meta.data))
  expect_true(is.numeric(result@meta.data$glycolysis))
})

test_that("computePathwayScores adds scores correctly", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  result <- computePathwayScores(
    seurat_obj,
    method = "mean",
    verbose = FALSE
  )
  
  expect_true(length(grep("glycolysis", colnames(result@meta.data))) > 0)
})

test_that("All scoring methods produce valid results", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 150)
  
  methods <- c("mean", "weighted_mean", "single_sample")
  
  for (m in methods) {
    result <- computePathwayScores(
      seurat_obj,
      gene_sets = metabolicGeneSets[1:2],
      method = m,
      verbose = FALSE
    )
    
    expect_true("glycolysis" %in% colnames(result@meta.data))
    expect_true(all(!is.na(result@meta.data$glycolysis)))
  }
})

# ============================================================================
# SECTION 4: ATP ESTIMATION TESTS
# ============================================================================

test_that("estimateATPProduction works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(100, 100)
  
  seurat_obj@meta.data$glycolysis <- runif(100, 0, 1)
  seurat_obj@meta.data$oxidative_phosphorylation <- runif(100, 0, 1)
  
  result <- estimateATPProduction(seurat_obj, normalize = TRUE, verbose = FALSE)
  
  expect_true("ATP_score" %in% colnames(result@meta.data))
  expect_true(all(result@meta.data$ATP_score >= 0, na.rm = TRUE))
  expect_true(all(result@meta.data$ATP_score <= 1, na.rm = TRUE))
})

test_that("ATP calculation accuracy", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(100, 100)
  
  seurat_obj@meta.data$glycolysis <- rep(1.0, 100)
  seurat_obj@meta.data$oxidative_phosphorylation <- rep(1.0, 100)
  
  result <- estimateATPProduction(seurat_obj, normalize = FALSE, method = "flux_based", verbose = FALSE)
  
  expect_equal(var(result@meta.data$ATP_score), 0, tolerance = 1e-10)
})

test_that("calculateGOXIndex works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  seurat_obj@meta.data$glycolysis <- rnorm(50, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(50, mean = 0.5, sd = 0.2)
  
  result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  
  expect_true("GOX_index" %in% colnames(result@meta.data))
  expect_true("GOX_index_category" %in% colnames(result@meta.data))
})

test_that("GOX index sign correctness", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  seurat_obj@meta.data$glycolysis <- rep(1.0, 50)
  seurat_obj@meta.data$oxidative_phosphorylation <- rep(0.0, 50)
  
  result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  expect_true(all(result@meta.data$GOX_index > 0))
  
  seurat_obj@meta.data$glycolysis <- rep(0.0, 50)
  seurat_obj@meta.data$oxidative_phosphorylation <- rep(1.0, 50)
  
  result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  expect_true(all(result@meta.data$GOX_index < 0))
})

test_that("calculateATPYieldRatio works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  result <- calculateATPYieldRatio(seurat_obj, verbose = FALSE)
  
  expect_true("ATP_yield_ratio" %in% colnames(result@meta.data))
})

test_that("bootstrapATPConfidenceIntervals works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  seurat_obj@meta.data$ATP_score <- rnorm(50, mean = 0.5, sd = 0.2)
  
  ci <- bootstrapATPConfidenceIntervals(seurat_obj, n_boot = 50)
  
  expect_s3_class(ci, "data.frame")
  expect_true("mean" %in% ci$metric)
  expect_true("CI_lower" %in% ci$metric)
})

# ============================================================================
# SECTION 5: PHENOTYPE CLASSIFICATION TESTS
# ============================================================================

test_that("classifyMetabolicPhenotype works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(100, 100)
  
  seurat_obj@meta.data$glycolysis <- rnorm(100, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(100, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$ATP_score <- runif(100, 0.3, 0.8)
  
  result <- classifyMetabolicPhenotype(seurat_obj, method = "quantile", verbose = FALSE)
  
  expect_true("metabolic_phenotype" %in% colnames(result@meta.data))
  expect_s3_class(result@meta.data$metabolic_phenotype, "factor")
})

test_that("classifyMetabolicPhenotype has valid levels", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(100, 100)
  seurat_obj@meta.data$glycolysis <- rnorm(100, mean = 0.5)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(100, mean = 0.5)
  seurat_obj@meta.data$ATP_score <- runif(100, 0.3, 0.8)
  
  result <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  
  valid_levels <- c("Glycolytic", "Oxidative", "Energetically Balanced",
                   "Energy-Stressed", "Hypermetabolic")
  expect_true(all(levels(result@meta.data$metabolic_phenotype) %in% valid_levels))
})

test_that("getPhenotypeStatistics works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100, include_phenotype = TRUE)
  
  stats <- getPhenotypeStatistics(seurat_obj)
  
  expect_type(stats, "list")
  expect_true("cell_counts" %in% names(stats))
})

test_that("getPhenotypeDistribution works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100, include_phenotype = TRUE)
  
  dist <- getPhenotypeDistribution(seurat_obj)
  
  expect_s3_class(dist, "data.frame")
  expect_true("phenotype" %in% colnames(dist))
})

# ============================================================================
# SECTION 6: FLUX APPROXIMATION TESTS
# ============================================================================

test_that("computeMetabolicFlux works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  
  result <- computeMetabolicFlux(
    seurat_obj,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    use_rate_limiting = TRUE,
    normalize = TRUE,
    verbose = FALSE
  )
  
  expect_true("glycolysis_flux" %in% colnames(result@meta.data))
  expect_true("oxidative_phosphorylation_flux" %in% colnames(result@meta.data))
})

test_that("calculateFluxContribution works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  seurat_obj@meta.data$glycolysis_flux <- runif(50, 0.1, 1)
  seurat_obj@meta.data$oxidative_phosphorylation_flux <- runif(50, 0.1, 1)
  
  result <- calculateFluxContribution(seurat_obj, verbose = FALSE)
  
  contrib_cols <- grep("^contribution_", colnames(result@meta.data), value = TRUE)
  expect_true(length(contrib_cols) > 0)
})

test_that("calculateFluxCoupling works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  seurat_obj@meta.data$glycolysis_flux <- rnorm(50)
  seurat_obj@meta.data$oxidative_phosphorylation_flux <- rnorm(50)
  
  coupling <- calculateFluxCoupling(seurat_obj)
  
  expect_true(is.matrix(coupling))
})

# ============================================================================
# SECTION 7: CELL TYPE AGGREGATION TESTS
# ============================================================================

test_that("aggregateByCellType works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(100, 100)
  
  seurat_obj@meta.data$glycolysis <- runif(100, 0, 1)
  seurat_obj@meta.data$ATP_score <- runif(100, 0, 1)
  
  result <- aggregateByCellType(
    seurat_obj,
    cell_type_col = "cell_type",
    score_cols = c("glycolysis", "ATP_score")
  )
  
  expect_s3_class(result, "data.frame")
  expect_true("cell_type" %in% colnames(result))
})

test_that("calculateCellTypeSimilarity works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  seurat_obj@meta.data$glycolysis <- runif(50)
  seurat_obj@meta.data$ATP_score <- runif(50)
  
  similarity <- calculateCellTypeSimilarity(seurat_obj, cell_type_col = "cell_type")
  
  expect_true(is.matrix(similarity))
})

# ============================================================================
# SECTION 8: DIFFERENTIAL ANALYSIS TESTS
# ============================================================================

test_that("differentialMetabolicAnalysis works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(60, 100)
  seurat_obj@meta.data$glycolysis <- rnorm(60, mean = 0.5)
  seurat_obj@meta.data$ATP_score <- rnorm(60, mean = 0.5)
  
  result <- differentialMetabolicAnalysis(
    seurat_obj,
    condition_col = "condition",
    control_group = "Control",
    case_group = "Disease",
    method = "t.test"
  )
  
  expect_type(result, "list")
  expect_true("results" %in% names(result))
  expect_s3_class(result$results, "data.frame")
})

test_that("pairwiseMetabolicComparison works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(90, 100)
  seurat_obj@meta.data$glycolysis <- rnorm(90)
  
  result <- pairwiseMetabolicComparison(
    seurat_obj,
    group_col = "condition",
    method = "t.test"
  )
  
  expect_s3_class(result, "data.frame")
})

test_that("calculateEffectSizes works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(60, 100)
  
  result <- calculateEffectSizes(
    seurat_obj,
    group1_cells = seurat_obj@meta.data$condition == "Control",
    group2_cells = seurat_obj@meta.data$condition == "Disease"
  )
  
  expect_s3_class(result, "data.frame")
})

# ============================================================================
# SECTION 9: VALIDATION TESTS
# ============================================================================

test_that("validateInputData validates Seurat object", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(20, 50)
  
  validation <- validateInputData(seurat_obj, check_genes = FALSE)
  
  expect_true(validation$valid)
})

test_that("validateInputData catches invalid input", {
  result <- validateInputData("not a seurat object")
  expect_false(result$valid)
})

test_that("hasMetabolicScores returns correct values", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100)
  expect_false(hasMetabolicScores(seurat_obj))
  
  seurat_obj@meta.data$glycolysis <- runif(50)
  seurat_obj@meta.data$oxidative_phosphorylation <- runif(50)
  seurat_obj@meta.data$ATP_score <- runif(50)
  
  expect_true(hasMetabolicScores(seurat_obj))
})

# ============================================================================
# SECTION 10: WORKFLOW TESTS
# ============================================================================

test_that("runMetabolicAnalysis wrapper works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 150, include_phenotype = FALSE)
  
  result <- runMetabolicAnalysis(
    seurat_obj,
    pathways_to_analyze = c("glycolysis", "tca_cycle"),
    scoring_method = "mean",
    classify_phenotype = TRUE,
    cell_type_col = "cell_type",
    verbose = FALSE
  )
  
  expect_true("glycolysis" %in% colnames(result@meta.data))
  expect_true("ATP_score" %in% colnames(result@meta.data))
  expect_true("metabolic_phenotype" %in% colnames(result@meta.data))
})

test_that("Complete workflow integration", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 200)
  
  seurat_obj <- computePathwayScores(seurat_obj, gene_sets = metabolicGeneSets[1:4],
                                   method = "mean", verbose = FALSE)
  expect_true("glycolysis" %in% colnames(seurat_obj@meta.data))
  
  seurat_obj <- estimateATPProduction(seurat_obj, verbose = FALSE)
  expect_true("ATP_score" %in% colnames(seurat_obj@meta.data))
  
  seurat_obj <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  expect_true("GOX_index" %in% colnames(seurat_obj@meta.data))
  
  seurat_obj <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  expect_true("metabolic_phenotype" %in% colnames(seurat_obj@meta.data))
})

# ============================================================================
# SECTION 11: VISUALIZATION TESTS
# ============================================================================

test_that("plotMetabolicUMAP runs without error", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100, include_phenotype = FALSE)
  seurat_obj@meta.data$glycolysis <- runif(30)
  seurat_obj@meta.data$ATP_score <- runif(30)
  
  expect_error(
    p <- plotMetabolicUMAP(seurat_obj, combine = FALSE),
    NA
  )
  expect_true(inherits(p[[1]], "ggplot"))
})

test_that("plotPhenotypeDistribution runs without error", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100, include_phenotype = TRUE)
  
  expect_error(p <- plotPhenotypeDistribution(seurat_obj), NA)
  expect_true(inherits(p, "ggplot"))
})

test_that("plotATPDistribution runs without error", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100, include_phenotype = FALSE)
  seurat_obj@meta.data$ATP_score <- runif(30)
  
  expect_error(p <- plotATPDistribution(seurat_obj, type = "violin"), NA)
  expect_true(inherits(p, "ggplot"))
})

test_that("plotPathwayCorrelation runs without error", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100, include_phenotype = FALSE)
  seurat_obj@meta.data$glycolysis <- runif(30)
  seurat_obj@meta.data$ATP_score <- runif(30)
  
  expect_error(p <- plotPathwayCorrelation(seurat_obj), NA)
  expect_true(inherits(p, "ggplot"))
})

# ============================================================================
# SECTION 12: EXPORT AND S4 CLASS TESTS
# ============================================================================

test_that("resetMetabolicResults removes scores", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100)
  seurat_obj@meta.data$glycolysis <- runif(30)
  seurat_obj@meta.data$ATP_score <- runif(30)
  
  result <- resetMetabolicResults(seurat_obj)
  
  expect_false("glycolysis" %in% colnames(result@meta.data))
  expect_false("ATP_score" %in% colnames(result@meta.data))
})

test_that("scMetaboFlux object creation works", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(30, 100)
  
  smf_obj <- createScMetaboFlux(seurat_obj, name = "test")
  
  expect_s4_class(smf_obj, "scMetaboFluxObject")
  expect_equal(smf_obj@name, "test")
})

# ============================================================================
# SECTION 13: EDGE CASE TESTS
# ============================================================================

test_that("NA and NaN values in metadata handled", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100, include_phenotype = FALSE)
  
  seurat_obj@meta.data$glycolysis[1:5] <- NA
  seurat_obj@meta.data$ATP_score[10:15] <- NaN
  
  result <- estimateATPProduction(seurat_obj, verbose = FALSE)
  
  expect_true("ATP_score" %in% colnames(result@meta.data))
  expect_true(all(!is.nan(result@meta.data$ATP_score)))
})

test_that("Missing pathway columns handled", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100, include_phenotype = FALSE)
  
  expect_error(
    estimateATPProduction(seurat_obj, glycolysis_col = "nonexistent"),
    "not found"
  )
})

test_that("Zero variance scores handled", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100, include_phenotype = FALSE)
  
  seurat_obj@meta.data$glycolysis <- rep(0.5, 50)
  seurat_obj@meta.data$oxidative_phosphorylation <- rep(0.5, 50)
  
  result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  expect_true("GOX_index" %in% colnames(result@meta.data))
})

test_that("Single cell type handled", {
  skip_if_missing()
  
  seurat_obj <- create_test_seurat(50, 100, include_phenotype = FALSE)
  seurat_obj@meta.data$cell_type <- rep("SingleType", 50)
  seurat_obj@meta.data$glycolysis <- rnorm(50)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(50)
  
  expect_no_error(
    result <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type")
  )
})

# ============================================================================
# FINISH
# ============================================================================

cat("\n[scMetaboFlux] All core tests completed successfully!\n")
