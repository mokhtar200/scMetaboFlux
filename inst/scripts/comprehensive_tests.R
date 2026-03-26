#!/usr/bin/env Rscript
# ============================================================================
# COMPREHENSIVE TEST SCRIPT FOR scMetaboFlux
# Tests the package with large, complex simulated datasets
# ============================================================================

cat("========================================\n")
cat("scMetaboFlux Comprehensive Test Suite\n")
cat("========================================\n\n")

# Set up test environment
suppressPackageStartupMessages({
  library(Seurat)
  library(scMetaboFlux)
})

# Helper function to print test results
print_test <- function(test_name, result, details = NULL) {
  status <- if (result) "PASS" else "FAIL"
  cat(sprintf("[%s] %s\n", status, test_name))
  if (!is.null(details) && !result) {
    cat(sprintf("       Details: %s\n", details))
  }
  return(result)
}

# Test 1: Generate Example Data - PBMC-like
cat("\n--- Test 1: Generate PBMC-like Dataset ---\n")
tryCatch({
  set.seed(42)
  pbmc_data <- generateExampleData(n_cells = 500, n_genes = 1000)
  has_data <- inherits(pbmc_data, "Seurat")
  print_test("Generate PBMC-like dataset", has_data)
  
  if (has_data) {
    n_cells <- ncol(pbmc_data)
    n_genes <- nrow(pbmc_data)
    print_test(sprintf("Dataset has %d cells", n_cells), n_cells >= 400)
    print_test(sprintf("Dataset has %d genes", n_genes), n_genes >= 800)
    print_test("Dataset has cell type annotations", "cell_type" %in% colnames(pbmc_data@meta.data))
  }
}, error = function(e) {
  print_test("Generate PBMC-like dataset", FALSE, e$message)
})

# Test 2: Generate Cancer-like Dataset
cat("\n--- Test 2: Generate Cancer-like Dataset ---\n")
tryCatch({
  set.seed(123)
  cancer_data <- generateCancerData(n_cells = 1000, n_genes = 2000)
  has_data <- inherits(cancer_data, "Seurat")
  print_test("Generate cancer-like dataset", has_data)
  
  if (has_data) {
    n_cells <- ncol(cancer_data)
    print_test(sprintf("Large cancer dataset: %d cells", n_cells), n_cells >= 800)
  }
}, error = function(e) {
  print_test("Generate cancer-like dataset", FALSE, e$message)
})

# Test 3: Metabolic Gene Sets
cat("\n--- Test 3: Metabolic Gene Sets ---\n")
tryCatch({
  gene_sets <- metabolicGeneSets
  print_test("metabolicGeneSets available", is.list(gene_sets))
  print_test("Has glycolysis pathway", "glycolysis" %in% names(gene_sets))
  print_test("Has oxidative_phosphorylation pathway", "oxidative_phosphorylation" %in% names(gene_sets))
  print_test("Has tca_cycle pathway", "tca_cycle" %in% names(gene_sets))
  print_test("Has fatty_acid_oxidation pathway", "fatty_acid_oxidation" %in% names(gene_sets))
  
  n_pathways <- length(gene_sets)
  print_test(sprintf("Has %d metabolic pathways", n_pathways), n_pathways >= 10)
  
  glycolysis_genes <- gene_sets$glycolysis
  print_test("Glycolysis has genes", length(glycolysis_genes) > 10)
  
  # Rate limiting weights
  rl_weights <- rateLimitingWeights
  print_test("rateLimitingWeights available", is.list(rl_weights))
  
  # ATP yield coefficients
  atp_coef <- atpYieldCoefficients
  print_test("atpYieldCoefficients available", is.numeric(atp_coef))
}, error = function(e) {
  print_test("Metabolic gene sets", FALSE, e$message)
})

# Test 4: Create Gene Set Collection
cat("\n--- Test 4: Create Gene Set Collection ---\n")
tryCatch({
  gsc <- createMetabolicGeneSetCollection(gene_sets = metabolicGeneSets[1:3])
  print_test("Create gene set collection", inherits(gsc, "GeneSetCollection"))
}, error = function(e) {
  print_test("Create gene set collection", FALSE, e$message)
})

# Test 5: Custom Gene Set
cat("\n--- Test 5: Custom Gene Set ---\n")
tryCatch({
  custom_gs <- createCustomGeneSet(
    genes = c("HK2", "PKM2", "LDHA", "GAPDH"),
    name = "test_pathway",
    description = "Test custom pathway"
  )
  print_test("Create custom gene set", inherits(custom_gs, "GeneSet"))
}, error = function(e) {
  print_test("Create custom gene set", FALSE, e$message)
})

# Test 6: Validate Input Data
cat("\n--- Test 6: Validate Input Data ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500)
  validation <- validateInputData(pbmc_data, check_genes = FALSE)
  print_test("Validate input data", is.logical(validation))
}, error = function(e) {
  print_test("Validate input data", FALSE, e$message)
})

# Test 7: Run Complete Metabolic Analysis
cat("\n--- Test 7: Run Complete Metabolic Analysis ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 300, n_genes = 500, seed = 42)
  
  result <- runMetabolicAnalysis(
    pbmc_data,
    pathways_to_analyze = c("glycolysis", "oxidative_phosphorylation", "tca_cycle"),
    scoring_method = "mean",
    classify_phenotype = FALSE
  )
  
  print_test("Run complete metabolic analysis", inherits(result, "Seurat"))
  
  if (inherits(result, "Seurat")) {
    score_cols <- grep("^(glycolysis|oxidative|ATP)", 
                      colnames(result@meta.data), value = TRUE)
    print_test(sprintf("Created %d metabolic score columns", length(score_cols)), 
               length(score_cols) >= 2)
  }
}, error = function(e) {
  print_test("Run complete metabolic analysis", FALSE, e$message)
})

# Test 8: Pathway Scoring Methods
cat("\n--- Test 8: Pathway Scoring Methods ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # Test AUCell
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "AUCell"
  )
  has_aucell <- "glycolysis_AUCell" %in% colnames(pbmc_data@meta.data)
  print_test("AUCell pathway scoring", has_aucell)
  
  # Test GSVA
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("tca_cycle")],
    method = "GSVA"
  )
  has_gsva <- "tca_cycle_GSVA" %in% colnames(pbmc_data@meta.data)
  print_test("GSVA pathway scoring", has_gsva)
  
  # Test mean
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("fatty_acid_oxidation")],
    method = "mean"
  )
  has_mean <- "fatty_acid_oxidation_mean" %in% colnames(pbmc_data@meta.data)
  print_test("Mean pathway scoring", has_mean)
  
}, error = function(e) {
  print_test("Pathway scoring methods", FALSE, e$message)
})

# Test 9: ATP Estimation
cat("\n--- Test 9: ATP Estimation ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # First compute pathway scores
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean"
  )
  
  # Estimate ATP
  pbmc_data <- estimateATPProduction(pbmc_data, glycolysis_weight = 0.3, oxphos_weight = 0.7)
  has_atp <- "ATP_score" %in% colnames(pbmc_data@meta.data)
  print_test("Estimate ATP production", has_atp)
  
  if (has_atp) {
    atp_values <- pbmc_data@meta.data$ATP_score
    print_test("ATP values are numeric", is.numeric(atp_values))
    print_test("ATP values are non-negative", all(atp_values >= 0, na.rm = TRUE))
  }
  
  # Calculate GOX Index
  pbmc_data <- calculateGOXIndex(pbmc_data)
  has_gox <- "GOX_index" %in% colnames(pbmc_data@meta.data)
  print_test("Calculate GOX index", has_gox)
  
  # ATP Yield Ratio
  pbmc_data <- calculateATPYieldRatio(pbmc_data)
  has_ratio <- "ATP_yield_ratio" %in% colnames(pbmc_data@meta.data)
  print_test("Calculate ATP yield ratio", has_ratio)
  
}, error = function(e) {
  print_test("ATP estimation", FALSE, e$message)
})

# Test 10: Flux Approximation
cat("\n--- Test 10: Flux Approximation ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # Compute pathway scores first
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation", "tca_cycle")],
    method = "mean"
  )
  
  # Compute flux
  pbmc_data <- computeMetabolicFlux(
    pbmc_data,
    flux_method = "pseudo_flux",
    use_rate_limiting = TRUE
  )
  
  flux_cols <- grep("flux", colnames(pbmc_data@meta.data), value = TRUE)
  print_test("Compute metabolic flux", length(flux_cols) > 0)
  
  # Flux coupling
  coupling <- calculateFluxCoupling(pbmc_data)
  print_test("Calculate flux coupling", is.list(coupling) || is.matrix(coupling))
  
  # Flux balance
  balance <- calculateFluxBalance(pbmc_data)
  print_test("Calculate flux balance", is.numeric(balance) || is.list(balance))
  
}, error = function(e) {
  print_test("Flux approximation", FALSE, e$message)
})

# Test 11: Phenotype Classification
cat("\n--- Test 11: Phenotype Classification ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # Compute scores first
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean"
  )
  pbmc_data <- estimateATPProduction(pbmc_data)
  
  # Quantile classification
  pbmc_data <- classifyMetabolicPhenotype(pbmc_data, method = "quantile")
  has_phenotype <- "metabolic_phenotype" %in% colnames(pbmc_data@meta.data)
  print_test("Classify by quantile", has_phenotype)
  
  if (has_phenotype) {
    phenotypes <- unique(pbmc_data@meta.data$metabolic_phenotype)
    print_test(sprintf("Identified %d phenotypes", length(phenotypes)), length(phenotypes) >= 2)
  }
  
  # Get phenotype distribution
  dist <- getPhenotypeDistribution(pbmc_data)
  print_test("Get phenotype distribution", is.table(dist) || is.data.frame(dist))
  
  # Get phenotype statistics
  stats <- getPhenotypeStatistics(pbmc_data)
  print_test("Get phenotype statistics", is.data.frame(stats))
  
}, error = function(e) {
  print_test("Phenotype classification", FALSE, e$message)
})

# Test 12: Cell-Type Aggregation
cat("\n--- Test 12: Cell-Type Aggregation ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 300, n_genes = 500, seed = 42)
  
  # Compute scores
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean"
  )
  
  # Aggregate by cell type
  celltype_summary <- aggregateByCellType(pbmc_data, celltype_col = "cell_type")
  print_test("Aggregate by cell type", is.data.frame(celltype_summary))
  
  if (is.data.frame(celltype_summary)) {
    print_test("Summary has cell type column", "cell_type" %in% colnames(celltype_summary))
    print_test(sprintf("Summary has %d cell types", nrow(celltype_summary)), 
               nrow(celltype_summary) >= 3)
  }
  
  # Compare cell types
  comparison <- compareCellTypesMetabolism(pbmc_data, celltype_col = "cell_type")
  print_test("Compare cell types", is.list(comparison) || is.data.frame(comparison))
  
  # Rank cell types
  ranked <- rankCellTypesByMetabolism(pbmc_data, celltype_col = "cell_type")
  print_test("Rank cell types", is.data.frame(ranked) || is.character(ranked))
  
}, error = function(e) {
  print_test("Cell-type aggregation", FALSE, e$message)
})

# Test 13: Differential Analysis
cat("\n--- Test 13: Differential Analysis ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 300, n_genes = 500, seed = 42)
  
  # Compute scores
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean"
  )
  
  # Get cell types
  cell_types <- unique(pbmc_data@meta.data$cell_type)
  
  if (length(cell_types) >= 2) {
    # Differential analysis
    diff_result <- differentialMetabolicAnalysis(
      pbmc_data,
      group_col = "cell_type",
      group1 = cell_types[1],
      group2 = cell_types[2]
    )
    print_test("Differential metabolic analysis", is.data.frame(diff_result) || is.list(diff_result))
    
    # Effect sizes
    effects <- calculateEffectSizes(pbmc_data, group_col = "cell_type")
    print_test("Calculate effect sizes", is.data.frame(effects) || is.list(effects))
  }
  
}, error = function(e) {
  print_test("Differential analysis", FALSE, e$message)
})

# Test 14: scMetaboFlux S4 Object
cat("\n--- Test 14: scMetaboFlux S4 Object ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # Create scMetaboFlux object
  smf_obj <- createScMetaboFlux(pbmc_data, name = "test_analysis")
  print_test("Create scMetaboFlux object", inherits(smf_obj, "scMetaboFluxObject"))
  
  # Test summary method
  summary_result <- summary(smf_obj)
  print_test("Summary method", TRUE)  # summary prints to console
  
  # Test print method
  print_result <- print(smf_obj)
  print_test("Print method", TRUE)
  
  # Extract results
  all_results <- extractResults(smf_obj, what = "all")
  print_test("Extract all results", is.data.frame(all_results))
  
}, error = function(e) {
  print_test("scMetaboFlux S4 object", FALSE, e$message)
})

# Test 15: Z-score Normalization
cat("\n--- Test 15: Z-score Normalization ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # Compute scores
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation")],
    method = "mean"
  )
  
  # Z-score
  pbmc_data <- zscorePathwayScores(pbmc_data)
  zscore_cols <- grep("_zscore", colnames(pbmc_data@meta.data), value = TRUE)
  print_test("Z-score normalization", length(zscore_cols) > 0)
  
}, error = function(e) {
  print_test("Z-score normalization", FALSE, e$message)
})

# Test 16: Combined Score
cat("\n--- Test 16: Combined Score ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 200, n_genes = 500, seed = 42)
  
  # Compute scores
  pbmc_data <- computePathwayScores(
    pbmc_data,
    gene_sets = metabolicGeneSets[c("glycolysis", "oxidative_phosphorylation", "tca_cycle")],
    method = "mean"
  )
  
  # Combined score
  pbmc_data <- computeCombinedScore(
    pbmc_data,
    score_cols = c("glycolysis_mean", "oxidative_phosphorylation_mean", "tca_cycle_mean"),
    weights = c(0.4, 0.4, 0.2)
  )
  
  has_combined <- "combined_score" %in% colnames(pbmc_data@meta.data)
  print_test("Compute combined score", has_combined)
  
}, error = function(e) {
  print_test("Combined score", FALSE, e$message)
})

# Test 17: Large Scale Dataset
cat("\n--- Test 17: Large Scale Dataset (1000 cells) ---\n")
tryCatch({
  large_data <- generateExampleData(n_cells = 1000, n_genes = 3000, seed = 42)
  print_test("Generate large dataset", inherits(large_data, "Seurat"))
  
  if (inherits(large_data, "Seurat")) {
    n_cells <- ncol(large_data)
    n_genes <- nrow(large_data)
    print_test(sprintf("Large dataset: %d cells x %d genes", n_cells, n_genes), 
               n_cells >= 800 && n_genes >= 2000)
    
    # Run full analysis on large dataset
    large_data <- runMetabolicAnalysis(
      large_data,
      pathways_to_analyze = c("glycolysis", "oxidative_phosphorylation", "tca_cycle", 
                             "fatty_acid_oxidation"),
      scoring_method = "mean",
      classify_phenotype = TRUE
    )
    
    score_cols <- grep("^(glycolysis|oxidative|ATP|metabolic)", 
                      colnames(large_data@meta.data), value = TRUE)
    print_test(sprintf("Large dataset analysis: %d scores", length(score_cols)), 
               length(score_cols) >= 5)
  }
}, error = function(e) {
  print_test("Large scale dataset", FALSE, e$message)
})

# Test 18: Hypoxia Response Analysis
cat("\n--- Test 18: Hypoxia Response Analysis ---\n")
tryCatch({
  pbmc_data <- generateExampleData(n_cells = 300, n_genes = 500, seed = 42)
  
  # Check if hypoxia pathway exists
  has_hypoxia <- "hypoxia_response" %in% names(metabolicGeneSets)
  print_test("Hypoxia response pathway available", has_hypoxia)
  
  if (has_hypoxia) {
    # Score hypoxia pathway
    pbmc_data <- computePathwayScores(
      pbmc_data,
      gene_sets = metabolicGeneSets["hypoxia_response"],
      method = "mean"
    )
    has_hypoxia_score <- "hypoxia_response_mean" %in% colnames(pbmc_data@meta.data)
    print_test("Score hypoxia pathway", has_hypoxia_score)
  }
  
}, error = function(e) {
  print_test("Hypoxia response analysis", FALSE, e$message)
})

# Test 19: Map Genes to Pathways
cat("\n--- Test 19: Map Genes to Pathways ---\n")
tryCatch({
  test_genes <- c("HK2", "PKM2", "LDHA", "NDUFA9", "SDHA", "COX4I1")
  mapping <- mapGenesToPathways(test_genes)
  print_test("Map genes to pathways", is.data.frame(mapping) || is.list(mapping))
  
}, error = function(e) {
  print_test("Map genes to pathways", FALSE, e$message)
})

# Test 20: Validate Package Integrity
cat("\n--- Test 20: Package Integrity ---\n")
tryCatch({
  # Check all exported functions exist
  exported_funcs <- c(
    "metabolicGeneSets", "rateLimitingWeights", "atpYieldCoefficients",
    "generateExampleData", "runMetabolicAnalysis", "computePathwayScores",
    "estimateATPProduction", "computeMetabolicFlux", "classifyMetabolicPhenotype",
    "aggregateByCellType", "createScMetaboFlux"
  )
  
  all_exist <- all(sapply(exported_funcs, exists))
  print_test("All core functions exported", all_exist)
  
  # Check package version
  pkg_info <- packageDescription("scMetaboFlux")
  has_version <- !is.null(pkg_info$Version)
  print_test("Package has version", has_version)
  if (has_version) {
    cat(sprintf("       Version: %s\n", pkg_info$Version))
  }
  
}, error = function(e) {
  print_test("Package integrity", FALSE, e$message)
})

# Summary
cat("\n========================================\n")
cat("Test Suite Complete\n")
cat("========================================\n")

# Count passed/failed
cat("\nNote: Check output above for individual test results.\n")
