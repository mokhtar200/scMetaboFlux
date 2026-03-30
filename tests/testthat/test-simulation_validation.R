# ============================================================================
# SIMULATION-BASED BIOLOGICAL VALIDATION TESTS FOR scMetaboFlux
# ============================================================================
# This file contains simulation-based tests that validate biological
# assumptions and expected behaviors of the metabolic analysis pipeline
# ============================================================================

library(testthat)
library(Seurat)

# Helper function to create test data without circular dependency
.create_test_seurat_for_sim <- function(n_cells = 100, seed = 42) {
  set.seed(seed)
  
  # Create expression matrix with metabolic genes
  glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH", "TPI1", "PFKM")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1", "NDUFS1", "SDHB")
  tca_genes <- c("CS", "IDH1", "IDH2", "OGDH", "FH", "MDH2", "ACO2")
  other_genes <- paste0("GENE", 1:100)
  
  all_genes <- c(glycolysis_genes, oxphos_genes, tca_genes, other_genes)
  n_genes <- length(all_genes)
  
  expr_matrix <- Matrix::Matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 5),
    nrow = n_genes,
    ncol = n_cells,
    sparse = TRUE
  )
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix)
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  
  seurat_obj@meta.data$cell_type <- sample(
    c("T_cell", "B_cell", "Macrophage"),
    n_cells, replace = TRUE
  )
  
  return(seurat_obj)
}

# ============================================================================
# TEST: Warburg Effect Simulation
# ============================================================================

test_that("Warburg effect simulation (cancer-like glycolytic phenotype)", {
  skip_if_missing()
  
  # Simulate cancer cells with high glycolysis, low OXPHOS
  n_cells <- 100
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  # Simulate Warburg phenotype
  glycolysis_signal <- rnorm(n_cells, mean = 2.5, sd = 0.5)
  oxphos_signal <- rnorm(n_cells, mean = 0.5, sd = 0.3)
  
  # Add to expression
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2", "PGK1")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  
  for (i in 1:n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- glycolysis_signal[i]
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- oxphos_signal[i]
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data", 
                            assay = "RNA", new.data = expr_data)
  
  # Score pathways
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = c("HK2", "LDHA", "PKM2", "PGK1"),
                                     oxphos = c("NDUFA9", "SDHA", "COX5A")
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Test expectations
  expect_true("glycolysis" %in% colnames(seurat_obj@meta.data))
  expect_true("oxphos" %in% colnames(seurat_obj@meta.data))
  
  # Glycolysis should be higher than OXPHOS
  expect_true(median(seurat_obj@meta.data$glycolysis) > 
              median(seurat_obj@meta.data$oxphos))
  
  # ATP score should reflect the metabolic state
  seurat_obj <- estimateATPProduction(seurat_obj, verbose = FALSE)
  expect_true("ATP_score" %in% colnames(seurat_obj@meta.data))
})

# ============================================================================
# TEST: Metabolic Switching Simulation
# ============================================================================

test_that("Metabolic switching from oxidative to glycolytic", {
  skip_if_missing()
  
  n_cells <- 150
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  # Simulate continuous transition
  progression <- seq(0, 1, length.out = n_cells)
  
  # Oxidative -> Glycolytic switch
  glycolysis_expr <- 0.5 + 2.5 * progression + rnorm(n_cells, 0, 0.2)
  oxphos_expr <- 2.5 - 2.0 * progression + rnorm(n_cells, 0, 0.2)
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  
  for (i in 1:n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- glycolysis_expr[i]
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- oxphos_expr[i]
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  # Score pathways
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = glycolytic_genes,
                                     oxphos = oxphos_genes
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Calculate GOX index
  seurat_obj <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  
  # GOX should increase with progression
  cor_with_progression <- cor(seurat_obj@meta.data$GOX_index, progression)
  expect_true(cor_with_progression > 0.9)
  
  # Classify phenotypes
  seurat_obj <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  
  # Early cells should be oxidative, late cells glycolytic
  early_gly <- sum(seurat_obj@meta.data$metabolic_phenotype[1:30] == "Oxidative")
  late_gly <- sum(seurat_obj@meta.data$metabolic_phenotype[121:150] == "Glycolytic")
  
  expect_true(early_gly > 10)
  expect_true(late_gly > 10)
})

# ============================================================================
# TEST: Cell Type-Specific Metabolism
# ============================================================================

test_that("Different cell types show distinct metabolic profiles", {
  skip_if_missing()
  
  n_t <- 30  # T cells - oxidative
  n_b <- 30  # B cells - glycolytic
  n_m <- 30  # Macrophages - mixed
  
  n_cells <- n_t + n_b + n_m
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  # T cells: high OXPHOS
  # B cells: high glycolysis  
  # Macrophages: mixed
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  
  # Set T cell expression
  for (i in 1:n_t) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 1.0, sd = 0.3)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 3.0, sd = 0.5)
    }
  }
  
  # Set B cell expression
  for (i in (n_t+1):(n_t+n_b)) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 3.5, sd = 0.5)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 0.8, sd = 0.3)
    }
  }
  
  # Set Macrophage expression
  for (i in (n_t+n_b+1):n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 2.0, sd = 0.4)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 2.0, sd = 0.4)
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = glycolytic_genes,
                                     oxphos = oxphos_genes
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Check that cell types have different profiles
  t_gly <- median(seurat_obj@meta.data$glycolysis[1:n_t])
  b_gly <- median(seurat_obj@meta.data$glycolysis[(n_t+1):(n_t+n_b)])
  
  expect_true(b_gly > t_gly)
  
  # T cells should have higher OXPHOS
  t_ox <- median(seurat_obj@meta.data$oxphos[1:n_t])
  b_ox <- median(seurat_obj@meta.data$oxphos[(n_t+1):(n_t+n_b)])
  
  expect_true(t_ox > b_ox)
})

# ============================================================================
# TEST: Disease vs Control Comparison
# ============================================================================

test_that("Disease shows metabolic dysregulation vs control", {
  skip_if_missing()
  
  n_control <- 50
  n_disease <- 50
  n_cells <- n_control + n_disease
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  # Control: normal metabolism
  # Disease: upregulated glycolysis, downregulated OXPHOS
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2", "PFKM")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1")
  
  # Control cells (first 50)
  for (i in 1:n_control) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 2.0, sd = 0.4)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 2.5, sd = 0.5)
    }
  }
  
  # Disease cells (next 50)
  for (i in (n_control+1):n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 4.0, sd = 0.6)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 1.2, sd = 0.4)
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  # Add condition metadata
  seurat_obj@meta.data$condition <- c(
    rep("Control", n_control),
    rep("Disease", n_disease)
  )
  
  # Score pathways
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = glycolytic_genes,
                                     oxphos = oxphos_genes
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Perform differential analysis
  results <- differentialMetabolicAnalysis(seurat_obj,
                                       condition_col = "condition",
                                       control_group = "Control",
                                       case_group = "Disease",
                                       method = "t.test")
  
  # Check results
  gly_result <- results$results[results$results$metric == "glycolysis", ]
  
  # Disease should have higher glycolysis
  expect_true(gly_result$case_mean > gly_result$control_mean)
  expect_true(gly_result$mean_diff > 0)
  
  # Effect size should be significant
  expect_true(gly_result$adj_p_value < 0.05)
})

# ============================================================================
# TEST: Flux Coupling Analysis
# ============================================================================

test_that("Metabolic fluxes are coupled correctly", {
  skip_if_missing()
  
  n_cells <- 80
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  # Simulate coupled glycolysis-TCA-OXPHOS
  glycolysis_signal <- rnorm(n_cells, mean = 2, sd = 0.5)
  tca_signal <- glycolysis_signal * 0.8 + rnorm(n_cells, 0, 0.2)
  oxphos_signal <- tca_signal * 0.9 + rnorm(n_cells, 0, 0.2)
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2")
  tca_genes <- c("CS", "IDH1", "IDH2")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  
  for (i in 1:n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- glycolysis_signal[i]
    }
    if (tca_genes[1] %in% rownames(expr_data)) {
      expr_data[tca_genes, i] <- tca_signal[i]
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- oxphos_signal[i]
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  # Compute flux
  seurat_obj <- computeMetabolicFlux(seurat_obj,
                                    gene_sets = list(
                                      glycolysis = glycolytic_genes,
                                      tca = tca_genes,
                                      oxphos = oxphos_genes
                                    ),
                                    verbose = FALSE)
  
  # Calculate coupling
  coupling <- calculateFluxCoupling(seurat_obj)
  
  # Expect positive correlation between coupled pathways
  expect_true(coupling["glycolysis_flux", "tca_flux"] > 0.5)
  expect_true(coupling["tca_flux", "oxphos_flux"] > 0.5)
})

# ============================================================================
# TEST: Energy Stress Response
# ============================================================================

test_that("Cells under energy stress show expected patterns", {
  skip_if_missing()
  
  n_cells <- 60
  
  # Normal cells (40)
  # Stressed cells (20) - low ATP
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  
  # Normal cells
  for (i in 1:40) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 2.0, sd = 0.3)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 2.0, sd = 0.3)
    }
  }
  
  # Stressed cells - low overall expression
  for (i in 41:n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 0.3, sd = 0.2)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 0.3, sd = 0.2)
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = glycolytic_genes,
                                     oxphos = oxphos_genes
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  seurat_obj <- estimateATPProduction(seurat_obj, verbose = FALSE)
  
  # Check ATP scores
  normal_atp <- median(seurat_obj@meta.data$ATP_score[1:40])
  stressed_atp <- median(seurat_obj@meta.data$ATP_score[41:n_cells])
  
  expect_true(normal_atp > stressed_atp)
})

# ============================================================================
# TEST: Hypoxia Response
# ============================================================================

test_that("Hypoxic conditions shift metabolism to glycolysis", {
  skip_if_missing()
  
  n_normoxic <- 40
  n_hypoxic <- 40
  n_cells <- n_normoxic + n_hypoxic
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2", "PFKM")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1")
  hypoxia_genes <- c("HIF1A", "VEGFA", "PDK1")
  
  # Normoxic: balanced metabolism
  for (i in 1:n_normoxic) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 2.0, sd = 0.4)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 2.5, sd = 0.5)
    }
  }
  
  # Hypoxic: high glycolysis, low OXPHOS
  for (i in (n_normoxic+1):n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 4.0, sd = 0.6)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 0.8, sd = 0.3)
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  seurat_obj@meta.data$oxygen <- c(
    rep("Normoxic", n_normoxic),
    rep("Hypoxic", n_hypoxic)
  )
  
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = glycolytic_genes,
                                     oxphos = oxphos_genes
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  # Compare groups
  normoxic_gly <- median(seurat_obj@meta.data$glycolysis[1:n_normoxic])
  hypoxic_gly <- median(seurat_obj@meta.data$glycolysis[(n_normoxic+1):n_cells])
  
  expect_true(hypoxic_gly > normoxic_gly)
  
  # Differential analysis
  results <- differentialMetabolicAnalysis(seurat_obj,
                                       condition_col = "oxygen",
                                       control_group = "Normoxic",
                                       case_group = "Hypoxic",
                                       method = "t.test")
  
  gly_result <- results$results[results$results$metric == "glycolysis", ]
  
  expect_true(gly_result$adj_p_value < 0.05)
})

# ============================================================================
# TEST: Batch Effect on Metabolism
# ============================================================================

test_that("Batch correction preserves biological signal", {
  skip_if_missing()
  
  n_batch1 <- 30
  n_batch2 <- 30
  n_cells <- n_batch1 + n_batch2
  
  seurat_obj <- .create_test_seurat_for_sim(n_cells)
  
  # Batch 1 has overall higher expression
  # Batch 2 has overall lower expression
  # But both should show same biological pattern
  
  expr_data <- GetAssayData(seurat_obj)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  
  batch1_factor <- 2.5
  batch2_factor <- 1.0
  
  for (i in 1:n_batch1) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 2.5 * batch1_factor, sd = 0.4)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 2.5 * batch1_factor, sd = 0.5)
    }
  }
  
  for (i in (n_batch1+1):n_cells) {
    if (glycolytic_genes[1] %in% rownames(expr_data)) {
      expr_data[glycolytic_genes, i] <- rnorm(1, mean = 2.5 * batch2_factor, sd = 0.4)
    }
    if (oxphos_genes[1] %in% rownames(expr_data)) {
      expr_data[oxphos_genes, i] <- rnorm(1, mean = 2.5 * batch2_factor, sd = 0.5)
    }
  }
  
  seurat_obj <- SetAssayData(seurat_obj, layer = "data",
                            assay = "RNA", new.data = expr_data)
  
  seurat_obj@meta.data$batch <- c(
    rep("Batch1", n_batch1),
    rep("Batch2", n_batch2)
  )
  
  seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = list(
                                     glycolysis = glycolytic_genes,
                                     oxphos = oxphos_genes
                                   ),
                                   method = "mean",
                                   verbose = FALSE)
  
  # After z-scoring, batches should be comparable
  seurat_obj <- zscorePathwayScores(seurat_obj)
  
  batch1_gly_zscore <- mean(seurat_obj@meta.data$zscore_glycolysis[1:n_batch1])
  batch2_gly_zscore <- mean(seurat_obj@meta.data$zscore_glycolysis[(n_batch1+1):n_cells])
  
  # Z-scores should be centered around 0
  expect_true(abs(batch1_gly_zscore) < 1.0)
  expect_true(abs(batch2_gly_zscore) < 1.0)
})

# ============================================================================
# FINISH
# ============================================================================

cat("\n[scMetaboFlux] Simulation validation tests completed!\n")
