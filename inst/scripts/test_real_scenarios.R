# ============================================================================
# COMPREHENSIVE REAL DATA TEST SCENARIOS FOR scMetaboFlux
# ============================================================================
# This script tests the package with various realistic scenarios
# ============================================================================

# Load required packages
library(Seurat)
library(scMetaboFlux)

# Set seed for reproducibility
set.seed(42)

# ============================================================================
# SCENARIO 1: PBMC Dataset (Standard 10X Data)
# ============================================================================
# This simulates a typical PBMC dataset from 10X Genomics

test_pbmc_scenario <- function() {
  cat("\n=== SCENARIO 1: PBMC Dataset ===\n")
  
  # Create synthetic PBMC-like data
  n_cells <- 500
  n_genes <- 1000
  
  # Cell type proportions typical of PBMCs
  cell_types <- c(
    rep("CD4_T", 150),
    rep("CD8_T", 100),
    rep("B_cell", 50),
    rep("NK", 75),
    rep("Monocyte", 75),
    rep("DC", 25),
    rep("Platelet", 25)
  )
  
  # Known metabolic genes in PBMCs
  glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH", "TPI1", "PFKM")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1", "NDUFS1", "SDHB", "COX6C")
  tca_genes <- c("CS", "IDH1", "IDH2", "OGDH", "FH", "MDH2", "ACO2", "IDH3A")
  
  all_metabolic_genes <- c(glycolysis_genes, oxphos_genes, tca_genes)
  other_genes <- paste0("GENE", 1:(n_genes - length(all_metabolic_genes)))
  all_genes <- c(all_metabolic_genes, other_genes)
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 3, prob = 0.5),
    nrow = n_genes,
    ncol = n_cells
  )
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("PBMC_", 1:n_cells)
  
  # Add cell-type specific metabolic signatures
  # T cells: higher OXPHOS
  t_cell_idx <- which(cell_types %in% c("CD4_T", "CD8_T"))
  for (gene in oxphos_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, t_cell_idx] <- rnbinom(length(t_cell_idx), size = 3, mu = 8, prob = 0.3)
    }
  }
  
  # B cells: moderate glycolysis
  b_cell_idx <- which(cell_types == "B_cell")
  for (gene in glycolysis_genes[1:4]) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, b_cell_idx] <- rnbinom(length(b_cell_idx), size = 3, mu = 6, prob = 0.4)
    }
  }
  
  # Monocytes: high glycolysis (known biology)
  mono_idx <- which(cell_types == "Monocyte")
  for (gene in glycolysis_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, mono_idx] <- rnbinom(length(mono_idx), size = 3, mu = 7, prob = 0.35)
    }
  }
  
  # Create Seurat object
  pbmc <- CreateSeuratObject(counts = expr_matrix, min.features = 100)
  pbmc <- NormalizeData(pbmc)
  pbmc <- ScaleData(pbmc)
  pbmc <- RunPCA(pbmc, npcs = 20, verbose = FALSE)
  pbmc <- RunUMAP(pbmc, dims = 1:20, verbose = FALSE)
  
  pbmc@meta.data$cell_type <- cell_types
  pbmc@meta.data$condition <- sample(c("Control", "Stimulated"), n_cells, replace = TRUE)
  
  cat("Created PBMC dataset with", n_cells, "cells and", n_genes, "genes\n")
  
  # Run analysis
  cat("Running metabolic analysis...\n")
  
  result <- tryCatch({
    pbmc <- computePathwayScores(pbmc, method = "mean", verbose = FALSE)
    pbmc <- estimateATPProduction(pbmc, verbose = FALSE)
    pbmc <- calculateGOXIndex(pbmc, verbose = FALSE)
    pbmc <- computeMetabolicFlux(pbmc, verbose = FALSE)
    pbmc <- classifyMetabolicPhenotype(pbmc, verbose = FALSE)
    pbmc <- aggregateByCellType(pbmc, cell_type_col = "cell_type")
    
    list(success = TRUE, object = pbmc)
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    list(success = FALSE, error = e)
  })
  
  return(result)
}

# ============================================================================
# SCENARIO 2: Cancer Dataset (Tumor vs Normal)
# ============================================================================
# Simulates tumor microenvironment with metabolic heterogeneity

test_cancer_scenario <- function() {
  cat("\n=== SCENARIO 2: Cancer Dataset ===\n")
  
  n_tumor <- 300
  n_normal <- 200
  n_cells <- n_tumor + n_normal
  n_genes <- 1200
  
  # Cell types in tumor microenvironment
  cell_types <- c(
    rep("Cancer_cell", 200),
    rep("T_cell", 100),
    rep("Macrophage", 80),
    rep("CAF", 60),
    rep("Endothelial", 40),
    rep("Normal_epithelial", 20)
  )
  
  glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH", "TPI1", "PFKM", "PFKP", "ALDOA")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1", "NDUFS1", "SDHB", "COX6C", "ATP5F1B", "NDUFV2")
  hypoxia_genes <- c("HIF1A", "VEGFA", "PDK1", "BNIP3", "LDHA")
  
  all_metabolic_genes <- c(glycolysis_genes, oxphos_genes, hypoxia_genes)
  other_genes <- paste0("GENE", 1:(n_genes - length(all_metabolic_genes)))
  all_genes <- c(all_metabolic_genes, other_genes)
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 3, prob = 0.5),
    nrow = n_genes,
    ncol = n_cells
  )
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Tumor_", 1:n_cells)
  
  # Tumor cells: Warburg effect (high glycolysis, low OXPHOS)
  tumor_idx <- 1:n_tumor
  for (gene in glycolysis_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, tumor_idx] <- rnbinom(n_tumor, size = 3, mu = 10, prob = 0.25)
    }
  }
  for (gene in oxphos_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, tumor_idx] <- rnbinom(n_tumor, size = 3, mu = 3, prob = 0.6)
    }
  }
  for (gene in hypoxia_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, tumor_idx] <- rnbinom(n_tumor, size = 3, mu = 7, prob = 0.3)
    }
  }
  
  # Normal cells: oxidative metabolism
  normal_idx <- (n_tumor + 1):n_cells
  for (gene in oxphos_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, normal_idx] <- rnbinom(n_normal, size = 3, mu = 8, prob = 0.3)
    }
  }
  for (gene in glycolysis_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, normal_idx] <- rnbinom(n_normal, size = 3, mu = 4, prob = 0.5)
    }
  }
  
  # Create Seurat object
  tumor <- CreateSeuratObject(counts = expr_matrix, min.features = 100)
  tumor <- NormalizeData(tumor)
  tumor <- ScaleData(tumor)
  tumor <- RunPCA(tumor, npcs = 20, verbose = FALSE)
  tumor <- RunUMAP(tumor, dims = 1:20, verbose = FALSE)
  
  tumor@meta.data$cell_type <- cell_types
  tumor@meta.data$condition <- c(rep("Tumor", n_tumor), rep("Normal", n_normal))
  tumor@meta.data$patient_id <- sample(paste0("Patient_", 1:5), n_cells, replace = TRUE)
  
  cat("Created tumor dataset with", n_cells, "cells\n")
  
  # Run analysis
  cat("Running metabolic analysis...\n")
  
  result <- tryCatch({
    tumor <- computePathwayScores(tumor, method = "mean", verbose = FALSE)
    tumor <- estimateATPProduction(tumor, verbose = FALSE)
    tumor <- calculateGOXIndex(tumor, verbose = FALSE)
    tumor <- computeMetabolicFlux(tumor, verbose = FALSE)
    tumor <- classifyMetabolicPhenotype(tumor, verbose = FALSE)
    
    # Differential analysis
    diff_results <- differentialMetabolicAnalysis(tumor,
                                                 condition_col = "condition",
                                                 control_group = "Normal",
                                                 case_group = "Tumor")
    
    # Cell type comparison
    comparison <- compareCellTypesMetabolism(tumor, cell_type_col = "cell_type")
    
    list(
      success = TRUE,
      object = tumor,
      differential = diff_results,
      comparison = comparison
    )
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    traceback()
    list(success = FALSE, error = e)
  })
  
  return(result)
}

# ============================================================================
# SCENARIO 3: Development/ Differentiation Time Course
# ============================================================================
# Simulates metabolic changes during differentiation

test_differentiation_scenario <- function() {
  cat("\n=== SCENARIO 3: Differentiation Time Course ===\n")
  
  n_timepoints <- 6
  n_cells_per_time <- 50
  n_total <- n_timepoints * n_cells_per_time
  n_genes <- 800
  
  timepoints <- rep(paste0("Day_", 0:(n_timepoints-1)), each = n_cells_per_time)
  
  glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1")
  
  all_metabolic_genes <- c(glycolysis_genes, oxphos_genes)
  other_genes <- paste0("GENE", 1:(n_genes - length(all_metabolic_genes)))
  all_genes <- c(all_metabolic_genes, other_genes)
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_genes * n_total, size = 2, mu = 3, prob = 0.5),
    nrow = n_genes,
    ncol = n_total
  )
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Diff_", 1:n_total)
  
  # Simulate metabolic switch during differentiation
  # Early: high glycolysis (proliferative state)
  # Late: high OXPHOS (mature state)
  
  progression <- rep(0:(n_timepoints-1), each = n_cells_per_time)
  
  for (i in 1:n_total) {
    prog <- progression[i]
    gly_level <- 8 - prog * 1.0  # Decreasing glycolysis
    ox_level <- 3 + prog * 1.0   # Increasing OXPHOS
    
    for (gene in glycolysis_genes) {
      if (gene %in% rownames(expr_matrix)) {
        expr_matrix[gene, i] <- rnbinom(1, size = 3, mu = max(2, gly_level), prob = 0.4)
      }
    }
    for (gene in oxphos_genes) {
      if (gene %in% rownames(expr_matrix)) {
        expr_matrix[gene, i] <- rnbinom(1, size = 3, mu = max(2, ox_level), prob = 0.4)
      }
    }
  }
  
  # Create Seurat object
  diff_data <- CreateSeuratObject(counts = expr_matrix, min.features = 50)
  diff_data <- NormalizeData(diff_data)
  diff_data <- ScaleData(diff_data)
  diff_data <- RunPCA(diff_data, npcs = 10, verbose = FALSE)
  diff_data <- RunUMAP(diff_data, dims = 1:10, verbose = FALSE)
  
  diff_data@meta.data$timepoint <- timepoints
  diff_data@meta.data$pseudotime <- rep(0:(n_timepoints-1), each = n_cells_per_time)
  
  cat("Created differentiation dataset with", n_total, "cells across", n_timepoints, "timepoints\n")
  
  # Run analysis
  cat("Running metabolic analysis...\n")
  
  result <- tryCatch({
    diff_data <- computePathwayScores(diff_data, method = "mean", verbose = FALSE)
    diff_data <- estimateATPProduction(diff_data, verbose = FALSE)
    diff_data <- calculateGOXIndex(diff_data, verbose = FALSE)
    diff_data <- computeMetabolicFlux(diff_data, verbose = FALSE)
    diff_data <- classifyMetabolicPhenotype(diff_data, verbose = FALSE)
    
    # Aggregate by timepoint
    agg_stats <- aggregateByCellType(diff_data, cell_type_col = "timepoint")
    
    list(
      success = TRUE,
      object = diff_data,
      aggregation = agg_stats
    )
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    list(success = FALSE, error = e)
  })
  
  return(result)
}

# ============================================================================
# SCENARIO 4: Hypoxia Treatment Experiment
# ============================================================================
# Simulates metabolic response to oxygen levels

test_hypoxia_scenario <- function() {
  cat("\n=== SCENARIO 4: Hypoxia Treatment ===\n")
  
  n_normoxic <- 100
  n_hypoxic <- 100
  n_cells <- n_normoxic + n_hypoxic
  n_genes <- 600
  
  glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PFKM", "PGK1", "ENO1", "TPI1")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1", "NDUFS1")
  hypoxia_genes <- c("HIF1A", "VEGFA", "PDK1", "BNIP3", "EGLN1", "SLC2A1")
  
  all_metabolic_genes <- c(glycolysis_genes, oxphos_genes, hypoxia_genes)
  other_genes <- paste0("GENE", 1:(n_genes - length(all_metabolic_genes)))
  all_genes <- c(all_metabolic_genes, other_genes)
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 3, prob = 0.5),
    nrow = n_genes,
    ncol = n_cells
  )
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Hypoxia_", 1:n_cells)
  
  # Normoxic cells: balanced metabolism
  for (gene in c(glycolysis_genes, oxphos_genes)) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, 1:n_normoxic] <- rnbinom(n_normoxic, size = 3, mu = 5, prob = 0.4)
    }
  }
  
  # Hypoxic cells: glycolytic shift + HIF1A activation
  for (gene in glycolysis_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, (n_normoxic+1):n_cells] <- rnbinom(n_hypoxic, size = 3, mu = 9, prob = 0.25)
    }
  }
  for (gene in oxphos_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, (n_normoxic+1):n_cells] <- rnbinom(n_hypoxic, size = 3, mu = 2, prob = 0.7)
    }
  }
  for (gene in hypoxia_genes) {
    if (gene %in% rownames(expr_matrix)) {
      expr_matrix[gene, (n_normoxic+1):n_cells] <- rnbinom(n_hypoxic, size = 3, mu = 10, prob = 0.2)
    }
  }
  
  # Create Seurat object
  hypoxia_exp <- CreateSeuratObject(counts = expr_matrix, min.features = 50)
  hypoxia_exp <- NormalizeData(hypoxia_exp)
  hypoxia_exp <- ScaleData(hypoxia_exp)
  hypoxia_exp <- RunPCA(hypoxia_exp, npcs = 15, verbose = FALSE)
  hypoxia_exp <- RunUMAP(hypoxia_exp, dims = 1:15, verbose = FALSE)
  
  hypoxia_exp@meta.data$condition <- c(rep("Normoxic", n_normoxic), rep("Hypoxic", n_hypoxic))
  hypoxia_exp@meta.data$cell_type <- sample(c("Type_A", "Type_B"), n_cells, replace = TRUE)
  
  cat("Created hypoxia dataset with", n_cells, "cells\n")
  
  # Run analysis
  cat("Running metabolic analysis...\n")
  
  result <- tryCatch({
    hypoxia_exp <- computePathwayScores(hypoxia_exp, method = "mean", verbose = FALSE)
    hypoxia_exp <- estimateATPProduction(hypoxia_exp, verbose = FALSE)
    hypoxia_exp <- calculateGOXIndex(hypoxia_exp, verbose = FALSE)
    hypoxia_exp <- computeMetabolicFlux(hypoxia_exp, verbose = FALSE)
    hypoxia_exp <- classifyMetabolicPhenotype(hypoxia_exp, verbose = FALSE)
    
    # Differential analysis
    diff_results <- differentialMetabolicAnalysis(hypoxia_exp,
                                                 condition_col = "condition",
                                                 control_group = "Normoxic",
                                                 case_group = "Hypoxic")
    
    list(
      success = TRUE,
      object = hypoxia_exp,
      differential = diff_results
    )
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    list(success = FALSE, error = e)
  })
  
  return(result)
}

# ============================================================================
# SCENARIO 5: Large Scale Dataset (Performance Test)
# ============================================================================
# Tests performance with larger datasets

test_large_scale_scenario <- function() {
  cat("\n=== SCENARIO 5: Large Scale Dataset ===\n")
  
  n_cells <- 2000
  n_genes <- 1500
  
  glycolysis_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1")
  
  all_metabolic_genes <- c(glycolysis_genes, oxphos_genes)
  other_genes <- paste0("GENE", 1:(n_genes - length(all_metabolic_genes)))
  all_genes <- c(all_metabolic_genes, other_genes)
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 3, prob = 0.5),
    nrow = n_genes,
    ncol = n_cells
  )
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Large_", 1:n_cells)
  
  # Create Seurat object
  large_data <- CreateSeuratObject(counts = expr_matrix, min.features = 100)
  large_data <- NormalizeData(large_data)
  large_data <- ScaleData(large_data)
  large_data <- RunPCA(large_data, npcs = 30, verbose = FALSE)
  large_data <- RunUMAP(large_data, dims = 1:30, verbose = FALSE)
  
  large_data@meta.data$cell_type <- sample(paste0("CellType_", 1:10), n_cells, replace = TRUE)
  
  cat("Created large dataset with", n_cells, "cells\n")
  
  # Time the analysis
  cat("Running metabolic analysis...\n")
  
  start_time <- Sys.time()
  
  result <- tryCatch({
    large_data <- computePathwayScores(large_data, method = "mean", verbose = FALSE)
    large_data <- estimateATPProduction(large_data, verbose = FALSE)
    large_data <- calculateGOXIndex(large_data, verbose = FALSE)
    large_data <- computeMetabolicFlux(large_data, verbose = FALSE)
    large_data <- classifyMetabolicPhenotype(large_data, verbose = FALSE)
    
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    list(
      success = TRUE,
      object = large_data,
      elapsed_time = elapsed
    )
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
    list(success = FALSE, error = e)
  })
  
  return(result)
}

# ============================================================================
# RUN ALL SCENARIOS
# ============================================================================

cat("========================================\n")
cat("scMetaboFlux Real Data Test Scenarios\n")
cat("========================================\n")

# Run all scenarios
scenarios <- list(
  "PBMC" = test_pbmc_scenario(),
  "Cancer" = test_cancer_scenario(),
  "Differentiation" = test_differentiation_scenario(),
  "Hypoxia" = test_hypoxia_scenario(),
  "LargeScale" = test_large_scale_scenario()
)

# Summary
cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")

for (name in names(scenarios)) {
  result <- scenarios[[name]]
  status <- ifelse(result$success, "PASSED", "FAILED")
  cat(sprintf("%-15s: %s\n", name, status))
  
  if (result$success && !is.null(result$elapsed_time)) {
    cat(sprintf("  Time: %.2f seconds\n", result$elapsed_time))
  }
  
  if (!result$success) {
    cat("  Error:", conditionMessage(result$error), "\n")
  }
}

# Final check
all_passed <- all(sapply(scenarios, function(x) x$success))

if (all_passed) {
  cat("\n*** ALL SCENARIOS PASSED ***\n")
} else {
  cat("\n*** SOME SCENARIOS FAILED ***\n")
}

invisible(scenarios)
