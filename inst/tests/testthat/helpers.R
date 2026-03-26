# ============================================================================
# TEST HELPER FUNCTIONS FOR scMetaboFlux
# ============================================================================
# Shared helper functions for creating test data
# ============================================================================

#' Create a mock Seurat object with specified properties
#' 
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @param include_reduction Include PCA/UMAP reductions
#' @param include_phenotype Include metabolic phenotype
#' @param seed Random seed for reproducibility
#' @return Seurat object
#' 
#' @keywords internal
create_test_seurat <- function(n_cells = 100, n_genes = 500, include_reduction = TRUE,
                              include_phenotype = FALSE, seed = 42) {
  set.seed(seed)
  
  # Include metabolic genes
  glycolysis_genes <- c("HK2", "HK3", "GCK", "GPI", "PFKL", "PFKM", "PFKP",
                        "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "GAPDHS", 
                        "PGK1", "PGK2", "PGAM1", "PGAM2", "ENO1", "ENO2", 
                        "ENO3", "PKM1", "PKM2", "PKLR", "LDHA", "LDHB", "LDHC")
  
  oxphos_genes <- c("NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6",
                    "NDUFA7", "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12",
                    "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6",
                    "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6",
                    "SDHA", "SDHB", "SDHC", "SDHD", "COX4I1", "COX5A", "COX5B",
                    "COX6A1", "COX6B1", "COX6C", "COX7A1", "COX7A2", "COX7B",
                    "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E")
  
  tca_genes <- c("CS", "ACO1", "ACO2", "IDH1", "IDH2", "IDH3A", "IDH3B", 
                 "OGDH", "SDHA", "SDHB", "FH", "MDH1", "MDH2", "ME1", "ME2",
                 "PCK1", "PCK2", "GOT1", "GOT2", "GLUD1")
  
  other_genes <- paste0("GENE", 1:(n_genes - length(c(glycolysis_genes, oxphos_genes, tca_genes))))
  all_genes <- c(glycolysis_genes, oxphos_genes, tca_genes, other_genes)
  all_genes <- all_genes[1:n_genes]
  
  # Create expression matrix with realistic properties
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 5, prob = 0.3),
    nrow = n_genes,
    ncol = n_cells
  )
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Cell_", sprintf("%04d", 1:n_cells))
  
  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, min.features = 10)
  
  # Add metadata
  seurat_obj@meta.data$cell_type <- sample(
    c("T_cell", "B_cell", "Macrophage", "NK_cell", "Dendritic"),
    n_cells,
    replace = TRUE,
    prob = c(0.3, 0.2, 0.25, 0.15, 0.1)
  )
  
  seurat_obj@meta.data$condition <- sample(
    c("Control", "Disease", "Treatment"),
    n_cells,
    replace = TRUE,
    prob = c(0.4, 0.35, 0.25)
  )
  
  seurat_obj@meta.data$patient_id <- sample(
    paste0("Patient_", 1:5),
    n_cells,
    replace = TRUE
  )
  
  seurat_obj@meta.data$batch <- sample(
    paste0("Batch_", 1:3),
    n_cells,
    replace = TRUE
  )
  
  # Add metabolic scores for advanced tests
  seurat_obj@meta.data$glycolysis <- rnorm(n_cells, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(n_cells, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$ATP_score <- runif(n_cells, min = 0.3, max = 0.9)
  seurat_obj@meta.data$GOX_index <- rnorm(n_cells, mean = 0, sd = 0.5)
  
  if (include_phenotype) {
    seurat_obj@meta.data$metabolic_phenotype <- sample(
      c("Glycolytic", "Oxidative", "Energetically Balanced", 
        "Energy-Stressed", "Hypermetabolic"),
      n_cells,
      replace = TRUE,
      prob = c(0.2, 0.25, 0.35, 0.1, 0.1)
    )
  }
  
  # Add flux scores
  seurat_obj@meta.data$glycolysis_flux <- rnorm(n_cells, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation_flux <- rnorm(n_cells, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$tca_cycle_flux <- rnorm(n_cells, mean = 0.4, sd = 0.15)
  
  # Normalize and scale
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  # Add PCA and UMAP reductions
  if (include_reduction) {
    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE, npcs = min(20, n_genes))
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:10, verbose = FALSE)
  }
  
  return(seurat_obj)
}

#' Create sparse matrix for testing
#' 
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @param sparsity Fraction of zeros (0-1)
#' @param seed Random seed
#' @return Seurat object with sparse data
#' 
#' @keywords internal
create_sparse_seurat <- function(n_cells = 100, n_genes = 200, sparsity = 0.95, seed = 42) {
  set.seed(seed)
  
  genes <- c("HK2", "LDHA", "SDHA", "COX5A", paste0("GENE", 1:(n_genes-4)))
  
  # Create sparse matrix
  counts <- Matrix::rsparsematrix(n_genes, n_cells, density = 1 - sparsity)
  counts <- as(counts, "dgCMatrix")
  rownames(counts) <- genes
  colnames(counts) <- paste0("Cell_", 1:n_cells)
  
  seurat_obj <- CreateSeuratObject(counts = counts)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  return(seurat_obj)
}

#' Create imbalanced dataset
#' 
#' @param n_cells Vector of cell counts per type
#' @param seed Random seed
#' @return Seurat object with imbalanced cell types
#' 
#' @keywords internal
create_imbalanced_seurat <- function(n_cells = c(10, 50, 200, 100, 5), seed = 42) {
  set.seed(seed)
  
  cell_types <- c("T_cell", "B_cell", "Macrophage", "NK_cell", "Neutrophil")
  total_cells <- sum(n_cells)
  
  genes <- c("HK2", "LDHA", "SDHA", "COX5A", paste0("GENE", 1:100))
  
  expr_matrix <- matrix(
    rnbinom(length(genes) * total_cells, size = 2, mu = 5),
    nrow = length(genes),
    ncol = total_cells
  )
  rownames(expr_matrix) <- genes
  colnames(expr_matrix) <- paste0("Cell_", 1:total_cells)
  
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  
  # Assign cell types with imbalanced distribution
  cell_types_vec <- rep(cell_types, times = n_cells)
  seurat_obj@meta.data$cell_type <- cell_types_vec[1:ncol(seurat_obj)]
  
  seurat_obj@meta.data$glycolysis <- rnorm(total_cells, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(total_cells, mean = 0.5, sd = 0.2)
  seurat_obj@meta.data$ATP_score <- runif(total_cells, 0.3, 0.9)
  
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  return(seurat_obj)
}

#' Skip tests if required packages are missing
#' 
#' @keywords internal
skip_if_missing <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    skip("Seurat not available")
  }
  if (!requireNamespace("AUCell", quietly = TRUE)) {
    skip("AUCell not available")
  }
  if (!requireNamespace("GSVA", quietly = TRUE)) {
    skip("GSVA not available")
  }
}

#' Create mock data for testing with specific metabolic signal
#' 
#' @param n_cells Number of cells
#' @param glycolysis_level Mean glycolysis expression level
#' @param oxphos_level Mean OXPHOS expression level
#' @param seed Random seed
#' @return Seurat object
#' 
#' @keywords internal
create_seurat_with_metabolism <- function(n_cells = 50, 
                                        glycolysis_level = 2,
                                        oxphos_level = 2,
                                        seed = 42) {
  set.seed(seed)
  
  glycolytic_genes <- c("HK2", "LDHA", "PKM2", "PGK1", "ENO1", "GAPDH", "TPI1", "PFKM")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A", "ATP5A1", "COX4I1", "NDUFS1", "SDHB")
  other_genes <- paste0("GENE", 1:50)
  
  all_genes <- c(glycolytic_genes, oxphos_genes, other_genes)
  
  n_genes <- length(all_genes)
  
  # Base expression
  expr_matrix <- matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 2),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Add metabolic signal
  for (gene in glycolytic_genes) {
    if (gene %in% all_genes) {
      idx <- which(all_genes == gene)
      expr_matrix[idx, ] <- rnorm(n_cells, mean = glycolysis_level, sd = 0.5)
      expr_matrix[idx, ] <- pmax(0, expr_matrix[idx, ])
    }
  }
  
  for (gene in oxphos_genes) {
    if (gene %in% all_genes) {
      idx <- which(all_genes == gene)
      expr_matrix[idx, ] <- rnorm(n_cells, mean = oxphos_level, sd = 0.5)
      expr_matrix[idx, ] <- pmax(0, expr_matrix[idx, ])
    }
  }
  
  rownames(expr_matrix) <- all_genes
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  return(seurat_obj)
}

#' Helper to check if ggplot is valid
#' 
#' @param p Plot object
#' @return Logical
#' 
#' @keywords internal
is_valid_ggplot <- function(p) {
  inherits(p, "ggplot")
}

#' Helper to create condition with specific effect
#' 
#' @param seurat_obj Seurat object
#' @param condition_col Column name for condition
#' @param case_value Case value
#' @param control_value Control value
#' @param effect Effect size (difference in means)
#' @param score_col Score column
#' @return Seurat object with modified scores
#' 
#' @keywords internal
add_condition_effect <- function(seurat_obj, condition_col = "condition",
                               case_value = "Case", control_value = "Control",
                               effect = 0.5, score_col = "glycolysis") {
  case_mask <- seurat_obj@meta.data[[condition_col]] == case_value
  control_mask <- seurat_obj@meta.data[[condition_col]] == control_value
  
  seurat_obj@meta.data[[score_col]][case_mask] <- 
    seurat_obj@meta.data[[score_col]][case_mask] + effect
  
  return(seurat_obj)
}

# Utility to check matrix properties
.check_matrix_properties <- function(matrix_obj) {
  list(
    is_sparse = inherits(matrix_obj, "sparseMatrix"),
    is_dgCMatrix = inherits(matrix_obj, "dgCMatrix"),
    n_rows = nrow(matrix_obj),
    n_cols = ncol(matrix_obj),
    min_val = min(matrix_obj),
    max_val = max(matrix_obj),
    mean_val = mean(matrix_obj),
    has_na = any(is.na(matrix_obj)),
    has_inf = any(is.infinite(matrix_obj))
  )
}
