#' @title Generate Example Single-Cell Dataset
#'
#' @description Creates a simulated single-cell dataset for demonstration purposes.
#' The dataset mimics PBMC-like data with distinct metabolic profiles across cell types.
#'
#' @param n_cells Number of cells (default 500)
#' @param n_genes Number of genes (default 1000)
#' @param seed Random seed for reproducibility
#'
#' @return A Seurat object with simulated expression data
#'
#' @export
#' @examples
#' \dontrun{
#' example_data <- generateExampleData(n_cells = 200, n_genes = 500)
#' }
generateExampleData <- function(n_cells = 500, n_genes = 1000, seed = 42) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required. Install with: BiocManager::install('Seurat')")
  }
  
  set.seed(seed)
  
  cell_types <- c("Naive_CD4_T", "Memory_CD4_T", "CD8_T", "B_cells", "Monocytes")
  n_cell_types <- length(cell_types)
  cells_per_type <- rep(floor(n_cells / n_cell_types), n_cell_types)
  cells_per_type[1] <- n_cells - sum(cells_per_type[-1])
  
  all_cell_types <- rep(cell_types, times = cells_per_type)
  
  metabolic_genes <- c(
    "HK2", "HK3", "HKDC1", "GPI", "PFKP", "ALDOA", "ALDOB", "ALDOC",
    "TPI1", "GAPDH", "PGK1", "PGK2", "PGAM1", "PGAM2", "ENO1", "ENO2",
    "ENO3", "PKM2", "PKM1", "LDHA", "LDHB", "LDHC", "PDHA1", "PDHA2",
    "PDHB", "DLAT", "DLD", "PDP1", "PDP2", "CS", "ACO2", "IDH3A", "IDH3B",
    "OGDH", "SUCLA2", "SUCLG1", "SUCLG2", "SDHA", "SDHB", "SDHC", "SDHD",
    "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7",
    "NDUFA8", "NDUFA9", "NDUFA10", "NDUFV1", "NDUFV2", "NDUFS1", "NDUFS2",
    "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFB1",
    "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C",
    "COX7A2", "COX7B", "COX7C", "COX8A", "COX8C", "ATP5F1A", "ATP5F1B",
    "ATP5F1C", "ATP5F1D", "ATP5F1E", "MTOR", "AMPK", "PIK3CA", "AKT1",
    "AKT2", "AKT3", "PTEN", "HIF1A", "VEGFA", "PDK1", "PFKFB3"
  )
  metabolic_genes <- unique(metabolic_genes)
  
  n_extra <- max(0, n_genes - length(metabolic_genes))
  if (n_extra > 0) {
    gene_names <- c(metabolic_genes, paste0("GENE", seq_len(n_extra)))
  } else {
    gene_names <- metabolic_genes[seq_len(n_genes)]
  }
  
  counts_matrix <- matrix(
    rnorm(n_genes * n_cells, mean = 2, sd = 1.5),
    nrow = n_genes,
    ncol = n_cells
  )
  counts_matrix <- pmax(counts_matrix, 0.01)
  rownames(counts_matrix) <- gene_names
  colnames(counts_matrix) <- paste0("Cell_", seq_len(n_cells))
  
  for (i in seq_along(cell_types)) {
    ct <- cell_types[i]
    start_idx <- sum(cells_per_type[seq_len(i - 1)]) + 1
    end_idx <- sum(cells_per_type[seq_len(i)])
    cell_idx <- start_idx:end_idx
    
    set.seed(seed + i)
    switch(ct,
      "Naive_CD4_T" = {
        for (g in c("LDHA", "PKM2", "HK2")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 2.5
          }
        }
        for (g in c("NDUFA1", "COX4I1", "ATP5F1A")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 3.0
          }
        }
      },
      "Memory_CD4_T" = {
        for (g in c("HK2", "PFKP", "GAPDH")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 2.8
          }
        }
      },
      "CD8_T" = {
        for (g in c("PGK1", "ENO1", "TPI1")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 3.0
          }
        }
      },
      "B_cells" = {
        for (g in c("PDHA1", "IDH3A", "CS")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 2.7
          }
        }
      },
      "Monocytes" = {
        for (g in c("HIF1A", "VEGFA", "PFKFB3")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 3.2
          }
        }
      }
    )
  }
  
  sparse_counts <- Matrix::Matrix(counts_matrix, sparse = TRUE)
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = sparse_counts,
    project = "scMetaboFlux_Example",
    assay = "RNA",
    min.cells = 3,
    min.features = 100
  )
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 500, verbose = FALSE)
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  
  seurat_obj@meta.data$cell_type <- all_cell_types
  seurat_obj@meta.data$cluster_id <- as.numeric(factor(all_cell_types))
  seurat_obj@meta.data$batch <- rep(c("Batch1", "Batch2"), length.out = n_cells)
  
  return(seurat_obj)
}

#' @title Generate Cancer-like Single-Cell Dataset
#'
#' @description Creates a simulated single-cell dataset mimicking cancer data
#' with hypoxic regions and metabolic heterogeneity.
#'
#' @param n_cells Number of cells (default 1000)
#' @param n_genes Number of genes (default 2000)
#' @param seed Random seed for reproducibility
#'
#' @return A Seurat object with simulated cancer expression data
#'
#' @export
#' @examples
#' \dontrun{
#' cancer_data <- generateCancerData(n_cells = 500, n_genes = 1000)
#' }
generateCancerData <- function(n_cells = 1000, n_genes = 2000, seed = 123) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }
  
  set.seed(seed)
  
  regions <- c("Normoxic_Tumor", "Hypoxic_Tumor", "Stromal")
  cells_per_region <- c(floor(n_cells * 0.4), floor(n_cells * 0.4), n_cells - 2 * floor(n_cells * 0.4))
  all_regions <- rep(regions, times = cells_per_region)
  
  metabolic_genes <- c(
    "HK2", "HK3", "PFKP", "ALDOA", "TPI1", "GAPDH", "PGK1", "ENO1", "ENO2", "PKM2",
    "LDHA", "LDHB", "PDHA1", "PDHB", "DLAT", "CS", "IDH3A", "OGDH",
    "SDHA", "SDHB", "NDUFA1", "NDUFA9", "COX4I1", "ATP5F1A",
    "HIF1A", "VEGFA", "PFKFB3", "PDK1", "MTOR", "AMPK"
  )
  metabolic_genes <- unique(metabolic_genes)
  
  n_extra <- max(0, n_genes - length(metabolic_genes))
  if (n_extra > 0) {
    gene_names <- c(metabolic_genes, paste0("GENE", seq_len(n_extra)))
  } else {
    gene_names <- metabolic_genes[seq_len(n_genes)]
  }
  
  counts_matrix <- matrix(
    rnorm(n_genes * n_cells, mean = 2.5, sd = 2.0),
    nrow = n_genes,
    ncol = n_cells
  )
  counts_matrix <- pmax(counts_matrix, 0.01)
  rownames(counts_matrix) <- gene_names
  colnames(counts_matrix) <- paste0("CancerCell_", seq_len(n_cells))
  
  for (i in seq_along(regions)) {
    reg <- regions[i]
    start_idx <- sum(cells_per_region[seq_len(i - 1)]) + 1
    end_idx <- sum(cells_per_region[seq_len(i)])
    cell_idx <- start_idx:end_idx
    
    set.seed(seed + i + 100)
    switch(reg,
      "Normoxic_Tumor" = {
        for (g in c("NDUFA1", "COX4I1", "ATP5F1A", "PDHA1", "CS")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 3.0
          }
        }
      },
      "Hypoxic_Tumor" = {
        for (g in c("HIF1A", "VEGFA", "PFKFB3", "LDHA", "PKM2", "PDK1")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 3.5
          }
        }
        for (g in c("COX4I1", "NDUFA1")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 0.3
          }
        }
      },
      "Stromal" = {
        for (g in c("MTOR", "AMPK", "PIK3CA")) {
          if (g %in% gene_names) {
            counts_matrix[which(gene_names == g), cell_idx] <- 
              counts_matrix[which(gene_names == g), cell_idx] * 2.5
          }
        }
      }
    )
  }
  
  sparse_counts <- Matrix::Matrix(counts_matrix, sparse = TRUE)
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = sparse_counts,
    project = "scMetaboFlux_Cancer",
    assay = "RNA",
    min.cells = 3,
    min.features = 100
  )
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 500, verbose = FALSE)
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  
  seurat_obj@meta.data$region <- all_regions
  seurat_obj@meta.data$cluster_id <- as.numeric(factor(all_regions))
  seurat_obj@meta.data$condition <- ifelse(grepl("Tumor", all_regions), "Tumor", "Stromal")
  
  return(seurat_obj)
}

#' @title Generate Differentiation Time Series Dataset
#'
#' @description Creates a simulated dataset mimicking stem cell differentiation
#' with metabolic transitions over pseudotime.
#'
#' @param n_cells Number of cells (default 800)
#' @param n_genes Number of genes (default 1500)
#' @param seed Random seed for reproducibility
#'
#' @return A Seurat object with simulated differentiation data
#'
#' @export
#' @examples
#' \dontrun{
#' diff_data <- generateDifferentiationData(n_cells = 400, n_genes = 800)
#' }
generateDifferentiationData <- function(n_cells = 800, n_genes = 1500, seed = 456) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }
  
  set.seed(seed)
  
  stages <- c("Stem", "Progenitor", "Early_Diff", "Late_Diff", "Mature")
  cells_per_stage <- rep(floor(n_cells / length(stages)), length(stages))
  cells_per_stage[1] <- n_cells - sum(cells_per_stage[-1])
  all_stages <- rep(stages, times = cells_per_stage)
  
  metabolic_genes <- c(
    "HK2", "PFKP", "GAPDH", "PGK1", "ENO1", "PKM2", "LDHA",
    "PDHA1", "CS", "IDH3A", "SDHA", "COX4I1", "ATP5F1A",
    "MTOR", "AMPK", "HIF1A"
  )
  metabolic_genes <- unique(metabolic_genes)
  
  n_extra <- max(0, n_genes - length(metabolic_genes))
  if (n_extra > 0) {
    gene_names <- c(metabolic_genes, paste0("GENE", seq_len(n_extra)))
  } else {
    gene_names <- metabolic_genes[seq_len(n_genes)]
  }
  
  counts_matrix <- matrix(
    rnorm(n_genes * n_cells, mean = 2, sd = 1.2),
    nrow = n_genes,
    ncol = n_cells
  )
  counts_matrix <- pmax(counts_matrix, 0.01)
  rownames(counts_matrix) <- gene_names
  colnames(counts_matrix) <- paste0("DiffCell_", seq_len(n_cells))
  
  for (i in seq_along(stages)) {
    stage <- stages[i]
    start_idx <- sum(cells_per_stage[seq_len(i - 1)]) + 1
    end_idx <- sum(cells_per_stage[seq_len(i)])
    cell_idx <- start_idx:end_idx
    
    set.seed(seed + i + 200)
    glycolytic_weight <- 0.5 + (i - 1) * 0.25
    oxphos_weight <- 1.0 - (i - 1) * 0.15
    
    for (g in c("HK2", "PFKP", "GAPDH", "LDHA", "PKM2")) {
      if (g %in% gene_names) {
        counts_matrix[which(gene_names == g), cell_idx] <- 
          counts_matrix[which(gene_names == g), cell_idx] * (1 + glycolytic_weight)
      }
    }
    for (g in c("PDHA1", "CS", "SDHA", "COX4I1", "ATP5F1A")) {
      if (g %in% gene_names) {
        counts_matrix[which(gene_names == g), cell_idx] <- 
          counts_matrix[which(gene_names == g), cell_idx] * (1 + oxphos_weight)
      }
    }
  }
  
  sparse_counts <- Matrix::Matrix(counts_matrix, sparse = TRUE)
  
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = sparse_counts,
    project = "scMetaboFlux_Differentiation",
    assay = "RNA",
    min.cells = 3,
    min.features = 100
  )
  
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = 500, verbose = FALSE)
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  
  seurat_obj@meta.data$differentiation_stage <- all_stages
  seurat_obj@meta.data$pseudotime <- rep(seq(0, 1, length.out = n_cells), each = 1)
  seurat_obj@meta.data$cluster_id <- as.numeric(factor(all_stages))
  
  return(seurat_obj)
}
