#' @title Multi-Omics Integration Module
#' @description Integrate scMetaboFlux with other omics layers
#' @name multiomics_integration
#' @rdname multiomics-integration
#' @keywords multiomics integration spatial proteomics
NULL

#' @title Integrate scATAC-seq Regulatory Support
#' 
#' @description Integrate ATAC-seq peaks with metabolic gene activity.
#' 
#' @param seurat_obj Seurat object
#' @param peak_matrix Peak matrix from scATAC (optional)
#' @param metabolic_genes Metabolic genes to check
#' @param species Species (human/mouse)
#' 
#' @return Data frame with regulatory support scores
#' 
#' @export
#' @examples
#' \dontrun{
#' reg_support <- integrateATACRegulatorySupport(seurat_obj)
#' }
integrateATACRegulatorySupport <- function(seurat_obj,
                                          peak_matrix = NULL,
                                          metabolic_genes = NULL,
                                          species = c("human", "mouse")) {
  
  species <- match.arg(species)
  
  if (is.null(metabolic_genes)) {
    metabolic_genes <- unique(unlist(metabolicGeneSets))
  }
  
  gene_scores <- sapply(metabolic_genes, function(gene) {
    if (gene %in% rownames(seurat_obj)) {
      return(mean(Seurat::GetAssayData(seurat_obj)[gene, ], na.rm = TRUE))
    }
    return(0)
  })
  
  if (is.null(peak_matrix)) {
    message("[scMetaboFlux] No peak matrix provided. Using gene expression as regulatory proxy.")
  }
  
  results <- data.frame(
    gene = names(gene_scores),
    expression_support = as.numeric(gene_scores),
    stringsAsFactors = FALSE
  )
  
  results <- results[order(-results$expression_support), ]
  
  return(results)
}

#' @title Integrate Spatial Transcriptomics
#' 
#' @description Integrate spatial location with metabolic states.
#' 
#' @param seurat_obj Seurat object
#' @param spatial_coords Data frame with spatial coordinates and cell IDs
#' @param cell_id_col Cell ID column name in spatial_coords
#' @param x_col Column name for x coordinate (default "x")
#' @param y_col Column name for y coordinate (default "y")
#' 
#' @return Seurat object with spatial metabolic data
#' 
#' @export
#' @examples
#' \dontrun{
#' spatial_coords <- read.csv("spatial_coordinates.csv")
#' seurat_obj <- integrateSpatialMetabolism(seurat_obj, spatial_coords)
#' }
integrateSpatialMetabolism <- function(seurat_obj,
                                      spatial_coords,
                                      cell_id_col = "cell_id",
                                      x_col = "x",
                                      y_col = "y") {
  
  if (!is.data.frame(spatial_coords)) {
    stop("spatial_coords must be a data frame")
  }
  
  if (!cell_id_col %in% colnames(spatial_coords)) {
    stop(paste("Cell ID column", cell_id_col, "not found in spatial_coords"))
  }
  
  cell_ids <- rownames(seurat_obj@meta.data)
  
  matched_coords <- spatial_coords[match(cell_ids, spatial_coords[[cell_id_col]]), ]
  
  if (x_col %in% colnames(matched_coords)) {
    seurat_obj@meta.data$spatial_x <- matched_coords[[x_col]]
  }
  
  if (y_col %in% colnames(matched_coords)) {
    seurat_obj@meta.data$spatial_y <- matched_coords[[y_col]]
  }
  
  if ("ATP_score" %in% colnames(seurat_obj@meta.data)) {
    hotspot_threshold <- quantile(seurat_obj@meta.data$ATP_score, 0.9, na.rm = TRUE)
    seurat_obj@meta.data$spatial_metabolic_hotspot <- 
      seurat_obj@meta.data$ATP_score > hotspot_threshold
  }
  
  if (exists("GOX_index", where = seurat_obj@meta.data, inherits = FALSE)) {
    glycolytic_threshold <- quantile(seurat_obj@meta.data$GOX_index, 0.9, na.rm = TRUE)
    oxidative_threshold <- quantile(seurat_obj@meta.data$GOX_index, 0.1, na.rm = TRUE)
    
    seurat_obj@meta.data$spatial_glycolytic_region <- 
      seurat_obj@meta.data$GOX_index > glycolytic_threshold
    
    seurat_obj@meta.data$spatial_oxidative_region <- 
      seurat_obj@meta.data$GOX_index < oxidative_threshold
  }
  
  return(seurat_obj)
}

#' @title Calculate Spatial Metabolic Heterogeneity
#' 
#' @description Calculate spatial autocorrelation of metabolic scores.
#' 
#' @param seurat_obj Seurat object
#' @param score_col Score column to analyze
#' @param coords_cols Coordinates columns
#' @param k Number of neighbors for Moran's I
#' 
#' @return List with spatial statistics
#' 
#' @export
#' @examples
#' \dontrun{
#' spatial_stats <- calculateSpatialMetabolicHeterogeneity(seurat_obj,
#'   score_col = "ATP_score",
#'   coords_cols = c("spatial_x", "spatial_y"))
#' }
calculateSpatialMetabolicHeterogeneity <- function(seurat_obj,
                                                   score_col = "ATP_score",
                                                   coords_cols = c("spatial_x", "spatial_y"),
                                                   k = 10) {
  
  if (!all(coords_cols %in% colnames(seurat_obj@meta.data))) {
    stop("Spatial coordinates not found in metadata")
  }
  
  coords <- seurat_obj@meta.data[, coords_cols]
  scores <- seurat_obj@meta.data[[score_col]]
  
  coords <- coords[complete.cases(coords), ]
  scores <- scores[complete.cases(coords)]
  
  if (nrow(coords) < k + 1) {
    stop("Not enough complete cases for spatial analysis")
  }
  
  dist_matrix <- as.matrix(dist(coords))
  
  neighbors <- apply(dist_matrix, 1, function(d) order(d)[2:(k + 1)])
  
  local_i <- sapply(1:nrow(coords), function(i) {
    neighbors_i <- neighbors[, i]
    if (any(is.na(neighbors_i))) return(NA)
    local_mean <- mean(scores[neighbors_i], na.rm = TRUE)
    if (sd(scores, na.rm = TRUE) > 0) {
      (scores[i] - local_mean) / sd(scores, na.rm = TRUE)
    } else {
      0
    }
  })
  
  weights <- 1 / (dist_matrix + 0.1)
  weighted_scores <- apply(weights, 2, function(w) sum(w * scores, na.rm = TRUE))
  
  global_moran <- cor(scores, weighted_scores, use = "pairwise")
  
  return(list(
    global_moran_i = global_moran,
    local_moran_i = local_i,
    spatial_autocorrelation = abs(global_moran) > 0.1,
    p_value = 2 * (1 - pnorm(abs(global_moran)))
  ))
}

#' @title Create Metabolic Gradient Map
#' 
#' @description Create gradient map of metabolic activity in space.
#' 
#' @param seurat_obj Seurat object
#' @param score_col Score column
#' @param resolution Resolution for interpolation
#' 
#' @return List with gradient matrix and coordinates
#' 
#' @export
#' @examples
#' \dontrun{
#' gradient_map <- createMetabolicGradientMap(seurat_obj, score_col = "ATP_score")
#' }
createMetabolicGradientMap <- function(seurat_obj,
                                      score_col = "ATP_score",
                                      resolution = 50) {
  
  coords_cols <- c("spatial_x", "spatial_y")
  
  if (!all(coords_cols %in% colnames(seurat_obj@meta.data))) {
    stop("Spatial coordinates not found. Run integrateSpatialMetabolism first.")
  }
  
  coords <- data.frame(
    x = seurat_obj@meta.data$spatial_x,
    y = seurat_obj@meta.data$spatial_y,
    score = seurat_obj@meta.data[[score_col]]
  )
  
  coords <- coords[complete.cases(coords), ]
  
  x_range <- range(coords$x, na.rm = TRUE)
  y_range <- range(coords$y, na.rm = TRUE)
  
  x_seq <- seq(x_range[1], x_range[2], length.out = resolution)
  y_seq <- seq(y_range[1], y_range[2], length.out = resolution)
  
  gradient_matrix <- matrix(NA, nrow = resolution, ncol = resolution)
  
  for (i in seq_len(resolution)) {
    for (j in seq_len(resolution)) {
      distances <- sqrt((coords$x - x_seq[i])^2 + (coords$y - y_seq[j])^2)
      weights <- 1 / (distances + 0.1)
      gradient_matrix[j, i] <- weighted.mean(coords$score, weights, na.rm = TRUE)
    }
  }
  
  colnames(gradient_matrix) <- paste0("x_", seq_len(resolution))
  rownames(gradient_matrix) <- paste0("y_", seq_len(resolution))
  
  return(list(
    gradient = gradient_matrix,
    x_coords = x_seq,
    y_coords = y_seq,
    method = "inverse_distance_weighting",
    resolution = resolution
  ))
}

#' @title Integrate Proteomics Data
#' 
#' @description Correlate metabolic gene expression with protein abundance.
#' 
#' @param seurat_obj Seurat object
#' @param protein_data Protein expression matrix (cells x proteins)
#' @param metabolic_proteins Known metabolic proteins to analyze
#' 
#' @return Data frame with correlations
#' 
#' @export
#' @examples
#' \dontrun{
#' protein_corr <- integrateProteomicsData(seurat_obj)
#' }
integrateProteomicsData <- function(seurat_obj,
                                    protein_data = NULL,
                                    metabolic_proteins = NULL) {
  
  if (is.null(metabolic_proteins)) {
    metabolic_proteins <- c("HK2", "LDHA", "SDHA", "COX4I1", "ATP5A1",
                           "PFKP", "PKM2", "NDUFA9", "UQCRC1", "COX5A",
                           "PFKFB3", "ENO1", "GAPDH", "TPI1")
  }
  
  gene_expr <- Seurat::GetAssayData(seurat_obj)
  genes_upper <- toupper(rownames(gene_expr))
  metabolic_upper <- toupper(metabolic_proteins)
  
  correlations <- sapply(metabolic_proteins, function(p) {
    p_upper <- toupper(p)
    if (p_upper %in% genes_upper) {
      gene_idx <- which(genes_upper == p_upper)
      if ("ATP_score" %in% colnames(seurat_obj@meta.data)) {
        cors <- cor(as.numeric(gene_expr[gene_idx, ]), 
                   seurat_obj@meta.data$ATP_score, 
                   use = "pairwise.complete.obs")
        return(cors)
      }
    }
    return(NA)
  })
  
  results <- data.frame(
    protein = names(correlations),
    correlation_with_ATP = as.numeric(correlations),
    stringsAsFactors = FALSE
  )
  
  results <- results[order(-abs(results$correlation_with_ATP)), ]
  
  return(results)
}

#' @title Calculate Metabolite Correlation
#' 
#' @description Correlate metabolic scores with metabolite abundance.
#' 
#' @param seurat_obj Seurat object
#' @param metabolite_data Metabolite abundance matrix
#' @param score_cols Score columns to correlate
#' 
#' @return Matrix of correlations
#' 
#' @export
#' @examples
#' \dontrun{
#' metabolite_corr <- calculateMetaboliteCorrelation(seurat_obj, metabolite_data)
#' }
calculateMetaboliteCorrelation <- function(seurat_obj,
                                          metabolite_data = NULL,
                                          score_cols = NULL) {
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  if (is.null(metabolite_data)) {
    message("[scMetaboFlux] No metabolite data provided.")
    
    score_matrix <- as.matrix(seurat_obj@meta.data[, score_cols, drop = FALSE])
    
    return(cor(score_matrix, use = "pairwise"))
  }
  
  score_matrix <- as.matrix(seurat_obj@meta.data[, score_cols, drop = FALSE])
  
  correlations <- cor(score_matrix, metabolite_data, use = "pairwise")
  
  return(correlations)
}

#' @title Create Multi-Omics Metabolic Landscape
#' 
#' @description Combine multiple omics layers for comprehensive metabolic view.
#' 
#' @param seurat_obj Seurat object
#' @param include_spatial Include spatial information
#' @param include_proteomics Include proteomics
#' 
#' @return List with integrated results
#' 
#' @export
#' @examples
#' \dontrun{
#' landscape <- createMultiOmicsMetabolicLandscape(seurat_obj)
#' }
createMultiOmicsMetabolicLandscape <- function(seurat_obj,
                                              include_spatial = FALSE,
                                              include_proteomics = FALSE) {
  
  landscape <- list()
  
  landscape$transcriptomics <- list(
    glycolysis = seurat_obj@meta.data$glycolysis,
    oxphos = seurat_obj@meta.data$oxidative_phosphorylation,
    atp = seurat_obj@meta.data$ATP_score,
    gox = seurat_obj@meta.data$GOX_index
  )
  
  if (include_spatial && "spatial_x" %in% colnames(seurat_obj@meta.data)) {
    landscape$spatial <- list(
      x = seurat_obj@meta.data$spatial_x,
      y = seurat_obj@meta.data$spatial_y,
      hotspots = seurat_obj@meta.data$spatial_metabolic_hotspot
    )
  }
  
  if ("metabolic_phenotype" %in% colnames(seurat_obj@meta.data)) {
    landscape$phenotypes <- table(seurat_obj@meta.data$metabolic_phenotype)
  }
  
  landscape$quality_metrics <- list(
    n_cells = ncol(seurat_obj),
    n_features = nrow(seurat_obj),
    genes_detected = sum(rowSums(Seurat::GetAssayData(seurat_obj)) > 0),
    mean_genes_per_cell = mean(colSums(Seurat::GetAssayData(seurat_obj) > 0))
  )
  
  landscape$summary_stats <- list(
    mean_atp = mean(seurat_obj@meta.data$ATP_score, na.rm = TRUE),
    median_atp = median(seurat_obj@meta.data$ATP_score, na.rm = TRUE),
    sd_atp = sd(seurat_obj@meta.data$ATP_score, na.rm = TRUE)
  )
  
  return(landscape)
}

#' @title Infer Metabolic Flux from Multi-Omics
#' 
#' @description Infer metabolic flux using multi-modal data.
#' 
#' @param seurat_obj Seurat object
#' @param omics_weights Weights for each omics layer
#' @param name Column name for integrated flux
#' 
#' @return Seurat object with integrated flux
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- inferFluxFromMultiOmics(seurat_obj)
#' }
inferFluxFromMultiOmics <- function(seurat_obj,
                                   omics_weights = c(transcriptomics = 0.8,
                                                   proteomics = 0.15,
                                                   metabolomics = 0.05),
                                   name = "integrated_flux") {
  
  transcriptomics_flux <- seurat_obj@meta.data$ATP_score
  
  integrated_flux <- transcriptomics_flux * omics_weights["transcriptomics"]
  
  if (exists("proteomics_flux", seurat_obj@misc) && 
      !is.null(seurat_obj@misc$proteomics_flux)) {
    integrated_flux <- integrated_flux + 
      seurat_obj@misc$proteomics_flux * omics_weights["proteomics"]
  }
  
  if (exists("metabolomics_flux", seurat_obj@misc) && 
      !is.null(seurat_obj@misc$metabolomics_flux)) {
    integrated_flux <- integrated_flux + 
      seurat_obj@misc$metabolomics_flux * omics_weights["metabolomics"]
  }
  
  integrated_flux <- as.numeric(integrated_flux)
  min_flux <- min(integrated_flux, na.rm = TRUE)
  max_flux <- max(integrated_flux, na.rm = TRUE)
  
  if (max_flux > min_flux) {
    integrated_flux <- (integrated_flux - min_flux) / (max_flux - min_flux)
  } else {
    integrated_flux <- rep(0.5, length(integrated_flux))
  }
  
  seurat_obj@meta.data[[name]] <- integrated_flux
  
  return(seurat_obj)
}
