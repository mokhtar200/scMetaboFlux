#' @title Metabolic Flux Approximation Module
#' @description Compute pseudo-flux scores for metabolic pathways
#' @name flux_approximation
#' @rdname flux-approximation
#' @keywords flux metabolic modeling
NULL

#' @title Compute Metabolic Flux Scores
#' 
#' @description Calculate pseudo-flux scores based on pathway gene expression.
#' 
#' The flux approximation formula is:
#' \code{Flux_pathway = sum(Expression_genes) / pathway_length}
#' 
#' Enhanced with rate-limiting enzyme weighting and normalization for 
#' improved biological relevance.
#' 
#' @param seurat_obj Seurat object with expression data
#' @param gene_sets Named list of pathway gene sets. If NULL, uses all metabolicGeneSets
#' @param use_rate_limiting Weight rate-limiting enzymes (default TRUE)
#' @param normalize Normalize flux scores across cells (default TRUE)
#' @param method Normalization method: "zscore", "minmax", "rank", or "none"
#' @param assay Assay to use
#' @param slot Slot for expression data
#' @param verbose Print progress messages
#' 
#' @return Seurat object with flux scores added to meta.data:
#'   \item{pathway_flux}{Raw flux score}
#'   \item{zscore_pathway_flux}{Z-score normalized (if normalize=TRUE)}
#'   \item{norm_pathway_flux}{Min-max normalized (if normalize=TRUE)}
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- computeMetabolicFlux(seurat_obj)
#' }
computeMetabolicFlux <- function(seurat_obj,
                                 gene_sets = NULL,
                                 use_rate_limiting = TRUE,
                                 normalize = TRUE,
                                 method = c("zscore", "minmax", "rank", "none"),
                                 assay = NULL,
                                 slot = "data",
                                 verbose = TRUE) {
  
  method <- match.arg(method)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (is.null(gene_sets)) {
    gene_sets <- metabolicGeneSets
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  if (verbose) {
    message("[scMetaboFlux] Computing metabolic flux scores...")
  }
  
  if (slot == "data") {
    expr_data <- as(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data"), "dgCMatrix")
  } else if (slot == "counts") {
    expr_data <- as(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts"), "dgCMatrix")
  } else if (slot == "scale.data") {
    expr_data <- as(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "scale.data"), "dgCMatrix")
  } else {
    expr_data <- as(Seurat::GetAssayData(seurat_obj, assay = assay), "dgCMatrix")
  }
  
  if (nrow(expr_data) == 0 || ncol(expr_data) == 0) {
    stop("Expression matrix is empty")
  }
  
  genes <- rownames(expr_data)
  genes_upper <- toupper(genes)
  
  flux_scores <- sapply(names(gene_sets), function(pathway) {
    pathway_genes <- toupper(gene_sets[[pathway]])
    
    matched_idx <- match(pathway_genes, genes_upper)
    matched_idx <- matched_idx[!is.na(matched_idx)]
    
    if (length(matched_idx) == 0) {
      return(rep(0, ncol(expr_data)))
    }
    
    gene_expr <- expr_data[matched_idx, , drop = FALSE]
    
    if (use_rate_limiting) {
      matched_genes <- genes[matched_idx]
      matched_genes_upper <- toupper(matched_genes)
      
      weights_upper <- toupper(names(rateLimitingWeights))
      
      weights <- sapply(matched_genes_upper, function(g) {
        idx <- match(g, weights_upper)
        if (!is.na(idx)) rateLimitingWeights[idx] else 1.0
      })
      
      if (length(weights) == 1) {
        weighted_expr <- gene_expr * weights
      } else {
        weighted_expr <- t(t(gene_expr) * weights)
      }
      
      flux <- Matrix::colSums(weighted_expr, na.rm = TRUE) / length(matched_idx)
    } else {
      flux <- Matrix::colSums(gene_expr, na.rm = TRUE) / length(matched_idx)
    }
    
    return(flux)
  }, simplify = TRUE)
  
  if (is.vector(flux_scores)) {
    flux_scores <- matrix(flux_scores, nrow = 1, 
                         dimnames = list(names(gene_sets)[1], NULL))
  }
  
  colnames(flux_scores) <- paste0(colnames(flux_scores), "_flux")
  
  for (col in colnames(flux_scores)) {
    seurat_obj@meta.data[[col]] <- as.numeric(flux_scores[, col])
  }
  
  if (normalize && method != "none") {
    flux_cols <- paste0(names(gene_sets), "_flux")
    flux_cols <- flux_cols[flux_cols %in% colnames(seurat_obj@meta.data)]
    
    for (col in flux_cols) {
      values <- seurat_obj@meta.data[[col]]
      
      if (method == "zscore") {
        mean_val <- mean(values, na.rm = TRUE)
        sd_val <- sd(values, na.rm = TRUE)
        if (sd_val > 0) {
          seurat_obj@meta.data[[paste0("zscore_", col)]] <- (values - mean_val) / sd_val
        } else {
          seurat_obj@meta.data[[paste0("zscore_", col)]] <- 0
        }
      } else if (method == "minmax") {
        min_val <- min(values, na.rm = TRUE)
        max_val <- max(values, na.rm = TRUE)
        if (max_val > min_val) {
          seurat_obj@meta.data[[paste0("norm_", col)]] <- (values - min_val) / (max_val - min_val)
        } else {
          seurat_obj@meta.data[[paste0("norm_", col)]] <- 0.5
        }
      } else if (method == "rank") {
        seurat_obj@meta.data[[paste0("rank_", col)]] <- 
          rank(values, ties.method = "average", na.last = "keep") / sum(!is.na(values))
      }
    }
  }
  
  if (verbose) {
    message("[scMetaboFlux] Flux scores computed for ", length(gene_sets), " pathways")
  }
  
  return(seurat_obj)
}

#' @title Calculate Flux Coupling
#' 
#' @description Calculate correlation between pathway fluxes to identify
#' coupled metabolic processes.
#' 
#' @param seurat_obj Seurat object with flux scores
#' @param flux_cols Columns with flux scores. If NULL, auto-detects.
#' @param method Correlation method: "pearson" or "spearman"
#' 
#' @return Matrix of flux correlations
#' 
#' @export
#' @examples
#' \dontrun{
#' cors <- calculateFluxCoupling(seurat_obj)
#' }
calculateFluxCoupling <- function(seurat_obj,
                                  flux_cols = NULL,
                                  method = c("pearson", "spearman")) {
  
  method <- match.arg(method)
  
  if (is.null(flux_cols)) {
    flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
    flux_cols <- flux_cols[!grepl("^(zscore|norm|rank)_", flux_cols)]
  }
  
  flux_cols <- intersect(flux_cols, colnames(seurat_obj@meta.data))
  
  if (length(flux_cols) < 2) {
    stop("Need at least 2 flux columns for correlation")
  }
  
  flux_matrix <- as.matrix(seurat_obj@meta.data[, flux_cols, drop = FALSE])
  flux_matrix[is.na(flux_matrix)] <- 0
  
  cors <- cor(flux_matrix, method = method, use = "pairwise")
  
  return(cors)
}

#' @title Calculate Metabolic Flux Balance
#' 
#' @description Calculate balance score between catabolic and anabolic fluxes.
#' 
#' @param seurat_obj Seurat object
#' @param catabolic_cols Catabolic pathway flux columns
#' @param anabolic_cols Anabolic pathway flux columns
#' @param name Column name for balance score
#' 
#' @return Seurat object with balance score
#' 
#' @export
calculateFluxBalance <- function(seurat_obj,
                                 catabolic_cols = c("glycolysis_flux", 
                                                    "oxidative_phosphorylation_flux"),
                                 anabolic_cols = c("tca_cycle_flux"),
                                 name = "flux_balance") {
  
  catabolic_cols <- intersect(catabolic_cols, colnames(seurat_obj@meta.data))
  anabolic_cols <- intersect(anabolic_cols, colnames(seurat_obj@meta.data))
  
  if (length(catabolic_cols) == 0) {
    stop("No catabolic columns found")
  }
  
  catabolic <- rowMeans(as.matrix(seurat_obj@meta.data[, catabolic_cols, drop = FALSE]), 
                        na.rm = TRUE)
  
  if (length(anabolic_cols) > 0) {
    anabolic <- rowMeans(as.matrix(seurat_obj@meta.data[, anabolic_cols, drop = FALSE]), 
                         na.rm = TRUE)
  } else {
    anabolic <- 0
  }
  
  balance <- (catabolic - anabolic) / (catabolic + anabolic + 1e-10)
  balance[is.infinite(balance) | is.nan(balance)] <- 0
  
  seurat_obj@meta.data[[name]] <- as.numeric(balance)
  
  return(seurat_obj)
}

#' @title Estimate Metabolic Flux Direction
#' 
#' @description Estimate predominant flux direction (glycolytic vs oxidative).
#' 
#' @param seurat_obj Seurat object
#' @param glycolysis_col Glycolysis flux column
#' @param oxphos_col OXPHOS flux column
#' @param name Column name
#' 
#' @return Seurat object with flux direction
#' 
#' @export
estimateFluxDirection <- function(seurat_obj,
                                 glycolysis_col = "glycolysis_flux",
                                 oxphos_col = "oxidative_phosphorylation_flux",
                                 name = "flux_direction") {
  
  glycolysis_flux <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos_flux <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  glycolysis_flux[is.na(glycolysis_flux)] <- 0
  oxphos_flux[is.na(oxphos_flux)] <- 0
  
  ratio <- glycolysis_flux / (oxphos_flux + 1e-10)
  ratio[is.infinite(ratio)] <- NA
  
  seurat_obj@meta.data[[name]] <- ratio
  
  seurat_obj@meta.data[[paste0(name, "_type")]] <- cut(
    ratio,
    breaks = c(-Inf, 0.5, 2, Inf),
    labels = c("OXPHOS_dominant", "Balanced", "Glycolysis_dominant")
  )
  
  return(seurat_obj)
}

#' @title Calculate Total Metabolic Flux
#' 
#' @description Calculate sum of all pathway fluxes.
#' 
#' @param seurat_obj Seurat object
#' @param flux_cols Columns to sum. If NULL, auto-detects.
#' @param name Column name
#' 
#' @return Seurat object with total flux
#' 
#' @export
calculateTotalMetabolicFlux <- function(seurat_obj,
                                       flux_cols = NULL,
                                       name = "total_flux") {
  
  if (is.null(flux_cols)) {
    flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
    flux_cols <- flux_cols[!grepl("^(zscore|norm|rank)_", flux_cols)]
  }
  
  flux_cols <- intersect(flux_cols, colnames(seurat_obj@meta.data))
  
  if (length(flux_cols) == 0) {
    stop("No flux columns found")
  }
  
  flux_matrix <- as.matrix(seurat_obj@meta.data[, flux_cols, drop = FALSE])
  total_flux <- rowSums(flux_matrix, na.rm = TRUE)
  
  seurat_obj@meta.data[[name]] <- as.numeric(total_flux)
  
  return(seurat_obj)
}

#' @title Calculate Flux Entropy
#' 
#' @description Calculate entropy of flux distribution across pathways.
#' Higher entropy indicates more balanced metabolic activity across pathways.
#' 
#' @param seurat_obj Seurat object
#' @param flux_cols Flux columns. If NULL, auto-detects.
#' @param name Column name
#' 
#' @return Seurat object with entropy
#' 
#' @export
calculateFluxEntropy <- function(seurat_obj,
                                flux_cols = NULL,
                                name = "flux_entropy") {
  
  if (is.null(flux_cols)) {
    flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
    flux_cols <- flux_cols[!grepl("^(zscore|norm|rank)_", flux_cols)]
  }
  
  flux_cols <- intersect(flux_cols, colnames(seurat_obj@meta.data))
  
  if (length(flux_cols) < 2) {
    stop("Need at least 2 flux columns for entropy calculation")
  }
  
  flux_matrix <- abs(as.matrix(seurat_obj@meta.data[, flux_cols, drop = FALSE]))
  flux_sum <- rowSums(flux_matrix, na.rm = TRUE)
  
  entropy <- apply(flux_matrix, 1, function(row) {
    p <- row / (sum(row) + 1e-10)
    p <- p[p > 1e-10]
    -sum(p * log(p))
  })
  
  entropy[is.na(entropy)] <- 0
  
  min_entropy <- min(entropy, na.rm = TRUE)
  max_entropy <- max(entropy, na.rm = TRUE)
  
  if (max_entropy > min_entropy) {
    entropy_norm <- (entropy - min_entropy) / (max_entropy - min_entropy)
  } else {
    entropy_norm <- rep(0.5, length(entropy))
  }
  
  seurat_obj@meta.data[[name]] <- as.numeric(entropy_norm)
  
  return(seurat_obj)
}

#' @title Detect Flux Switching Points
#' 
#' @description Identify cells undergoing metabolic switching between
#' glycolysis and OXPHOS.
#' 
#' @param seurat_obj Seurat object
#' @param glycolysis_col Glycolysis score
#' @param oxphos_col OXPHOS score
#' @param threshold Threshold for switching detection
#' @param window_size Sliding window size for smoothing
#' 
#' @return Data frame with switching information
#' 
#' @export
detectFluxSwitchingPoints <- function(seurat_obj,
                                     glycolysis_col = "glycolysis",
                                     oxphos_col = "oxidative_phosphorylation",
                                     threshold = 0.1,
                                     window_size = 50) {
  
  glycolysis <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  glycolysis[is.na(glycolysis)] <- 0
  oxphos[is.na(oxphos)] <- 0
  
  ratio <- glycolysis / (oxphos + 1e-10)
  ratio[is.infinite(ratio) | is.nan(ratio)] <- median(ratio, na.rm = TRUE)
  
  if (requireNamespace("pracma", quietly = TRUE)) {
    ratio_smooth <- pracma::movavg(ratio, window_size, type = "s")
    derivative <- diff(ratio_smooth)
  } else {
    ratio_smooth <- ratio
    derivative <- diff(ratio)
  }
  
  switching_idx <- which(abs(derivative) > threshold)
  
  if (length(switching_idx) == 0) {
    return(data.frame(
      cell_id = character(0),
      position = integer(0),
      ratio = numeric(0),
      derivative = numeric(0),
      direction = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  results <- data.frame(
    cell_id = rownames(seurat_obj@meta.data)[switching_idx],
    position = switching_idx,
    ratio = ratio[switching_idx],
    derivative = derivative[switching_idx],
    direction = ifelse(derivative[switching_idx] > 0, "Glycolytic", "Oxidative"),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' @title Calculate Pathway Flux Contribution
#' 
#' @description Calculate percentage contribution of each pathway to total flux.
#' 
#' @param seurat_obj Seurat object
#' @param flux_cols Flux columns. If NULL, auto-detects.
#' @param prefix Prefix for contribution columns
#' 
#' @return Seurat object with contribution percentages
#' 
#' @export
calculateFluxContribution <- function(seurat_obj,
                                     flux_cols = NULL,
                                     prefix = "contribution_") {
  
  if (is.null(flux_cols)) {
    flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
    flux_cols <- flux_cols[!grepl("^(zscore|norm|rank)_", flux_cols)]
  }
  
  flux_cols <- intersect(flux_cols, colnames(seurat_obj@meta.data))
  
  if (length(flux_cols) == 0) {
    stop("No flux columns found")
  }
  
  flux_matrix <- as.matrix(seurat_obj@meta.data[, flux_cols, drop = FALSE])
  total_flux <- rowSums(flux_matrix, na.rm = TRUE)
  total_flux[total_flux == 0] <- 1
  
  for (col in flux_cols) {
    pathway_name <- gsub("_flux$", "", col)
    contribution_col <- paste0(prefix, pathway_name)
    
    contribution <- (as.numeric(seurat_obj@meta.data[[col]]) / total_flux) * 100
    contribution[is.na(contribution) | is.infinite(contribution)] <- 0
    
    seurat_obj@meta.data[[contribution_col]] <- contribution
  }
  
  return(seurat_obj)
}

#' @title Simulate Flux Distribution
#' 
#' @description Simulate metabolic flux distribution with optional noise.
#' 
#' @param seurat_obj Seurat object
#' @param condition_col Condition column
#' @param conditions Conditions to compare
#' @param n_simulations Number of simulations
#' @param noise_sd Standard deviation of noise to add
#' 
#' @return List with simulated flux distributions
#' 
#' @export
simulateFluxDistribution <- function(seurat_obj,
                                     condition_col = NULL,
                                     conditions = NULL,
                                     n_simulations = 100,
                                     noise_sd = 0.1) {
  
  flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
  flux_cols <- flux_cols[!grepl("^(zscore|norm|rank)_", flux_cols)]
  flux_cols <- intersect(flux_cols, colnames(seurat_obj@meta.data))
  
  if (is.null(condition_col)) {
    flux_matrix <- as.matrix(seurat_obj@meta.data[, flux_cols, drop = FALSE])
    flux_matrix[is.na(flux_matrix)] <- 0
    
    sim_results <- lapply(1:n_simulations, function(i) {
      noise <- matrix(rnorm(length(flux_matrix), 0, noise_sd), 
                      nrow = nrow(flux_matrix), ncol = ncol(flux_matrix))
      flux_matrix + noise
    })
    
    return(list(
      mean_flux = colMeans(flux_matrix),
      simulated_flux = Reduce(`+`, sim_results) / n_simulations,
      sd_flux = apply(flux_matrix, 2, sd),
      n_simulations = n_simulations
    ))
  }
  
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop("Condition column not found")
  }
  
  if (is.null(conditions)) {
    conditions <- unique(seurat_obj@meta.data[[condition_col]])
  }
  
  results <- lapply(conditions, function(cond) {
    cond_cells <- seurat_obj@meta.data[[condition_col]] == cond
    flux_matrix <- as.matrix(seurat_obj@meta.data[cond_cells, flux_cols, drop = FALSE])
    flux_matrix[is.na(flux_matrix)] <- 0
    
    sim_results <- lapply(1:n_simulations, function(i) {
      noise <- matrix(rnorm(length(flux_matrix), 0, noise_sd), 
                      nrow = nrow(flux_matrix), ncol = ncol(flux_matrix))
      flux_matrix + noise
    })
    
    list(
      condition = cond,
      mean_flux = colMeans(flux_matrix),
      simulated_flux = Reduce(`+`, sim_results) / n_simulations,
      sd_flux = apply(flux_matrix, 2, sd),
      n_cells = sum(cond_cells)
    )
  })
  
  names(results) <- conditions
  
  return(results)
}
