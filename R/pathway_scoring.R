#' @title Pathway Activity Scoring Module
#' @description Compute pathway activity scores using multiple methods
#' @name pathway_scoring
#' @rdname pathway-scoring
#' @keywords pathway scoring gene set enrichment
NULL

#' @title Compute Pathway Activity Scores
#' 
#' @description Calculate metabolic pathway activity scores for each cell using
#' various enrichment methods including AUCell, GSVA, and mean expression.
#' 
#' @param seurat_obj Seurat object with expression data
#' @param gene_sets Named list of gene sets (pathway genes). If NULL, uses
#'   all predefined metabolicGeneSets.
#' @param method Scoring method: "AUCell", "GSVA", "mean", "weighted_mean", "single_sample"
#' @param aucell_nbin AUCell: Number of bins for AUC calculation (default 24)
#' @param gsva_k GSVA: Number of genes for kernel width (default 100)
#' @param gsva_tau GSVA: Tau parameter for gene set scoring (default 0.25)
#' @param parallel_cores Number of cores for parallel computation (default 1)
#' @param assay Assay name to use (default active assay)
#' @param slot Slot to use for expression data (default "data")
#' @param verbose Print progress messages
#' 
#' @return Seurat object with pathway scores added to meta.data
#'   Each score is named according to the pathway (e.g., "glycolysis", "oxidative_phosphorylation")
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- computePathwayScores(seurat_obj, 
#'   gene_sets = metabolicGeneSets[c("glycolysis", "tca_cycle", 
#'   "oxidative_phosphorylation")],
#'   method = "AUCell")
#' }
#' 
#' \dontrun{
#' seurat_obj <- computePathwayScores(seurat_obj, 
#'   method = "mean",
#'   verbose = TRUE)
#' }
computePathwayScores <- function(seurat_obj,
                                 gene_sets = NULL,
                                 method = c("AUCell", "GSVA", "mean", "weighted_mean", "single_sample"),
                                 aucell_nbin = 24,
                                 gsva_k = 100,
                                 gsva_tau = 0.25,
                                 parallel_cores = 1,
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
  
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list")
  }
  
  if (verbose) {
    message("[scMetaboFlux] Computing pathway scores using '", method, "' method...")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  if (slot == "data") {
    expr_matrix <- as(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "data"), "dgCMatrix")
  } else if (slot == "counts") {
    expr_matrix <- as(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "counts"), "dgCMatrix")
  } else if (slot == "scale.data") {
    expr_matrix <- as(Seurat::GetAssayData(seurat_obj, assay = assay, layer = "scale.data"), "dgCMatrix")
  } else {
    expr_matrix <- as(Seurat::GetAssayData(seurat_obj, assay = assay), "dgCMatrix")
  }
  
  if (ncol(expr_matrix) == 0 || nrow(expr_matrix) == 0) {
    stop("Expression matrix is empty")
  }
  
  genes_rownames <- rownames(expr_matrix)
  cells_colnames <- colnames(expr_matrix)
  
  switch(method,
    "AUCell" = {
      scores <- computeAUCellScores(expr_matrix, gene_sets, aucell_nbin, 
                                   parallel_cores, verbose, genes_rownames)
    },
    "GSVA" = {
      scores <- computeGSVAScores(expr_matrix, gene_sets, gsva_k, gsva_tau,
                                 parallel_cores, verbose, genes_rownames)
    },
    "mean" = {
      scores <- computeMeanScores(expr_matrix, gene_sets, verbose, genes_rownames)
    },
    "weighted_mean" = {
      scores <- computeWeightedMeanScores(expr_matrix, gene_sets, verbose, genes_rownames)
    },
    "single_sample" = {
      scores <- computeSingleSampleScores(expr_matrix, gene_sets, verbose, genes_rownames)
    }
  )
  
  if (is.null(scores) || nrow(scores) == 0) {
    stop("Scoring failed to produce results")
  }
  
  colnames(scores) <- cells_colnames
  
  score_names <- rownames(scores)
  for (i in seq_along(score_names)) {
    score_name <- score_names[i]
    if (score_name %in% colnames(seurat_obj@meta.data)) {
      score_name <- paste0(score_name, "_", method)
    }
    seurat_obj@meta.data[[score_name]] <- as.numeric(scores[i, ])
  }
  
  if (verbose) {
    message("[scMetaboFlux] Pathway scores computed successfully!")
    message("[scMetaboFlux] Added ", length(score_names), " scores to meta.data")
  }
  
  return(seurat_obj)
}

#' @title AUCell Scoring
#' @description Internal function for AUCell-based scoring
#' @keywords internal
computeAUCellScores <- function(expr_matrix, gene_sets, nbin, nCores, verbose, genes_rownames) {
  tryCatch({
    if (verbose) message("[scMetaboFlux] Building AUCell rankings...")
    
    rankings <- AUCell::AUCell_buildRankings(
      expr_matrix, 
      nCores = nCores,
      plotStats = FALSE,
      verbose = FALSE
    )
    
    if (verbose) message("[scMetaboFlux] Calculating AUC scores...")
    
    gsc <- GSEABase::GeneSetCollection(lapply(names(gene_sets), function(name) {
      genes_in_set <- toupper(gene_sets[[name]])
      genes_in_data <- genes_rownames[toupper(genes_rownames) %in% genes_in_set]
      GSEABase::GeneSet(setName = name, geneIds = genes_in_data)
    }))
    
    if (length(gsc) == 0) {
      stop("No valid gene sets for AUCell analysis")
    }
    
    cells_AUC <- AUCell::AUCell_calcAUC(
      gsc, 
      rankings, 
      aucMaxRank = ceiling(0.05 * nrow(expr_matrix)),
      nCores = nCores, 
      verbose = FALSE
    )
    
    scores <- as.matrix(SummarizedExperiment::assay(cells_AUC))
    return(scores)
    
  }, error = function(e) {
    if (verbose) message("[scMetaboFlux] AUCell failed: ", e$message)
    if (verbose) message("[scMetaboFlux] Falling back to mean scoring...")
    return(computeMeanScores(expr_matrix, gene_sets, verbose, genes_rownames))
  })
}

#' @title GSVA Scoring
#' @description Internal function for GSVA-based scoring
#' @keywords internal
computeGSVAScores <- function(expr_matrix, gene_sets, k, tau, nCores, verbose, genes_rownames) {
  tryCatch({
    if (verbose) message("[scMetaboFlux] Running GSVA analysis...")
    
    gene_sets_filtered <- lapply(gene_sets, function(gs) {
      gs_upper <- toupper(gs)
      genes_rownames_upper <- toupper(genes_rownames)
      gs[gs_upper %in% genes_rownames_upper]
    })
    
    gene_sets_filtered <- gene_sets_filtered[sapply(gene_sets_filtered, length) >= 2]
    
    if (length(gene_sets_filtered) == 0) {
      if (verbose) message("[scMetaboFlux] No valid gene sets for GSVA, using mean scoring...")
      return(computeMeanScores(expr_matrix, gene_sets, verbose, genes_rownames))
    }
    
    scores <- GSVA::gsva(
      expr_matrix, 
      gene_sets_filtered, 
      method = "gsva",
      k = k, 
      tau = tau, 
      mx.diff = TRUE,
      parallel.sz = nCores, 
      verbose = FALSE
    )
    
    return(as.matrix(scores))
    
  }, error = function(e) {
    if (verbose) message("[scMetaboFlux] GSVA failed: ", e$message)
    if (verbose) message("[scMetaboFlux] Falling back to mean scoring...")
    return(computeMeanScores(expr_matrix, gene_sets, verbose, genes_rownames))
  })
}

#' @title Mean Expression Scoring
#' @description Internal function for simple mean expression scoring
#' @keywords internal
computeMeanScores <- function(expr_matrix, gene_sets, verbose, genes_rownames) {
  if (verbose) message("[scMetaboFlux] Computing mean expression scores...")
  
  genes_upper <- toupper(genes_rownames)
  
  scores <- sapply(names(gene_sets), function(pathway) {
    pathway_genes <- toupper(gene_sets[[pathway]])
    genes_in_data_idx <- match(pathway_genes, genes_upper)
    genes_in_data_idx <- genes_in_data_idx[!is.na(genes_in_data_idx)]
    
    if (length(genes_in_data_idx) == 0) {
      return(rep(0, ncol(expr_matrix)))
    }
    
    gene_expr <- expr_matrix[genes_in_data_idx, , drop = FALSE]
    mean_expr <- Matrix::colMeans(gene_expr, na.rm = TRUE)
    return(mean_expr)
  }, simplify = TRUE)
  
  if (is.vector(scores)) {
    scores <- matrix(scores, nrow = 1, dimnames = list(names(gene_sets)[1], NULL))
  }
  
  return(t(scores))
}

#' @title Weighted Mean Scoring
#' @description Internal function for weighted mean expression scoring
#' @keywords internal
computeWeightedMeanScores <- function(expr_matrix, gene_sets, verbose, genes_rownames) {
  if (verbose) message("[scMetaboFlux] Computing weighted mean scores...")
  
  genes_upper <- toupper(genes_rownames)
  weights_all <- rateLimitingWeights
  weights_upper <- toupper(names(weights_all))
  
  scores <- sapply(names(gene_sets), function(pathway) {
    pathway_genes <- toupper(gene_sets[[pathway]])
    genes_in_data_idx <- match(pathway_genes, genes_upper)
    genes_in_data_idx <- genes_in_data_idx[!is.na(genes_in_data_idx)]
    
    if (length(genes_in_data_idx) == 0) {
      return(rep(0, ncol(expr_matrix)))
    }
    
    gene_names <- genes_rownames[genes_in_data_idx]
    gene_names_upper <- toupper(gene_names)
    
    weights <- sapply(gene_names_upper, function(g) {
      idx <- match(g, weights_upper)
      if (!is.na(idx)) weights_all[idx] else 1.0
    })
    
    gene_expr <- expr_matrix[genes_in_data_idx, , drop = FALSE]
    
    weighted_expr <- t(t(gene_expr) * weights)
    weighted_sum <- Matrix::colSums(weighted_expr, na.rm = TRUE)
    total_weight <- sum(weights)
    
    weighted_mean <- weighted_sum / total_weight
    return(weighted_mean)
  }, simplify = TRUE)
  
  if (is.vector(scores)) {
    scores <- matrix(scores, nrow = 1, dimnames = list(names(gene_sets)[1], NULL))
  }
  
  return(t(scores))
}

#' @title Single Sample Scoring
#' @description Internal function for single sample gene set enrichment
#' @keywords internal
computeSingleSampleScores <- function(expr_matrix, gene_sets, verbose, genes_rownames) {
  if (verbose) message("[scMetaboFlux] Computing single sample scores...")
  
  genes_upper <- toupper(genes_rownames)
  
  scores <- sapply(names(gene_sets), function(pathway) {
    pathway_genes <- toupper(gene_sets[[pathway]])
    genes_in_data_idx <- match(pathway_genes, genes_upper)
    genes_in_data_idx <- genes_in_data_idx[!is.na(genes_in_data_idx)]
    
    if (length(genes_in_data_idx) == 0) {
      return(rep(0, ncol(expr_matrix)))
    }
    
    gene_expr <- as.matrix(expr_matrix[genes_in_data_idx, , drop = FALSE])
    
    gene_means <- Matrix::rowMeans(gene_expr, na.rm = TRUE)
    gene_sds <- apply(gene_expr, 1, sd, na.rm = TRUE)
    gene_sds[gene_sds == 0] <- 1
    
    zscore_matrix <- (gene_expr - gene_means) / gene_sds
    
    cell_zscore <- Matrix::colMeans(zscore_matrix, na.rm = TRUE)
    
    return(cell_zscore)
  }, simplify = TRUE)
  
  if (is.vector(scores)) {
    scores <- matrix(scores, nrow = 1, dimnames = list(names(gene_sets)[1], NULL))
  }
  
  return(t(scores))
}

#' @title Compute Combined Pathway Score
#' 
#' @description Combine multiple pathway scores into a single metric using
#' weighted or unweighted combinations.
#' 
#' @param seurat_obj Seurat object with pathway scores in meta.data
#' @param score_cols Character vector of score column names to combine
#' @param weights Optional numeric weights for each score (must sum to 1 or same length as score_cols)
#' @param name Name for combined score column
#' @param method Combination method: "mean", "sum", "weighted_sum", "geometric_mean"
#' 
#' @return Seurat object with combined score added to meta.data
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- computeCombinedScore(seurat_obj, 
#'   score_cols = c("glycolysis", "tca_cycle", "oxidative_phosphorylation"),
#'   weights = c(0.3, 0.3, 0.4),
#'   method = "weighted_sum")
#' }
computeCombinedScore <- function(seurat_obj, score_cols, weights = NULL,
                                name = "combined_score",
                                method = c("mean", "sum", "weighted_sum", "geometric_mean")) {
  
  method <- match.arg(method)
  
  missing_cols <- setdiff(score_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop(paste("Score columns not found in meta.data:", paste(missing_cols, collapse = ", ")))
  }
  
  score_matrix <- as.matrix(seurat_obj@meta.data[, score_cols, drop = FALSE])
  
  if (any(is.na(score_matrix))) {
    score_matrix[is.na(score_matrix)] <- 0
  }
  
  combined <- switch(method,
    "mean" = rowMeans(score_matrix, na.rm = TRUE),
    "sum" = rowSums(score_matrix, na.rm = TRUE),
    "weighted_sum" = {
      if (is.null(weights)) {
        weights <- rep(1, length(score_cols))
      }
      weights <- weights / sum(weights)
      rowSums(t(t(score_matrix) * weights), na.rm = TRUE)
    },
    "geometric_mean" = {
      score_matrix[score_matrix <= 0] <- 1e-10
      exp(rowMeans(log(score_matrix), na.rm = TRUE))
    }
  )
  
  seurat_obj@meta.data[[name]] <- as.numeric(combined)
  
  return(seurat_obj)
}

#' @title Z-Score Pathway Scores
#' 
#' @description Z-score normalize pathway scores across cells for better
#' comparison and to remove batch effects.
#' 
#' @param seurat_obj Seurat object with pathway scores
#' @param score_cols Columns to z-score. If NULL, auto-detects metabolic scores.
#' @param prefix Prefix for z-score columns (default "zscore_")
#' 
#' @return Seurat object with z-scored columns added
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- zscorePathwayScores(seurat_obj)
#' }
zscorePathwayScores <- function(seurat_obj, score_cols = NULL,
                                prefix = "zscore_") {
  
  if (is.null(score_cols)) {
    all_pathways <- c(names(metabolicGeneSets), 
                     paste0(names(metabolicGeneSets), "_AUCell"),
                     paste0(names(metabolicGeneSets), "_GSVA"),
                     paste0(names(metabolicGeneSets), "_mean"))
    score_cols <- intersect(all_pathways, colnames(seurat_obj@meta.data))
  }
  
  for (col in score_cols) {
    zscore_col <- paste0(prefix, col)
    values <- seurat_obj@meta.data[[col]]
    
    if (is.numeric(values)) {
      mean_val <- mean(values, na.rm = TRUE)
      sd_val <- sd(values, na.rm = TRUE)
      
      if (sd_val > 0) {
        seurat_obj@meta.data[[zscore_col]] <- (values - mean_val) / sd_val
      } else {
        seurat_obj@meta.data[[zscore_col]] <- 0
      }
    }
  }
  
  return(seurat_obj)
}

#' @title Get Pathway Score Matrix
#' 
#' @description Extract pathway scores as a matrix for downstream analysis.
#' 
#' @param seurat_obj Seurat object
#' @param score_pattern Pattern to match score columns (regex). If NULL,
#'   extracts all metabolic score columns.
#' 
#' @return Matrix of pathway scores (rows = scores, columns = cells)
#' 
#' @export
#' @examples
#' \dontrun{
#' score_matrix <- getPathwayScoreMatrix(seurat_obj)
#' dim(score_matrix)
#' }
getPathwayScoreMatrix <- function(seurat_obj, score_pattern = NULL) {
  
  if (is.null(score_pattern)) {
    all_pathways <- c(
      names(metabolicGeneSets),
      paste0(names(metabolicGeneSets), "_AUCell"),
      paste0(names(metabolicGeneSets), "_GSVA"),
      paste0(names(metabolicGeneSets), "_mean")
    )
    score_cols <- intersect(all_pathways, colnames(seurat_obj@meta.data))
  } else {
    score_cols <- grep(score_pattern, colnames(seurat_obj@meta.data), value = TRUE)
  }
  
  if (length(score_cols) == 0) {
    stop("No matching score columns found")
  }
  
  score_matrix <- as.matrix(seurat_obj@meta.data[, score_cols, drop = FALSE])
  rownames(score_matrix) <- score_cols
  
  return(score_matrix)
}

#' @title Compare Pathway Methods
#' 
#' @description Compare different pathway scoring methods by computing
#' pairwise correlations between methods.
#' 
#' @param seurat_obj Seurat object
#' @param gene_sets Gene sets to score
#' @param methods Character vector of methods to compare
#' @param assay Assay to use
#' 
#' @return Data frame with pairwise correlations between methods
#' 
#' @export
#' @examples
#' \dontrun{
#' comparison <- comparePathwayMethods(seurat_obj,
#'   gene_sets = metabolicGeneSets[1:3],
#'   methods = c("AUCell", "GSVA", "mean"))
#' head(comparison)
#' }
comparePathwayMethods <- function(seurat_obj, 
                                 gene_sets = metabolicGeneSets[1:3],
                                 methods = c("AUCell", "GSVA", "mean"),
                                 assay = NULL) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  results <- lapply(methods, function(m) {
    temp_obj <- computePathwayScores(
      seurat_obj, 
      gene_sets = gene_sets, 
      method = m, 
      assay = assay,
      verbose = FALSE
    )
    getPathwayScoreMatrix(temp_obj)
  })
  names(results) <- methods
  
  method_names <- names(results)
  correlations <- list()
  
  for (i in 1:(length(method_names) - 1)) {
    for (j in (i + 1):length(method_names)) {
      m1 <- method_names[i]
      m2 <- method_names[j]
      
      if (nrow(results[[m1]]) != nrow(results[[m2]])) next
      
      for (k in 1:nrow(results[[m1]])) {
        vals1 <- results[[m1]][k, ]
        vals2 <- results[[m2]][k, ]
        
        if (sd(vals1) > 0 && sd(vals2) > 0) {
          cors <- cor(vals1, vals2, use = "pairwise.complete.obs", method = "pearson")
        } else {
          cors <- NA
        }
        
        correlations[[length(correlations) + 1]] <- data.frame(
          method1 = m1,
          method2 = m2,
          pathway = rownames(results[[m1]])[k],
          correlation = cors,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  corr_df <- do.call(rbind, correlations)
  return(corr_df)
}
