#' @title Metabolic Phenotype Classification Module
#' @description Classify cells into metabolic phenotypes
#' @name phenotype_classification
#' @rdname phenotype-classification
#' @keywords phenotype classification metabolism
NULL

#' @title Classify Metabolic Phenotypes
#' 
#' @description Assign metabolic phenotype to each cell based on pathway activity.
#' 
#' Phenotypes:
#' \itemize{
#'   \item Glycolytic: High glycolysis, low OXPHOS (Warburg-like)
#'   \item Oxidative: High OXPHOS, low glycolysis
#'   \item Energetically Balanced: Both moderate
#'   \item Energy-Stressed: Low total ATP
#'   \item Hypermetabolic: High total ATP
#' }
#' 
#' @param seurat_obj Seurat object with pathway and ATP scores
#' @param glycolysis_col Glycolysis score column
#' @param oxphos_col OXPHOS score column
#' @param atp_col ATP score column. If NULL, estimates from glycolysis/oxphos
#' @param glycolysis_threshold Quantile thresholds for glycolysis (c(low, high))
#' @param oxphos_threshold Quantile thresholds for OXPHOS (c(low, high))
#' @param atp_low_threshold Quantile threshold for low ATP (default 0.25)
#' @param atp_high_threshold Quantile threshold for high ATP (default 0.75)
#' @param method Classification method: "quantile", "kmeans", or "hierarchical"
#' @param name Column name for phenotype
#' @param verbose Print progress messages
#' 
#' @return Seurat object with phenotype classification added to meta.data
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- classifyMetabolicPhenotype(seurat_obj)
#' table(seurat_obj@meta.data$metabolic_phenotype)
#' }
classifyMetabolicPhenotype <- function(seurat_obj,
                                       glycolysis_col = "glycolysis",
                                       oxphos_col = "oxidative_phosphorylation",
                                       atp_col = NULL,
                                       glycolysis_threshold = c(0.33, 0.66),
                                       oxphos_threshold = c(0.33, 0.66),
                                       atp_low_threshold = 0.25,
                                       atp_high_threshold = 0.75,
                                       method = c("quantile", "kmeans", "hierarchical"),
                                       name = "metabolic_phenotype",
                                       verbose = TRUE) {
  
  method <- match.arg(method)
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (!glycolysis_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", glycolysis_col, "not found in meta.data"))
  }
  
  if (!oxphos_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", oxphos_col, "not found in meta.data"))
  }
  
  if (verbose) {
    message("[scMetaboFlux] Classifying metabolic phenotypes...")
  }
  
  glycolysis <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  glycolysis[is.na(glycolysis)] <- median(glycolysis, na.rm = TRUE)
  oxphos[is.na(oxphos)] <- median(oxphos, na.rm = TRUE)
  
  if (is.null(atp_col)) {
    if ("ATP_score" %in% colnames(seurat_obj@meta.data)) {
      atp_col <- "ATP_score"
    } else {
      atp <- (glycolysis + oxphos) / 2
    }
  }
  
  if (!is.null(atp_col) && atp_col %in% colnames(seurat_obj@meta.data)) {
    atp <- as.numeric(seurat_obj@meta.data[[atp_col]])
    atp[is.na(atp)] <- median(atp, na.rm = TRUE)
  } else {
    atp <- (glycolysis + oxphos) / 2
  }
  
  phenotype <- switch(method,
    "quantile" = classifyByQuantile(glycolysis, oxphos, atp,
                                   glycolysis_threshold,
                                   oxphos_threshold,
                                   atp_low_threshold,
                                   atp_high_threshold),
    "kmeans" = classifyByKmeans(glycolysis, oxphos, atp),
    "hierarchical" = classifyByHierarchical(glycolysis, oxphos, atp)
  )
  
  phenotype_levels <- c("Glycolytic", "Oxidative", "Energetically Balanced",
                        "Energy-Stressed", "Hypermetabolic")
  
  phenotype <- factor(phenotype, levels = phenotype_levels)
  
  seurat_obj@meta.data[[name]] <- phenotype
  
  if (verbose) {
    message("[scMetaboFlux] Phenotype distribution:")
    print(table(phenotype, useNA = "ifany"))
  }
  
  return(seurat_obj)
}

#' @title Classify by Quantile Method
#' @keywords internal
classifyByQuantile <- function(glycolysis, oxphos, atp,
                              glycolysis_threshold,
                              oxphos_threshold,
                              atp_low_threshold,
                              atp_high_threshold) {
  
  gly_q <- quantile(glycolysis, probs = glycolysis_threshold, na.rm = TRUE)
  oxphos_q <- quantile(oxphos, probs = oxphos_threshold, na.rm = TRUE)
  atp_low_q <- quantile(atp, probs = atp_low_threshold, na.rm = TRUE)
  atp_high_q <- quantile(atp, probs = atp_high_threshold, na.rm = TRUE)
  
  phenotype <- rep("Energetically Balanced", length(glycolysis))
  
  high_gly_low_oxphos <- glycolysis >= gly_q[2] & oxphos <= oxphos_q[1]
  phenotype[high_gly_low_oxphos] <- "Glycolytic"
  
  low_gly_high_oxphos <- glycolysis <= gly_q[1] & oxphos >= oxphos_q[2]
  phenotype[low_gly_high_oxphos] <- "Oxidative"
  
  both_high <- glycolysis >= gly_q[2] & oxphos >= oxphos_q[2]
  phenotype[both_high] <- "Hypermetabolic"
  
  both_low <- glycolysis <= gly_q[1] & oxphos <= oxphos_q[1]
  phenotype[both_low] <- "Energy-Stressed"
  
  low_atp <- atp <= atp_low_q & !phenotype %in% c("Hypermetabolic")
  phenotype[low_atp] <- "Energy-Stressed"
  
  high_atp <- atp >= atp_high_q & !phenotype %in% c("Energy-Stressed")
  phenotype[high_atp] <- "Hypermetabolic"
  
  return(phenotype)
}

#' @title Classify by K-means
#' @keywords internal
classifyByKmeans <- function(glycolysis, oxphos, atp) {
  
  features <- cbind(glycolysis, oxphos, atp)
  features <- scale(features)
  
  n_cells <- nrow(features)
  n_centers <- min(5, n_cells)
  
  set.seed(42)
  kmeans_result <- kmeans(features, centers = n_centers, nstart = 100, 
                          iter.max = 100)
  
  cluster_phenotypes <- assignClusterPhenotypes(kmeans_result$centers)
  
  phenotype <- cluster_phenotypes[kmeans_result$cluster]
  
  return(phenotype)
}

#' @title Assign Cluster Phenotypes
#' @keywords internal
assignClusterPhenotypes <- function(centers) {
  phenotypes <- rep("Energetically Balanced", nrow(centers))
  
  if (ncol(centers) >= 3) {
    gox_ratio <- centers[, 1] - centers[, 2]
    atp_vals <- centers[, 3]
  } else {
    gox_ratio <- centers[, 1]
    atp_vals <- centers[, 2]
  }
  
  atp_median <- median(atp_vals)
  atp_75 <- quantile(atp_vals, 0.75)
  atp_25 <- quantile(atp_vals, 0.25)
  
  for (i in 1:nrow(centers)) {
    if (gox_ratio[i] > 0.5 && atp_vals[i] > atp_median) {
      phenotypes[i] <- "Glycolytic"
    } else if (gox_ratio[i] < -0.5 && atp_vals[i] > atp_median) {
      phenotypes[i] <- "Oxidative"
    } else if (gox_ratio[i] > 0.5 && atp_vals[i] <= atp_median) {
      phenotypes[i] <- "Energy-Stressed"
    } else if (gox_ratio[i] < -0.5 && atp_vals[i] <= atp_median) {
      phenotypes[i] <- "Energy-Stressed"
    } else {
      phenotypes[i] <- "Energetically Balanced"
    }
  }
  
  high_atp_clusters <- which(atp_vals >= atp_75)
  phenotypes[high_atp_clusters] <- "Hypermetabolic"
  
  low_atp_clusters <- which(atp_vals <= atp_25)
  phenotypes[low_atp_clusters] <- "Energy-Stressed"
  
  return(phenotypes)
}

#' @title Classify by Hierarchical Clustering
#' @keywords internal
classifyByHierarchical <- function(glycolysis, oxphos, atp) {
  
  features <- cbind(glycolysis, oxphos, atp)
  features <- scale(features)
  
  d <- dist(features, method = "euclidean")
  hc <- hclust(d, method = "ward.D2")
  
  n_clusters <- min(5, nrow(features) / 10)
  n_clusters <- max(3, n_clusters)
  
  clusters <- cutree(hc, k = n_clusters)
  
  cluster_means <- aggregate(features, by = list(clusters), mean)
  cluster_phenotypes <- assignClusterPhenotypes(as.matrix(cluster_means[, -1]))
  
  phenotype <- cluster_phenotypes[clusters]
  
  return(phenotype)
}

#' @title Classify Using GOX Index
#' 
#' @description Classify phenotypes using pre-computed GOX index.
#' 
#' @param seurat_obj Seurat object
#' @param gox_col GOX index column
#' @param atp_col ATP score column
#' @param gox_threshold Threshold for GOX classification (default 0)
#' @param atp_threshold Threshold for ATP classification
#' @param name Column name
#' 
#' @return Seurat object with phenotype classification
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- calculateGOXIndex(seurat_obj)
#' seurat_obj <- classifyByGOXIndex(seurat_obj)
#' }
classifyByGOXIndex <- function(seurat_obj,
                               gox_col = "GOX_index",
                               atp_col = "ATP_score",
                               gox_threshold = 0,
                               atp_threshold = NULL,
                               name = "metabolic_phenotype") {
  
  gox <- as.numeric(seurat_obj@meta.data[[gox_col]])
  gox[is.na(gox)] <- 0
  
  if (is.null(atp_threshold)) {
    atp_threshold <- median(seurat_obj@meta.data[[atp_col]], na.rm = TRUE)
  }
  
  atp <- as.numeric(seurat_obj@meta.data[[atp_col]])
  atp[is.na(atp)] <- median(atp, na.rm = TRUE)
  
  phenotype <- rep("Energetically Balanced", length(gox))
  
  glycolytic_mask <- gox > gox_threshold & atp > atp_threshold
  phenotype[glycolytic_mask] <- "Glycolytic"
  
  oxidative_mask <- gox < -gox_threshold & atp > atp_threshold
  phenotype[oxidative_mask] <- "Oxidative"
  
  stressed_mask <- atp <= atp_threshold
  phenotype[stressed_mask] <- "Energy-Stressed"
  
  hyper_mask <- gox > gox_threshold & atp > quantile(atp, 0.9, na.rm = TRUE)
  phenotype[hyper_mask] <- "Hypermetabolic"
  
  seurat_obj@meta.data[[name]] <- factor(phenotype,
                                         levels = c("Glycolytic", "Oxidative",
                                                    "Energetically Balanced",
                                                    "Energy-Stressed", "Hypermetabolic"))
  
  return(seurat_obj)
}

#' @title Get Phenotype Statistics
#' 
#' @description Get statistics for each metabolic phenotype.
#' 
#' @param seurat_obj Seurat object
#' @param phenotype_col Phenotype column
#' @param score_cols Score columns to include
#' 
#' @return List containing:
#'   \item{statistics}{Data frame with phenotype statistics}
#'   \item{cell_counts}{Data frame with cell counts per phenotype}
#' 
#' @export
#' @examples
#' \dontrun{
#' stats <- getPhenotypeStatistics(seurat_obj)
#' print(stats$cell_counts)
#' }
getPhenotypeStatistics <- function(seurat_obj,
                                  phenotype_col = "metabolic_phenotype",
                                  score_cols = NULL) {
  
  if (!phenotype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", phenotype_col, "not found in meta.data"))
  }
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score", "GOX_index")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  stats_list <- lapply(score_cols, function(col) {
    if (!col %in% colnames(seurat_obj@meta.data)) return(NULL)
    
    seurat_obj@meta.data %>%
      dplyr::group_by(.data[[phenotype_col]]) %>%
      dplyr::summarize(
        metric = col,
        mean = mean(.data[[col]], na.rm = TRUE),
        sd = sd(.data[[col]], na.rm = TRUE),
        median = median(.data[[col]], na.rm = TRUE),
        min = min(.data[[col]], na.rm = TRUE),
        max = max(.data[[col]], na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  stats_list <- stats_list[!sapply(stats_list, is.null)]
  stats_df <- do.call(rbind, stats_list)
  
  counts <- table(seurat_obj@meta.data[[phenotype_col]], useNA = "ifany")
  counts_df <- data.frame(
    phenotype = names(counts),
    n_cells = as.numeric(counts),
    pct = as.numeric(counts) / sum(counts) * 100,
    stringsAsFactors = FALSE
  )
  
  return(list(
    statistics = stats_df,
    cell_counts = counts_df
  ))
}

#' @title Calculate Phenotype Transition Probability
#' 
#' @description Calculate probability distribution across phenotypes.
#' 
#' @param seurat_obj Seurat object
#' @param phenotype_col Phenotype column
#' @param group_col Grouping column for transition matrix
#' 
#' @return List with transition probabilities
#' 
#' @export
calculatePhenotypeTransitionProbability <- function(seurat_obj,
                                                   phenotype_col = "metabolic_phenotype",
                                                   group_col = NULL) {
  
  phenotypes <- seurat_obj@meta.data[[phenotype_col]]
  phenotypes[is.na(phenotypes)] <- "Unknown"
  
  if (is.null(group_col)) {
    trans_matrix <- table(phenotypes) / length(phenotypes)
    return(list(
      transition_matrix = as.matrix(trans_matrix),
      by_group = NULL
    ))
  }
  
  groups <- seurat_obj@meta.data[[group_col]]
  unique_groups <- unique(groups[!is.na(groups)])
  
  by_group <- lapply(unique_groups, function(g) {
    g_mask <- groups == g & !is.na(groups)
    tab <- table(phenotypes[g_mask]) / sum(g_mask)
    list(group = g, matrix = as.matrix(tab), n_cells = sum(g_mask))
  })
  
  names(by_group) <- unique_groups
  
  return(list(
    transition_matrix = NULL,
    by_group = by_group
  ))
}

#' @title Get Phenotype Distribution
#' 
#' @description Get phenotype distribution summary.
#' 
#' @param seurat_obj Seurat object
#' @param phenotype_col Phenotype column
#' @param group_col Grouping column for comparison
#' 
#' @return Data frame with phenotype distribution
#' 
#' @export
#' @examples
#' \dontrun{
#' dist <- getPhenotypeDistribution(seurat_obj, group_col = "cell_type")
#' head(dist)
#' }
getPhenotypeDistribution <- function(seurat_obj,
                                    phenotype_col = "metabolic_phenotype",
                                    group_col = NULL) {
  
  if (!phenotype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", phenotype_col, "not found in meta.data"))
  }
  
  if (is.null(group_col)) {
    dist <- data.frame(
      phenotype = seurat_obj@meta.data[[phenotype_col]],
      stringsAsFactors = FALSE
    )
    
    dist <- dist %>%
      dplyr::group_by(phenotype) %>%
      dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
      dplyr::mutate(pct = n / sum(n) * 100)
    
    return(as.data.frame(dist))
  }
  
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Column", group_col, "not found in meta.data"))
  }
  
  dist <- seurat_obj@meta.data %>%
    dplyr::group_by(.data[[group_col]], .data[[phenotype_col]]) %>%
    dplyr::summarize(n = dplyr::n(), .groups = "drop") %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(pct = n / sum(n) * 100)
  
  return(as.data.frame(dist))
}

#' @title Calculate Phenotype Purity Score
#' 
#' @description Calculate purity of phenotype assignment within clusters.
#' 
#' @param seurat_obj Seurat object
#' @param phenotype_col Phenotype column
#' @param cluster_col Cluster column
#' 
#' @return Data frame with purity scores
#' 
#' @export
#' @examples
#' \dontrun{
#' purity <- calculatePhenotypePurity(seurat_obj)
#' print(purity)
#' }
calculatePhenotypePurity <- function(seurat_obj,
                                    phenotype_col = "metabolic_phenotype",
                                    cluster_col = "seurat_clusters") {
  
  if (!phenotype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Phenotype column", phenotype_col, "not found"))
  }
  
  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cluster column", cluster_col, "not found"))
  }
  
  contingency <- table(seurat_obj@meta.data[[cluster_col]],
                      seurat_obj@meta.data[[phenotype_col]])
  
  row_totals <- rowSums(contingency)
  max_per_row <- apply(contingency, 1, max)
  
  purity_df <- data.frame(
    cluster = rownames(contingency),
    purity = round(max_per_row / row_totals, 3),
    n_cells = as.numeric(row_totals),
    dominant_phenotype = apply(contingency, 1, function(x) names(which.max(x))),
    stringsAsFactors = FALSE
  )
  
  return(purity_df)
}

#' @title Refine Phenotype Classification
#' 
#' @description Refine phenotype using nearest neighbor smoothing for
#' more robust classification.
#' 
#' @param seurat_obj Seurat object
#' @param phenotype_col Phenotype column
#' @param k Number of neighbors for smoothing
#' @param name Column name for refined phenotype
#' @param reduction Reduction to use for neighbors
#' 
#' @return Seurat object with refined phenotype
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- refinePhenotypeClassification(seurat_obj, k = 10)
#' }
refinePhenotypeClassification <- function(seurat_obj,
                                         phenotype_col = "metabolic_phenotype",
                                         k = 10,
                                         name = "metabolic_phenotype_refined",
                                         reduction = "pca") {
  
  if (!phenotype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Phenotype column", phenotype_col, "not found"))
  }
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found"))
  }
  
  metadata <- seurat_obj@meta.data
  phenotype <- as.numeric(metadata[[phenotype_col]])
  
  tryCatch({
    snn <- Seurat::FindNeighbors(seurat_obj, reduction = reduction, dims = 1:50,
                                k.param = k, verbose = FALSE)$snn
    
    refined_phenotype <- sapply(1:nrow(metadata), function(i) {
      neighbors <- as.numeric(snn[i, 1:k])
      neighbor_phenotypes <- phenotype[neighbors]
      
      if (all(is.na(neighbor_phenotypes))) {
        return(phenotype[i])
      }
      
      neighbor_phenotypes <- neighbor_phenotypes[!is.na(neighbor_phenotypes)]
      
      modal_phenotype <- as.numeric(names(which.max(table(neighbor_phenotypes))))
      
      consensus_ratio <- sum(neighbor_phenotypes == phenotype[i]) / length(neighbor_phenotypes)
      
      if (consensus_ratio < 0.5) {
        return(modal_phenotype)
      }
      return(phenotype[i])
    })
    
  }, error = function(e) {
    message("[scMetaboFlux] Refinement failed, using original classification")
    refined_phenotype <- phenotype
  })
  
  phenotype_levels <- c("Glycolytic", "Oxidative", "Energetically Balanced",
                        "Energy-Stressed", "Hypermetabolic")
  
  seurat_obj@meta.data[[name]] <- factor(refined_phenotype, levels = 1:5,
                                         labels = phenotype_levels)
  
  return(seurat_obj)
}
