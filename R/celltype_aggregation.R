#' @title Cell-Type Level Aggregation Module
#' @description Aggregate metabolic scores at cell-type level
#' @name celltype_aggregation
#' @rdname celltype-aggregation
#' @keywords cell type aggregation metabolism
NULL

#' @title Aggregate Scores by Cell Type
#' 
#' @description Calculate metabolic statistics aggregated by cell type.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type annotation column
#' @param score_cols Score columns to aggregate. If NULL, uses default metabolic scores
#' @param include_sd Include standard deviation
#' @param include_sem Include standard error of mean
#' @param include_ci Include 95% confidence interval
#' @param verbose Print progress messages
#' 
#' @return Data frame with aggregated statistics per cell type
#' 
#' @export
#' @examples
#' \dontrun{
#' cell_type_stats <- aggregateByCellType(seurat_obj, 
#'   cell_type_col = "cell_type",
#'   score_cols = c("glycolysis", "oxidative_phosphorylation", "ATP_score"))
#' head(cell_type_stats)
#' }
aggregateByCellType <- function(seurat_obj,
                                cell_type_col = "cell_type",
                                score_cols = NULL,
                                include_sd = TRUE,
                                include_sem = TRUE,
                                include_ci = FALSE,
                                verbose = TRUE) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cell type column", cell_type_col, "not found in metadata"))
  }
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score", 
                     "GOX_index", "metabolic_power")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  if (length(score_cols) == 0) {
    stop("No valid score columns found")
  }
  
  stats_list <- lapply(score_cols, function(col) {
    if (!col %in% colnames(seurat_obj@meta.data)) return(NULL)
    
    col_data <- seurat_obj@meta.data[[col]]
    
    if (!is.numeric(col_data)) return(NULL)
    
    stats <- seurat_obj@meta.data %>%
      dplyr::group_by(.data[[cell_type_col]]) %>%
      dplyr::summarize(
        metric = col,
        mean = mean(.data[[col]], na.rm = TRUE),
        median = median(.data[[col]], na.rm = TRUE),
        .groups = "drop"
      )
    
    if (include_sd) {
      sd_stats <- seurat_obj@meta.data %>%
        dplyr::group_by(.data[[cell_type_col]]) %>%
        dplyr::summarize(
          sd = sd(.data[[col]], na.rm = TRUE),
          .groups = "drop"
        )
      stats <- merge(stats, sd_stats, by = cell_type_col)
    }
    
    if (include_sem || include_ci) {
      n_stats <- seurat_obj@meta.data %>%
        dplyr::group_by(.data[[cell_type_col]]) %>%
        dplyr::summarize(
          n = dplyr::n(),
          .groups = "drop"
        )
      stats <- merge(stats, n_stats, by = cell_type_col)
      
      if (include_sem) {
        stats$sem <- stats$sd / sqrt(stats$n)
      }
      
      if (include_ci) {
        stats$ci_lower <- stats$mean - 1.96 * stats$sem
        stats$ci_upper <- stats$mean + 1.96 * stats$sem
      }
    }
    
    return(stats)
  })
  
  stats_list <- stats_list[!sapply(stats_list, is.null)]
  
  if (length(stats_list) == 0) {
    stop("No valid statistics could be computed")
  }
  
  aggregated <- do.call(rbind, stats_list)
  
  cell_counts <- seurat_obj@meta.data %>%
    dplyr::count(.data[[cell_type_col]], name = "total_cells")
  
  colnames(cell_counts)[1] <- cell_type_col
  
  aggregated <- merge(aggregated, cell_counts, by = cell_type_col, all.x = TRUE)
  
  aggregated$mean <- round(aggregated$mean, 4)
  aggregated$median <- round(aggregated$median, 4)
  if (include_sd) aggregated$sd <- round(aggregated$sd, 4)
  
  return(aggregated)
}

#' @title Compare Cell Types by Metabolic Profile
#' 
#' @description Compare metabolic profiles between cell types using statistical tests.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param score_cols Score columns
#' @param method Statistical test: "t.test", "wilcox.test", "anova"
#' @param p.adjust.method P-value adjustment method
#' 
#' @return List with comparison results
#' 
#' @export
#' @examples
#' \dontrun{
#' results <- compareCellTypesMetabolism(seurat_obj, cell_type_col = "cell_type")
#' }
compareCellTypesMetabolism <- function(seurat_obj,
                                      cell_type_col = "cell_type",
                                      score_cols = NULL,
                                      method = c("t.test", "wilcox.test", "anova"),
                                      p.adjust.method = "BH") {
  
  method <- match.arg(method)
  
  if (!cell_type_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Cell type column not found"))
  }
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  cell_types <- unique(seurat_obj@meta.data[[cell_type_col]])
  
  if (length(cell_types) < 2) {
    stop("Need at least 2 cell types for comparison")
  }
  
  results <- lapply(score_cols, function(col) {
    if (!col %in% colnames(seurat_obj@meta.data)) return(NULL)
    
    if (method == "anova") {
      formula <- as.formula(paste0("`", col, "`", " ~ ", cell_type_col))
      
      anova_result <- aov(formula, data = seurat_obj@meta.data)
      anova_summary <- summary(anova_result)[[1]]
      
      tukey <- TukeyHSD(anova_result)
      
      return(list(
        test = "ANOVA",
        p_value = anova_summary$"Pr(>F)"[1],
        f_statistic = anova_summary$"F value"[1],
        tukey_hsd = tukey
      ))
    } else {
      comparisons <- combn(cell_types, 2, simplify = FALSE)
      
      pairwise_results <- lapply(comparisons, function(comp) {
        group1_data <- seurat_obj@meta.data[[col]][seurat_obj@meta.data[[cell_type_col]] == comp[1]]
        group2_data <- seurat_obj@meta.data[[col]][seurat_obj@meta.data[[cell_type_col]] == comp[2]]
        
        group1_data <- group1_data[!is.na(group1_data)]
        group2_data <- group2_data[!is.na(group2_data)]
        
        if (length(group1_data) < 2 || length(group2_data) < 2) {
          return(data.frame(
            group1 = comp[1],
            group2 = comp[2],
            mean_diff = NA,
            p_value = NA,
            adj_p_value = NA,
            stringsAsFactors = FALSE
          ))
        }
        
        if (method == "t.test") {
          test_result <- tryCatch({
            t.test(group1_data, group2_data)
          }, error = function(e) NULL)
        } else {
          test_result <- tryCatch({
            wilcox.test(group1_data, group2_data)
          }, error = function(e) NULL)
        }
        
        if (is.null(test_result)) {
          return(data.frame(
            group1 = comp[1],
            group2 = comp[2],
            mean_diff = mean(group2_data, na.rm = TRUE) - mean(group1_data, na.rm = TRUE),
            p_value = NA,
            adj_p_value = NA,
            stringsAsFactors = FALSE
          ))
        }
        
        data.frame(
          group1 = comp[1],
          group2 = comp[2],
          mean_diff = mean(group2_data, na.rm = TRUE) - mean(group1_data, na.rm = TRUE),
          p_value = test_result$p.value,
          adj_p_value = NA,
          stringsAsFactors = FALSE
        )
      })
      
      pairwise_df <- do.call(rbind, pairwise_results)
      pairwise_df$adj_p_value <- p.adjust(pairwise_df$p_value, method = p.adjust.method)
      pairwise_df$significant <- pairwise_df$adj_p_value < 0.05
      
      return(list(
        test = method,
        pairwise_results = pairwise_df
      ))
    }
  })
  
  names(results) <- score_cols
  return(results)
}

#' @title Calculate Cell Type Metabolic Diversity
#' 
#' @description Calculate metabolic diversity within cell types.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param score_cols Score columns
#' 
#' @return Data frame with diversity metrics
#' 
#' @export
#' @examples
#' \dontrun{
#' diversity <- calculateCellTypeDiversity(seurat_obj, cell_type_col = "cell_type")
#' }
calculateCellTypeDiversity <- function(seurat_obj,
                                      cell_type_col = "cell_type",
                                      score_cols = NULL) {
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  diversity <- seurat_obj@meta.data %>%
    dplyr::group_by(.data[[cell_type_col]]) %>%
    dplyr::summarize(
      metabolic_cv = mean(sapply(score_cols, function(col) {
        vals <- .data[[col]]
        if (all(is.na(vals)) || mean(vals, na.rm = TRUE) == 0) return(0)
        sd(vals, na.rm = TRUE) / mean(vals, na.rm = TRUE)
      }), na.rm = TRUE),
      phenotype_diversity = length(unique(metabolic_phenotype)),
      .groups = "drop"
    )
  
  colnames(diversity)[1] <- cell_type_col
  
  return(diversity)
}

#' @title Rank Cell Types by Metabolic Activity
#' 
#' @description Rank cell types by various metabolic metrics.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param rank_by Column to rank by
#' @param decreasing Sort order
#' 
#' @return Data frame with rankings
#' 
#' @export
#' @examples
#' \dontrun{
#' rankings <- rankCellTypesByMetabolism(seurat_obj, 
#'   cell_type_col = "cell_type",
#'   rank_by = "ATP_score")
#' }
rankCellTypesByMetabolism <- function(seurat_obj,
                                     cell_type_col = "cell_type",
                                     rank_by = "ATP_score",
                                     decreasing = TRUE) {
  
  rankings <- aggregateByCellType(seurat_obj, 
                                 cell_type_col = cell_type_col,
                                 score_cols = rank_by, 
                                 include_sd = FALSE,
                                 include_sem = FALSE)
  
  if (nrow(rankings) == 0) {
    stop("No valid rankings computed")
  }
  
  rankings <- rankings[order(rankings$mean, decreasing = decreasing), ]
  rankings$rank <- seq_len(nrow(rankings))
  
  return(rankings)
}

#' @title Get Dominant Phenotype per Cell Type
#' 
#' @description Get dominant metabolic phenotype for each cell type.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param phenotype_col Phenotype column
#' 
#' @return Data frame with dominant phenotypes
#' 
#' @export
#' @examples
#' \dontrun{
#' dominant <- getDominantPhenotypePerCellType(seurat_obj, cell_type_col = "cell_type")
#' print(dominant)
#' }
getDominantPhenotypePerCellType <- function(seurat_obj,
                                           cell_type_col = "cell_type",
                                           phenotype_col = "metabolic_phenotype") {
  
  if (!phenotype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Phenotype column not found"))
  }
  
  phenotype_dist <- getPhenotypeDistribution(seurat_obj, 
                                          phenotype_col = phenotype_col,
                                          group_col = cell_type_col)
  
  colnames(phenotype_dist)[1] <- cell_type_col
  
  dominant <- phenotype_dist %>%
    dplyr::group_by(.data[[cell_type_col]]) %>%
    dplyr::slice(which.max(pct)) %>%
    dplyr::rename(dominant_phenotype = 2, pct_dominant = pct) %>%
    dplyr::ungroup()
  
  return(as.data.frame(dominant))
}

#' @title Calculate Cell Type Metabolic Similarity
#' 
#' @description Calculate similarity between cell types based on metabolic profiles.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param score_cols Score columns
#' @param method Distance/similarity method: "correlation", "euclidean", "cosine"
#' 
#' @return Matrix of cell type similarities
#' 
#' @export
#' @examples
#' \dontrun{
#' similarity <- calculateCellTypeSimilarity(seurat_obj, cell_type_col = "cell_type")
#' }
calculateCellTypeSimilarity <- function(seurat_obj,
                                       cell_type_col = "cell_type",
                                       score_cols = NULL,
                                       method = c("correlation", "euclidean", "cosine")) {
  
  method <- match.arg(method)
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  if (length(score_cols) == 0) {
    stop("No valid score columns found")
  }
  
  aggregated <- aggregateByCellType(seurat_obj, cell_type_col = cell_type_col,
                                   score_cols = score_cols, include_sd = FALSE,
                                   include_sem = FALSE)
  
  if (nrow(aggregated) == 0) {
    stop("No aggregated data available")
  }
  
  cell_types <- unique(aggregated[[cell_type_col]])
  
  profile_matrix <- reshape2::dcast(aggregated, 
                                   as.formula(paste(cell_type_col, "~ metric")),
                                   value.var = "mean")
  
  rownames(profile_matrix) <- profile_matrix[[cell_type_col]]
  profile_matrix[[cell_type_col]] <- NULL
  
  profile_matrix <- as.matrix(profile_matrix)
  rownames(profile_matrix) <- cell_types
  
  if (method == "correlation") {
    similarity <- cor(t(profile_matrix), use = "pairwise.complete.obs")
  } else if (method == "euclidean") {
    dist_matrix <- as.matrix(dist(profile_matrix, method = "euclidean"))
    similarity <- 1 / (1 + dist_matrix)
  } else {
    similarity <- sapply(seq_len(nrow(profile_matrix)), function(i) {
      sapply(seq_len(nrow(profile_matrix)), function(j) {
        dot_product <- sum(profile_matrix[i, ] * profile_matrix[j, ], na.rm = TRUE)
        norm_i <- sqrt(sum(profile_matrix[i, ]^2, na.rm = TRUE))
        norm_j <- sqrt(sum(profile_matrix[j, ]^2, na.rm = TRUE))
        if (norm_i > 0 && norm_j > 0) {
          dot_product / (norm_i * norm_j)
        } else {
          0
        }
      })
    })
    rownames(similarity) <- rownames(profile_matrix)
    colnames(similarity) <- rownames(profile_matrix)
  }
  
  return(similarity)
}

#' @title Export Cell Type Metabolic Summary
#' 
#' @description Create comprehensive summary of cell type metabolism.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param output_dir Output directory
#' @param prefix File prefix
#' @param include_heatmap Save heatmap as PNG
#' 
#' @return List with summary components
#' 
#' @export
#' @examples
#' \dontrun{
#' summary <- exportCellTypeSummary(seurat_obj, 
#'   cell_type_col = "cell_type",
#'   output_dir = "results")
#' }
exportCellTypeSummary <- function(seurat_obj,
                                 cell_type_col = "cell_type",
                                 output_dir = ".",
                                 prefix = "celltype_metabolism",
                                 include_heatmap = FALSE) {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  aggregated <- aggregateByCellType(seurat_obj, cell_type_col = cell_type_col,
                                   include_sd = TRUE, include_sem = TRUE)
  
  aggregated_file <- file.path(output_dir, paste0(prefix, "_aggregated.csv"))
  write.csv(aggregated, aggregated_file, row.names = FALSE)
  
  phenotype_summary <- NULL
  if ("metabolic_phenotype" %in% colnames(seurat_obj@meta.data)) {
    phenotype_summary <- getDominantPhenotypePerCellType(seurat_obj, 
                                                        cell_type_col = cell_type_col)
    phenotype_file <- file.path(output_dir, paste0(prefix, "_phenotypes.csv"))
    write.csv(phenotype_summary, phenotype_file, row.names = FALSE)
  }
  
  similarity <- calculateCellTypeSimilarity(seurat_obj, cell_type_col = cell_type_col)
  similarity_file <- file.path(output_dir, paste0(prefix, "_similarity.csv"))
  write.csv(similarity, similarity_file)
  
  if (include_heatmap && requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    heatmap_file <- file.path(output_dir, paste0(prefix, "_heatmap.png"))
    png(heatmap_file, width = 800, height = 600)
    
    heatmap_matrix <- reshape2::dcast(aggregated, 
                                     as.formula(paste(cell_type_col, "~ metric")),
                                     value.var = "mean")
    rownames(heatmap_matrix) <- heatmap_matrix[[cell_type_col]]
    heatmap_matrix[[cell_type_col]] <- NULL
    heatmap_matrix <- as.matrix(heatmap_matrix)
    
    heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
    
    ComplexHeatmap::Heatmap(heatmap_matrix_scaled,
                           name = "Z-score",
                           cluster_rows = TRUE,
                           cluster_columns = TRUE) %>%
      ComplexHeatmap::draw()
    
    dev.off()
  }
  
  return(list(
    aggregated = aggregated,
    phenotypes = phenotype_summary,
    similarity = similarity,
    files_written = list(
      aggregated = aggregated_file,
      phenotypes = if (!is.null(phenotype_summary)) phenotype_file else NULL,
      similarity = similarity_file
    )
  ))
}
