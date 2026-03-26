#' @title Visualization Module
#' @description Create visualizations for metabolic analysis results
#' @name visualization
#' @rdname visualization
#' @keywords visualization plotting
NULL

#' @title Plot Metabolic Scores on UMAP
#' 
#' @description Overlay metabolic scores on UMAP visualization.
#' 
#' @param seurat_obj Seurat object with UMAP and scores
#' @param score_cols Score columns to plot. If NULL, uses default scores
#' @param reduction Reduction to use (umap, tsne, pca)
#' @param ncol Number of columns for multiple plots
#' @param combine Combine plots using patchwork
#' @param palette Color palette for scores
#' 
#' @return ggplot or list of ggplots
#' 
#' @export
#' @examples
#' \dontrun{
#' plotMetabolicUMAP(seurat_obj, 
#'   score_cols = c("glycolysis", "oxidative_phosphorylation", "ATP_score"))
#' }
plotMetabolicUMAP <- function(seurat_obj,
                             score_cols = NULL,
                             reduction = "umap",
                             ncol = 2,
                             combine = TRUE,
                             palette = "plasma") {
  
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction", reduction, "not found. Available:",
              paste(names(seurat_obj@reductions), collapse = ", ")))
  }
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score", "GOX_index")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  if (length(score_cols) == 0) {
    stop("No valid score columns found")
  }
  
  emb <- Seurat::Embeddings(seurat_obj, reduction = reduction)
  
  if (!is.matrix(emb) && !is.data.frame(emb)) {
    emb <- as.matrix(emb)
  }
  
  emb_df <- as.data.frame(emb)
  colnames(emb_df) <- c("dim1", "dim2")
  
  plot_list <- lapply(score_cols, function(col) {
    emb_df$score <- seurat_obj@meta.data[[col]]
    
    p <- ggplot2::ggplot(emb_df, ggplot2::aes_string(x = "dim1", y = "dim2", color = "score")) +
      ggplot2::geom_point(size = 0.5, alpha = 0.6) +
      ggplot2::scale_color_viridis_c(name = col, option = palette) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = col, 
                   x = paste0(toupper(reduction), "_1"), 
                   y = paste0(toupper(reduction), "_2")) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::guides(color = ggplot2::guide_colorbar(title.position = "top"))
    
    return(p)
  })
  
  if (combine && length(plot_list) > 1) {
    ncol_actual <- min(ncol, length(plot_list))
    nrow_actual <- ceiling(length(plot_list) / ncol_actual)
    
    combined <- patchwork::wrap_plots(plot_list, ncol = ncol_actual) +
      patchwork::plot_annotation(
        title = "Metabolic Score Distribution",
        theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold"))
      )
    return(combined)
  }
  
  return(plot_list)
}

#' @title Plot Metabolic Phenotype Distribution
#' 
#' @description Plot distribution of metabolic phenotypes.
#' 
#' @param seurat_obj Seurat object
#' @param phenotype_col Phenotype column
#' @param group_col Optional grouping column
#' @param colors Color palette for phenotypes
#' @param position Position for bars: "stack" or "fill"
#' 
#' @return ggplot
#' 
#' @export
#' @examples
#' \dontrun{
#' plotPhenotypeDistribution(seurat_obj, group_col = "cell_type")
#' }
plotPhenotypeDistribution <- function(seurat_obj,
                                     phenotype_col = "metabolic_phenotype",
                                     group_col = NULL,
                                     colors = NULL,
                                     position = c("stack", "fill")) {
  
  position <- match.arg(position)
  
  if (!phenotype_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Phenotype column not found"))
  }
  
  if (is.null(colors)) {
    colors <- c(
      "Glycolytic" = "#E64B35",
      "Oxidative" = "#4DBBD5",
      "Energetically Balanced" = "#00A087",
      "Energy-Stressed" = "#F39B7F",
      "Hypermetabolic" = "#8491B4"
    )
  }
  
  if (is.null(group_col)) {
    dist_data <- data.frame(
      phenotype = seurat_obj@meta.data[[phenotype_col]],
      stringsAsFactors = FALSE
    ) %>%
      dplyr::group_by(phenotype) %>%
      dplyr::summarize(Freq = dplyr::n(), .groups = "drop")
    
    dist_data$phenotype <- factor(dist_data$phenotype,
                                  levels = c("Glycolytic", "Oxidative", 
                                            "Energetically Balanced",
                                            "Energy-Stressed", "Hypermetabolic"))
    dist_data <- dist_data[!is.na(dist_data$phenotype), ]
    
    p <- ggplot2::ggplot(dist_data, ggplot2::aes(x = reorder(phenotype, -Freq), 
                                                  y = Freq, fill = phenotype)) +
      ggplot2::geom_bar(stat = "identity") +
      ggplot2::scale_fill_manual(values = colors, na.translate = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "Metabolic Phenotype", y = "Number of Cells",
                   title = "Metabolic Phenotype Distribution") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none",
        panel.grid.minor = ggplot2::element_blank()
      )
    
  } else {
    if (!group_col %in% colnames(seurat_obj@meta.data)) {
      stop(paste("Group column not found"))
    }
    
    dist_data <- getPhenotypeDistribution(seurat_obj, 
                                        phenotype_col = phenotype_col,
                                        group_col = group_col)
    colnames(dist_data)[1] <- "group"
    
    p <- ggplot2::ggplot(dist_data, ggplot2::aes_string(x = "group", y = "pct", 
                                                          fill = phenotype_col)) +
      ggplot2::geom_bar(stat = "identity", position = position) +
      ggplot2::scale_fill_manual(values = colors, na.translate = FALSE) +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = group_col, 
                   y = ifelse(position == "fill", "Percentage", "Cell Count"),
                   title = "Metabolic Phenotype Distribution by Cell Type") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "right",
        panel.grid.minor = ggplot2::element_blank()
      )
  }
  
  return(p)
}

#' @title Plot Glycolysis vs OXPHOS
#' 
#' @description Create scatter plot of glycolysis vs OXPHOS.
#' 
#' @param seurat_obj Seurat object
#' @param glycolysis_col Glycolysis score column
#' @param oxphos_col OXPHOS score column
#' @param color_by Color points by this variable
#' @param density Add density contours
#' @param sample_n Sample number of points for large datasets
#' 
#' @return ggplot
#' 
#' @export
#' @examples
#' \dontrun{
#' plotGlycolysisVsOxphos(seurat_obj, color_by = "cell_type")
#' }
plotGlycolysisVsOxphos <- function(seurat_obj,
                                   glycolysis_col = "glycolysis",
                                   oxphos_col = "oxidative_phosphorylation",
                                   color_by = NULL,
                                   density = TRUE,
                                   sample_n = NULL) {
  
  plot_data <- data.frame(
    glycolysis = seurat_obj@meta.data[[glycolysis_col]],
    oxphos = seurat_obj@meta.data[[oxphos_col]]
  )
  
  if (!is.null(sample_n) && nrow(plot_data) > sample_n) {
    set.seed(42)
    plot_data <- plot_data[sample(nrow(plot_data), sample_n), ]
  }
  
  if (!is.null(color_by) && color_by %in% colnames(seurat_obj@meta.data)) {
    if (!is.null(sample_n) && nrow(seurat_obj@meta.data) > sample_n) {
      set.seed(42)
      idx <- sample(nrow(seurat_obj@meta.data), sample_n)
      plot_data$group <- seurat_obj@meta.data[[color_by]][idx]
    } else {
      plot_data$group <- seurat_obj@meta.data[[color_by]]
    }
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = glycolysis, y = oxphos, color = group))
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = glycolysis, y = oxphos))
  }
  
  p <- p + ggplot2::geom_point(alpha = 0.3, size = 0.5)
  
  if (density) {
    p <- p + ggplot2::stat_density_2d(color = "gray30", alpha = 0.5, size = 0.5)
  }
  
  p <- p + 
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", 
                         color = "red", size = 0.8) +
    ggplot2::geom_smooth(method = "lm", color = "blue", se = FALSE, 
                        linetype = "dashed", size = 0.8) +
    ggplot2::scale_color_brewer(palette = "Set2") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Glycolysis Score", 
      y = "OXPHOS Score",
      title = "Glycolysis vs OXPHOS Activity",
      subtitle = "Points above line: Glycolytic, Below: Oxidative"
    ) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}

#' @title Plot ATP Distribution
#' 
#' @description Plot ATP score distribution.
#' 
#' @param seurat_obj Seurat object
#' @param atp_col ATP score column
#' @param group_col Grouping variable
#' @param type Plot type: "violin", "boxplot", "histogram", "density"
#' 
#' @return ggplot
#' 
#' @export
#' @examples
#' \dontrun{
#' plotATPDistribution(seurat_obj, group_col = "cell_type", type = "violin")
#' }
plotATPDistribution <- function(seurat_obj,
                               atp_col = "ATP_score",
                               group_col = NULL,
                               type = c("violin", "boxplot", "histogram", "density")) {
  
  type <- match.arg(type)
  
  if (!atp_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("ATP column not found"))
  }
  
  plot_data <- data.frame(
    atp = seurat_obj@meta.data[[atp_col]],
    cell_id = rownames(seurat_obj@meta.data)
  )
  
  if (!is.null(group_col) && group_col %in% colnames(seurat_obj@meta.data)) {
    plot_data$group <- seurat_obj@meta.data[[group_col]]
  }
  
  if (type == "histogram") {
    if (!is.null(group_col)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = atp, fill = group)) +
        ggplot2::geom_histogram(bins = 50, alpha = 0.7, position = "dodge") +
        ggplot2::labs(x = "ATP Score", y = "Count", 
                     title = "ATP Score Distribution")
    } else {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = atp)) +
        ggplot2::geom_histogram(bins = 50, fill = "#8491B4", alpha = 0.7) +
        ggplot2::labs(x = "ATP Score", y = "Count", 
                     title = "ATP Score Distribution")
    }
  } else if (type == "density") {
    if (!is.null(group_col)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = atp, fill = group)) +
        ggplot2::geom_density(alpha = 0.5) +
        ggplot2::labs(x = "ATP Score", y = "Density", 
                     title = "ATP Score Density")
    } else {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = atp)) +
        ggplot2::geom_density(fill = "#8491B4", alpha = 0.5) +
        ggplot2::labs(x = "ATP Score", y = "Density", 
                     title = "ATP Score Density")
    }
  } else {
    if (!is.null(group_col)) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "group", y = "atp", fill = "group"))
    } else {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(y = atp))
    }
    
    if (type == "violin") {
      p <- p + 
        ggplot2::geom_violin(fill = "#8491B4", alpha = 0.7) +
        ggplot2::geom_boxplot(width = 0.1, fill = "white", alpha = 0.8)
    } else {
      p <- p + ggplot2::geom_boxplot(fill = "#8491B4", alpha = 0.7)
    }
    
    p <- p + ggplot2::labs(y = "ATP Score", title = "ATP Score by Cell Type") +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
        legend.position = "none"
      )
  }
  
  p <- p + ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  
  return(p)
}

#' @title Plot GOX Index Distribution
#' 
#' @description Plot GOX (Glycolytic-OXPHOS) index distribution.
#' 
#' @param seurat_obj Seurat object
#' @param gox_col GOX index column
#' @param group_col Grouping variable
#' 
#' @return ggplot
#' 
#' @export
#' @examples
#' \dontrun{
#' plotGOXIndex(seurat_obj, group_col = "cell_type")
#' }
plotGOXIndex <- function(seurat_obj,
                        gox_col = "GOX_index",
                        group_col = NULL) {
  
  if (!gox_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("GOX index column not found"))
  }
  
  plot_data <- data.frame(
    gox = seurat_obj@meta.data[[gox_col]],
    cell_id = rownames(seurat_obj@meta.data)
  )
  
  if (!is.null(group_col) && group_col %in% colnames(seurat_obj@meta.data)) {
    plot_data$group <- seurat_obj@meta.data[[group_col]]
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = gox)) +
    ggplot2::geom_histogram(bins = 50, ggplot2::aes(fill = ..x..), alpha = 0.8) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "red", size = 1) +
    ggplot2::scale_fill_gradient2(
      low = "#4DBBD5", 
      mid = "#F0F0F0", 
      high = "#E64B35",
      name = "GOX Index"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "GOX Index", 
      y = "Count",
      title = "Glycolytic-OXPHOS (GOX) Index Distribution",
      subtitle = "Positive = Glycolytic, Negative = Oxidative"
    ) +
    ggplot2::theme(
      legend.position = "none",
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.null(group_col)) {
    p <- p + ggplot2::facet_wrap(~group, scales = "free_y")
  }
  
  return(p)
}

#' @title Plot Pathway Correlation
#' 
#' @description Plot correlation matrix of pathway scores.
#' 
#' @param seurat_obj Seurat object
#' @param score_cols Score columns
#' @param method Correlation method: "pearson" or "spearman"
#' @param label Include correlation values
#' 
#' @return ggplot
#' 
#' @export
#' @examples
#' \dontrun{
#' plotPathwayCorrelation(seurat_obj)
#' }
plotPathwayCorrelation <- function(seurat_obj,
                                  score_cols = NULL,
                                  method = c("pearson", "spearman"),
                                  label = TRUE) {
  
  method <- match.arg(method)
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score", 
                     "tca_cycle", "fatty_acid_oxidation", "GOX_index")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  if (length(score_cols) < 2) {
    stop("Need at least 2 score columns")
  }
  
  score_matrix <- as.matrix(seurat_obj@meta.data[, score_cols, drop = FALSE])
  cor_matrix <- cor(score_matrix, method = method, use = "pairwise")
  
  cor_df <- reshape2::melt(cor_matrix)
  colnames(cor_df) <- c("pathway1", "pathway2", "correlation")
  
  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = pathway1, y = pathway2, fill = correlation)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
      low = "#E64B35", 
      mid = "#FFFFFF", 
      high = "#4DBBD5",
      midpoint = 0
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Pathway Score Correlation (", toupper(method), ")", sep = ""),
      x = "", 
      y = ""
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 9),
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (label) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = round(correlation, 2)), 
      size = 3,
      color = ifelse(abs(cor_df$correlation) > 0.5, "white", "black")
    )
  }
  
  return(p)
}

#' @title Plot Flux Distribution
#' 
#' @description Plot flux score distribution.
#' 
#' @param seurat_obj Seurat object
#' @param flux_cols Flux columns to plot. If NULL, auto-detects
#' @param group_col Grouping variable
#' 
#' @return ggplot
#' 
#' @export
#' @examples
#' \dontrun{
#' plotFluxDistribution(seurat_obj, group_col = "cell_type")
#' }
plotFluxDistribution <- function(seurat_obj,
                                flux_cols = NULL,
                                group_col = NULL) {
  
  if (is.null(flux_cols)) {
    flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
    flux_cols <- flux_cols[!grepl("^(zscore|norm|rank)_", flux_cols)]
    flux_cols <- head(flux_cols, 4)
  }
  
  flux_cols <- intersect(flux_cols, colnames(seurat_obj@meta.data))
  
  if (length(flux_cols) == 0) {
    stop("No flux columns found")
  }
  
  plot_data <- seurat_obj@meta.data[, flux_cols, drop = FALSE]
  plot_data$cell_id <- rownames(seurat_obj@meta.data)
  
  if (!is.null(group_col) && group_col %in% colnames(seurat_obj@meta.data)) {
    plot_data$group <- seurat_obj@meta.data[[group_col]]
  }
  
  plot_long <- reshape2::melt(plot_data, id.vars = c("cell_id", group_col),
                              variable.name = "pathway", value.name = "flux")
  
  if (!is.null(group_col)) {
    p <- ggplot2::ggplot(plot_long, ggplot2::aes_string(x = "pathway", y = "flux", fill = "group")) +
      ggplot2::geom_boxplot(alpha = 0.7)
  } else {
    p <- ggplot2::ggplot(plot_long, ggplot2::aes_string(x = "pathway", y = "flux")) +
      ggplot2::geom_boxplot(fill = "#8491B4", alpha = 0.7)
  }
  
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Pathway", y = "Flux Score", title = "Metabolic Flux Distribution") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  if (!is.null(group_col)) {
    p <- p + ggplot2::theme(legend.position = "bottom")
  } else {
    p <- p + ggplot2::theme(legend.position = "none")
  }
  
  return(p)
}

#' @title Create Summary Dashboard
#' 
#' @description Create a comprehensive summary dashboard.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param reduction Reduction to use for UMAP
#' 
#' @return patchwork layout
#' 
#' @export
#' @examples
#' \dontrun{
#' dashboard <- createMetabolicDashboard(seurat_obj, cell_type_col = "cell_type")
#' print(dashboard)
#' }
createMetabolicDashboard <- function(seurat_obj,
                                    cell_type_col = "cell_type",
                                    reduction = "umap") {
  
  p1_list <- plotMetabolicUMAP(seurat_obj, 
                               score_cols = "ATP_score",
                               reduction = reduction,
                               combine = FALSE)
  p1 <- p1_list[[1]] + ggplot2::theme(legend.position = "none", 
                                      plot.title = ggplot2::element_blank())
  
  p2 <- plotPhenotypeDistribution(seurat_obj, group_col = cell_type_col) +
    ggplot2::theme(legend.position = "right", plot.title = ggplot2::element_blank())
  
  p3 <- plotGlycolysisVsOxphos(seurat_obj, density = TRUE) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_blank()
    )
  
  p4 <- plotATPDistribution(seurat_obj, group_col = cell_type_col, type = "violin") +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_blank()
    )
  
  p5 <- plotPathwayCorrelation(seurat_obj) + 
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y = ggplot2::element_text(size = 7)
    )
  
  layout <- "
  AAAB
  CCDE
  FFE#
  "
  
  dashboard <- p1 + p2 + p3 + p4 + p5 +
    patchwork::plot_layout(
      design = layout,
      widths = c(1, 1, 1, 1),
      heights = c(1, 0.8, 0.8)
    ) +
    patchwork::plot_annotation(
      title = "scMetaboFlux: Metabolic Analysis Dashboard",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 16, face = "bold")
      )
    )
  
  return(dashboard)
}

#' @title Plot Cell Type Metabolic Profiles
#' 
#' @description Plot metabolic profiles by cell type.
#' 
#' @param seurat_obj Seurat object
#' @param cell_type_col Cell type column
#' @param score_cols Score columns
#' @param top_n Number of top cell types to show
#' 
#' @return ggplot (bar chart)
#' 
#' @export
#' @examples
#' \dontrun{
#' plotCellTypeProfiles(seurat_obj, cell_type_col = "cell_type")
#' }
plotCellTypeProfiles <- function(seurat_obj,
                                cell_type_col = "cell_type",
                                score_cols = NULL,
                                top_n = 10) {
  
  if (is.null(score_cols)) {
    score_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
  }
  
  aggregated <- aggregateByCellType(seurat_obj, cell_type_col = cell_type_col,
                                   score_cols = score_cols, include_sd = FALSE,
                                   include_sem = FALSE)
  
  cell_type_counts <- table(seurat_obj@meta.data[[cell_type_col]])
  top_types <- names(sort(cell_type_counts, decreasing = TRUE))[1:min(top_n, length(cell_type_counts))]
  
  profile_matrix <- reshape2::dcast(aggregated, 
                                   as.formula(paste(cell_type_col, "~ metric")),
                                   value.var = "mean")
  
  rownames(profile_matrix) <- profile_matrix[[cell_type_col]]
  profile_matrix[[cell_type_col]] <- NULL
  profile_matrix <- as.matrix(profile_matrix)
  
  profile_matrix <- profile_matrix[rownames(profile_matrix) %in% top_types, , drop = FALSE]
  
  profile_long <- reshape2::melt(profile_matrix)
  colnames(profile_long) <- c("cell_type", "pathway", "value")
  
  profile_long <- profile_long %>%
    dplyr::group_by(cell_type) %>%
    dplyr::mutate(value_scaled = (value - min(value)) / (max(value) - min(value) + 1e-10))
  
  p <- ggplot2::ggplot(profile_long, ggplot2::aes(x = cell_type, y = pathway, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_viridis_c(name = "Score", option = "plasma") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = cell_type_col,
      y = "",
      title = "Metabolic Profile by Cell Type"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  return(p)
}
