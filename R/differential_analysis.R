#' @title Differential Metabolism Analysis Module
#' @description Compare metabolic states between conditions
#' @name differential_analysis
#' @rdname differential-analysis
#' @keywords differential analysis metabolism
NULL

#' @title Differential Metabolic Analysis
#' 
#' @description Compare metabolic scores between conditions (e.g., disease vs control).
#' 
#' @param seurat_obj Seurat object
#' @param condition_col Condition column (e.g., disease vs control)
#' @param control_group Control group name
#' @param case_group Case group name
#' @param score_cols Score columns to compare
#' @param method Statistical test method: "wilcox.test" or "t.test"
#' @param adjust_method P-value adjustment method (default "BH")
#' 
#' @return List with:
#'   \item{results}{Data frame with differential results}
#'   \item{condition_col}{Condition column used}
#'   \item{control_group}{Control group}
#'   \item{case_group}{Case group}
#' 
#' @export
#' @examples
#' \dontrun{
#' diff_results <- differentialMetabolicAnalysis(seurat_obj,
#'   condition_col = "disease_status",
#'   control_group = "Control",
#'   case_group = "Disease")
#' print(diff_results$results)
#' }
differentialMetabolicAnalysis <- function(seurat_obj,
                                         condition_col = "condition",
                                         control_group = NULL,
                                         case_group = NULL,
                                         score_cols = NULL,
                                         method = c("wilcox.test", "t.test"),
                                         adjust_method = "BH") {
  
  method <- match.arg(method)
  
  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop(paste("Condition column", condition_col, "not found in metadata"))
  }
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score",
                     "GOX_index", "metabolic_efficiency")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  if (length(score_cols) == 0) {
    stop("No valid score columns found")
  }
  
  metadata <- seurat_obj@meta.data
  
  if (!is.null(control_group) && !is.null(case_group)) {
    control_mask <- metadata[[condition_col]] == control_group
    case_mask <- metadata[[condition_col]] == case_group
    
    if (sum(control_mask) == 0) {
      stop(paste("No cells found in control group:", control_group))
    }
    if (sum(case_mask) == 0) {
      stop(paste("No cells found in case group:", case_group))
    }
  } else {
    unique_conditions <- unique(metadata[[condition_col]])
    if (length(unique_conditions) < 2) {
      stop("Need at least 2 conditions for comparison")
    }
    control_group <- unique_conditions[1]
    case_group <- unique_conditions[2]
    control_mask <- metadata[[condition_col]] == control_group
    case_mask <- metadata[[condition_col]] == case_group
  }
  
  results <- lapply(score_cols, function(col) {
    control_vals <- as.numeric(metadata[[col]][control_mask])
    case_vals <- as.numeric(metadata[[col]][case_mask])
    
    control_vals <- control_vals[!is.na(control_vals)]
    case_vals <- case_vals[!is.na(case_vals)]
    
    if (length(control_vals) < 2 || length(case_vals) < 2) {
      return(data.frame(
        metric = col,
        control_mean = NA,
        case_mean = NA,
        mean_diff = NA,
        control_sd = NA,
        case_sd = NA,
        effect_size = NA,
        p_value = NA,
        adj_p_value = NA,
        significance = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    if (method == "t.test") {
      test_result <- tryCatch({
        t.test(case_vals, control_vals)
      }, error = function(e) NULL)
    } else {
      test_result <- tryCatch({
        wilcox.test(case_vals, control_vals)
      }, error = function(e) NULL)
    }
    
    if (is.null(test_result)) {
      return(data.frame(
        metric = col,
        control_mean = mean(control_vals, na.rm = TRUE),
        case_mean = mean(case_vals, na.rm = TRUE),
        mean_diff = mean(case_vals, na.rm = TRUE) - mean(control_vals, na.rm = TRUE),
        control_sd = sd(control_vals, na.rm = TRUE),
        case_sd = sd(case_vals, na.rm = TRUE),
        effect_size = NA,
        p_value = NA,
        adj_p_value = NA,
        significance = NA,
        stringsAsFactors = FALSE
      ))
    }
    
    pooled_var <- (var(control_vals, na.rm = TRUE) + var(case_vals, na.rm = TRUE)) / 2
    effect_size <- if (pooled_var > 0) {
      (mean(case_vals, na.rm = TRUE) - mean(control_vals, na.rm = TRUE)) / sqrt(pooled_var)
    } else 0
    
    data.frame(
      metric = col,
      control_mean = mean(control_vals, na.rm = TRUE),
      case_mean = mean(case_vals, na.rm = TRUE),
      mean_diff = mean(case_vals, na.rm = TRUE) - mean(control_vals, na.rm = TRUE),
      control_sd = sd(control_vals, na.rm = TRUE),
      case_sd = sd(case_vals, na.rm = TRUE),
      effect_size = effect_size,
      p_value = test_result$p.value,
      adj_p_value = NA,
      significance = NA,
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  results_df$adj_p_value <- p.adjust(results_df$p_value, method = adjust_method)
  results_df$significance <- ifelse(results_df$adj_p_value < 0.05, "sig", "ns")
  
  results_df$control_mean <- round(results_df$control_mean, 4)
  results_df$case_mean <- round(results_df$case_mean, 4)
  results_df$mean_diff <- round(results_df$mean_diff, 4)
  results_df$effect_size <- round(results_df$effect_size, 4)
  
  return(list(
    results = results_df,
    condition_col = condition_col,
    control_group = control_group,
    case_group = case_group,
    n_control = sum(control_mask),
    n_case = sum(case_mask)
  ))
}

#' @title Test Metabolic Pathway Enrichment
#' 
#' @description Test for metabolic pathway enrichment in differential genes.
#' 
#' @param seurat_obj Seurat object
#' @param group_col Grouping column
#' @param group1 First group
#' @param group2 Second group
#' @param fc_threshold Fold change threshold (log scale)
#' @param padj_threshold Adjusted p-value threshold
#' @param assay Assay to use
#' 
#' @return Data frame with enrichment results
#' 
#' @export
#' @examples
#' \dontrun{
#' enrichment <- testMetabolicPathwayEnrichment(seurat_obj,
#'   group_col = "condition",
#'   group1 = "Control",
#'   group2 = "Disease")
#' }
testMetabolicPathwayEnrichment <- function(seurat_obj,
                                          group_col = "condition",
                                          group1 = "Control",
                                          group2 = "Case",
                                          fc_threshold = 0.5,
                                          padj_threshold = 0.05,
                                          assay = NULL) {
  
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop("Group column not found")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  metadata <- seurat_obj@meta.data
  expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay)
  
  group1_mask <- metadata[[group_col]] == group1
  group2_mask <- metadata[[group_col]] == group2
  
  if (sum(group1_mask) < 3 || sum(group2_mask) < 3) {
    stop("Need at least 3 cells per group")
  }
  
  group1_mean <- rowMeans(as.matrix(expr_data[, group1_mask]))
  group2_mean <- rowMeans(as.matrix(expr_data[, group2_mask]))
  
  fc <- group2_mean - group1_mean
  
  results <- lapply(names(metabolicGeneSets), function(pathway) {
    pathway_genes <- toupper(metabolicGeneSets[[pathway]])
    genes_upper <- toupper(rownames(expr_data))
    
    pathway_genes_in_data <- pathway_genes[pathway_genes %in% genes_upper]
    
    if (length(pathway_genes_in_data) == 0) return(NULL)
    
    pathway_fc <- fc[genes_upper %in% pathway_genes_in_data]
    
    n_up <- sum(pathway_fc > fc_threshold, na.rm = TRUE)
    n_down <- sum(pathway_fc < -fc_threshold, na.rm = TRUE)
    n_total <- length(pathway_genes_in_data)
    
    mean_fc <- mean(pathway_fc, na.rm = TRUE)
    
    data.frame(
      pathway = pathway,
      n_genes = n_total,
      n_up = n_up,
      n_down = n_down,
      pct_up = round(n_up / n_total * 100, 1),
      pct_down = round(n_down / n_total * 100, 1),
      mean_fc = round(mean_fc, 4),
      direction = ifelse(mean_fc > 0, "up", "down"),
      stringsAsFactors = FALSE
    )
  })
  
  results_df <- do.call(rbind, results)
  results_df <- results_df[!sapply(results_df, is.null), ]
  
  return(results_df)
}

#' @title Calculate Effect Sizes
#' 
#' @description Calculate Cohen's d and Glass's delta effect sizes.
#' 
#' @param seurat_obj Seurat object
#' @param group1_cells Logical vector for group 1
#' @param group2_cells Logical vector for group 2
#' @param score_cols Score columns
#' 
#' @return Data frame with effect sizes
#' 
#' @export
#' @examples
#' \dontrun{
#' effects <- calculateEffectSizes(seurat_obj,
#'   group1_cells = seurat_obj@meta.data$condition == "Control",
#'   group2_cells = seurat_obj@meta.data$condition == "Disease")
#' }
calculateEffectSizes <- function(seurat_obj,
                                 group1_cells,
                                 group2_cells,
                                 score_cols = NULL) {
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  results <- lapply(score_cols, function(col) {
    vals1 <- as.numeric(seurat_obj@meta.data[[col]][group1_cells])
    vals2 <- as.numeric(seurat_obj@meta.data[[col]][group2_cells])
    
    vals1 <- vals1[!is.na(vals1)]
    vals2 <- vals2[!is.na(vals2)]
    
    mean_diff <- mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE)
    
    pooled_sd <- sqrt((var(vals1, na.rm = TRUE) + var(vals2, na.rm = TRUE)) / 2)
    
    cohens_d <- if (pooled_sd > 0) mean_diff / pooled_sd else 0
    
    glass_delta <- if (sd(vals1, na.rm = TRUE) > 0) {
      (mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE)) / sd(vals1, na.rm = TRUE)
    } else 0
    
    data.frame(
      metric = col,
      mean_diff = round(mean_diff, 4),
      pooled_sd = round(pooled_sd, 4),
      cohens_d = round(cohens_d, 4),
      glass_delta = round(glass_delta, 4),
      stringsAsFactors = FALSE
    )
  })
  
  return(do.call(rbind, results))
}

#' @title Pairwise Metabolic Comparison
#' 
#' @description Perform pairwise comparisons across multiple groups.
#' 
#' @param seurat_obj Seurat object
#' @param group_col Grouping column
#' @param score_cols Score columns
#' @param method Statistical test: "wilcox.test" or "t.test"
#' @param p.adjust.method P-value adjustment method
#' 
#' @return Data frame with all pairwise comparisons
#' 
#' @export
#' @examples
#' \dontrun{
#' pairwise <- pairwiseMetabolicComparison(seurat_obj,
#'   group_col = "treatment",
#'   method = "t.test")
#' }
pairwiseMetabolicComparison <- function(seurat_obj,
                                       group_col = "condition",
                                       score_cols = NULL,
                                       method = c("wilcox.test", "t.test"),
                                       p.adjust.method = "BH") {
  
  method <- match.arg(method)
  
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop("Group column not found")
  }
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  groups <- unique(seurat_obj@meta.data[[group_col]])
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  results <- lapply(comparisons, function(comp) {
    group1 <- comp[1]
    group2 <- comp[2]
    
    group1_mask <- seurat_obj@meta.data[[group_col]] == group1
    group2_mask <- seurat_obj@meta.data[[group_col]] == group2
    
    comp_results <- lapply(score_cols, function(col) {
      vals1 <- as.numeric(seurat_obj@meta.data[[col]][group1_mask])
      vals2 <- as.numeric(seurat_obj@meta.data[[col]][group2_mask])
      
      vals1 <- vals1[!is.na(vals1)]
      vals2 <- vals2[!is.na(vals2)]
      
      if (length(vals1) < 2 || length(vals2) < 2) {
        return(data.frame(
          comparison = paste0(group2, "_vs_", group1),
          group1 = group1,
          group2 = group2,
          metric = col,
          mean_diff = NA,
          p_value = NA,
          adj_p_value = NA,
          significant = NA,
          stringsAsFactors = FALSE
        ))
      }
      
      if (method == "t.test") {
        test <- tryCatch({
          t.test(vals2, vals1)
        }, error = function(e) NULL)
      } else {
        test <- tryCatch({
          wilcox.test(vals2, vals1)
        }, error = function(e) NULL)
      }
      
      if (is.null(test)) {
        return(data.frame(
          comparison = paste0(group2, "_vs_", group1),
          group1 = group1,
          group2 = group2,
          metric = col,
          mean_diff = mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE),
          p_value = NA,
          adj_p_value = NA,
          significant = NA,
          stringsAsFactors = FALSE
        ))
      }
      
      data.frame(
        comparison = paste0(group2, "_vs_", group1),
        group1 = group1,
        group2 = group2,
        metric = col,
        mean_diff = round(mean(vals2, na.rm = TRUE) - mean(vals1, na.rm = TRUE), 4),
        p_value = test$p.value,
        adj_p_value = NA,
        significant = NA,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, comp_results)
  })
  
  results_df <- do.call(rbind, results)
  results_df$adj_p_value <- p.adjust(results_df$p_value, method = p.adjust.method)
  results_df$significant <- results_df$adj_p_value < 0.05
  
  return(results_df)
}

#' @title Create Metabolic Volcano Data
#' 
#' @description Prepare data for volcano plot visualization.
#' 
#' @param diff_results Differential analysis results
#' @param fc_threshold Fold change threshold for coloring
#' @param pval_threshold P-value threshold for significance
#' 
#' @return Data frame for volcano plot
#' 
#' @export
#' @examples
#' \dontrun{
#' volcano_data <- createMetabolicVolcanoData(diff_results)
#' }
createMetabolicVolcanoData <- function(diff_results,
                                       fc_threshold = 0.1,
                                       pval_threshold = 0.05) {
  
  if (!is.list(diff_results) || !"results" %in% names(diff_results)) {
    stop("diff_results must be output from differentialMetabolicAnalysis")
  }
  
  volcano_df <- diff_results$results
  
  volcano_df$log2FC <- log2(exp(volcano_df$mean_diff))
  volcano_df$neg_log10pval <- -log10(volcano_df$adj_p_value + 1e-100)
  
  volcano_df$category <- "ns"
  volcano_df$category[volcano_df$mean_diff > fc_threshold & 
                        volcano_df$adj_p_value < pval_threshold] <- "up"
  volcano_df$category[volcano_df$mean_diff < -fc_threshold & 
                        volcano_df$adj_p_value < pval_threshold] <- "down"
  
  return(volcano_df)
}

#' @title Test Metabolic Heterogeneity
#' 
#' @description Test for metabolic heterogeneity differences between groups.
#' 
#' @param seurat_obj Seurat object
#' @param group_col Grouping column
#' @param test_var Test for equality of variance: "fligner.test" or "leveneTest"
#' @param score_cols Score columns
#' 
#' @return Data frame with heterogeneity test results
#' 
#' @export
#' @examples
#' \dontrun{
#' heterogeneity <- testMetabolicHeterogeneity(seurat_obj, group_col = "condition")
#' }
testMetabolicHeterogeneity <- function(seurat_obj,
                                     group_col = "condition",
                                     test_var = c("fligner.test", "leveneTest"),
                                     score_cols = NULL) {
  
  test_var <- match.arg(test_var)
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  results <- lapply(score_cols, function(col) {
    formula <- as.formula(paste0("`", col, "`", " ~ ", group_col))
    
    if (test_var == "fligner.test") {
      test_result <- tryCatch({
        fligner.test(formula, data = seurat_obj@meta.data)
      }, error = function(e) NULL)
      
      if (is.null(test_result)) {
        return(data.frame(
          metric = col,
          test_statistic = NA,
          p_value = NA,
          var_ratio = NA,
          stringsAsFactors = FALSE
        ))
      }
      
      test_stat <- test_result$statistic
      p_val <- test_result$p.value
    } else {
      if (requireNamespace("car", quietly = TRUE)) {
        test_result <- tryCatch({
          car::leveneTest(formula, data = seurat_obj@meta.data)
        }, error = function(e) NULL)
        
        if (is.null(test_result)) {
          return(data.frame(
            metric = col,
            test_statistic = NA,
            p_value = NA,
            var_ratio = NA,
            stringsAsFactors = FALSE
          ))
        }
        
        test_stat <- test_result$`F value`[1]
        p_val <- test_result$`Pr(>F)`[1]
      } else {
        test_result <- tryCatch({
          fligner.test(formula, data = seurat_obj@meta.data)
        }, error = function(e) NULL)
        
        if (is.null(test_result)) {
          return(data.frame(
            metric = col,
            test_statistic = NA,
            p_value = NA,
            var_ratio = NA,
            stringsAsFactors = FALSE
          ))
        }
        
        test_stat <- test_result$statistic
        p_val <- test_result$p.value
      }
    }
    
    groups <- unique(seurat_obj@meta.data[[group_col]])
    variances <- sapply(groups, function(g) {
      var(seurat_obj@meta.data[[col]][seurat_obj@meta.data[[group_col]] == g], 
          na.rm = TRUE)
    })
    
    var_ratio <- if (min(variances, na.rm = TRUE) > 0) {
      max(variances, na.rm = TRUE) / min(variances, na.rm = TRUE)
    } else NA
    
    data.frame(
      metric = col,
      test_statistic = round(test_stat, 4),
      p_value = p_val,
      var_ratio = round(var_ratio, 4),
      stringsAsFactors = FALSE
    )
  })
  
  return(do.call(rbind, results))
}

#' @title Calculate Metabolic Dysregulation Score
#' 
#' @description Calculate composite score of metabolic dysregulation relative
#' to a control group.
#' 
#' @param seurat_obj Seurat object
#' @param control_group Control group name
#' @param condition_col Condition column
#' @param score_cols Score columns to include
#' @param name Column name for dysregulation score
#' 
#' @return Seurat object with dysregulation scores
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- calculateMetabolicDysregulationScore(seurat_obj,
#'   control_group = "Control",
#'   condition_col = "condition")
#' }
calculateMetabolicDysregulationScore <- function(seurat_obj,
                                                 control_group = "Control",
                                                 condition_col = "condition",
                                                 score_cols = NULL,
                                                 name = "dysregulation_score") {
  
  if (is.null(score_cols)) {
    default_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
    score_cols <- intersect(default_cols, colnames(seurat_obj@meta.data))
  }
  
  score_cols <- intersect(score_cols, colnames(seurat_obj@meta.data))
  
  if (length(score_cols) == 0) {
    stop("No valid score columns found")
  }
  
  control_mask <- seurat_obj@meta.data[[condition_col]] == control_group
  
  if (sum(control_mask) < 3) {
    stop("Need at least 3 control cells for reference")
  }
  
  reference_values <- sapply(score_cols, function(col) {
    mean(seurat_obj@meta.data[[col]][control_mask], na.rm = TRUE)
  })
  
  dysregulation <- sapply(seq_len(nrow(seurat_obj@meta.data)), function(i) {
    mean(sapply(score_cols, function(col) {
      val <- seurat_obj@meta.data[[col]][i]
      ref <- reference_values[col]
      if (ref != 0) {
        abs(val - ref) / abs(ref)
      } else {
        abs(val)
      }
    }), na.rm = TRUE)
  })
  
  seurat_obj@meta.data[[name]] <- as.numeric(dysregulation)
  
  return(seurat_obj)
}
