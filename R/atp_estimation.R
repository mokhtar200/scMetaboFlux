#' @title ATP Production Approximation Module
#' @description Core novelty: Convert transcriptomics to quantitative ATP estimates
#' @name atp_estimation
#' @rdname atp-estimation
#' @keywords ATP energy metabolism
NULL

#' @title Estimate ATP Production Rate
#' 
#' @description Estimate cellular ATP production from metabolic pathway scores.
#' This is the core innovation of scMetaboFlux - converting transcriptomics
#' data into quantitative ATP production estimates.
#' 
#' The ATP model is based on the formula:
#' \code{ATP_total = w1 * Glycolysis + w2 * OXPHOS + w3 * FAO}
#' 
#' Where weights reflect relative ATP yield per pathway:
#' \itemize{
#'   \item Glycolysis: ~2 ATP/glucose (fast, inefficient)
#'   \item OXPHOS: ~30-32 ATP/glucose (slow, efficient)
#'   \item FAO: ~14 ATP per 2-carbon unit
#' }
#' 
#' @param seurat_obj Seurat object with pathway scores
#' @param glycolysis_col Column name for glycolysis score (default "glycolysis")
#' @param oxphos_col Column name for OXPHOS score (default "oxidative_phosphorylation")
#' @param fao_col Column name for fatty acid oxidation score (default "fatty_acid_oxidation")
#' @param weights Named numeric vector of pathway weights. If NULL, uses 
#'   atpYieldCoefficients normalized.
#' @param normalize Normalize ATP score to 0-1 range (default TRUE)
#' @param method ATP estimation method: 
#'   \itemize{
#'     \item "weighted_sum": Standard weighted sum of pathway scores
#'     \item "flux_based": Use theoretical ATP yields directly
#'     \item "relative": Relative to maximum observed values
#'   }
#' @param name Column name for ATP score (default "ATP_score")
#' @param verbose Print progress messages
#' 
#' @return Seurat object with ATP scores added to meta.data:
#'   \item{ATP_score}{Estimated ATP production score}
#'   \item{ATP_score_glycolytic_contribution}{Percentage from glycolysis}
#'   \item{ATP_score_oxidative_contribution}{Percentage from OXPHOS}
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- estimateATPProduction(seurat_obj)
#' seurat_obj <- estimateATPProduction(seurat_obj, 
#'   weights = c(glycolysis = 0.3, oxphos = 0.7),
#'   normalize = TRUE,
#'   method = "flux_based")
#' }
estimateATPProduction <- function(seurat_obj,
                                 glycolysis_col = "glycolysis",
                                 oxphos_col = "oxidative_phosphorylation",
                                 fao_col = "fatty_acid_oxidation",
                                 weights = NULL,
                                 normalize = TRUE,
                                 method = c("weighted_sum", "flux_based", "relative"),
                                 name = "ATP_score",
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
    message("[scMetaboFlux] Estimating ATP production...")
  }
  
  if (is.null(weights)) {
    glycolysis_yield <- atpYieldCoefficients["glycolysis"]
    oxphos_yield <- atpYieldCoefficients["oxidative_phosphorylation"]
    fao_yield <- atpYieldCoefficients["fatty_acid_oxidation"]
    
    total_yield <- glycolysis_yield + oxphos_yield + fao_yield
    
    weights <- c(
      glycolysis = glycolysis_yield / total_yield,
      oxphos = oxphos_yield / total_yield,
      fao = fao_yield / total_yield
    )
  } else {
    if (!all(c("glycolysis", "oxphos") %in% names(weights))) {
      stop("weights must contain 'glycolysis' and 'oxphos' named elements")
    }
  }
  
  glycolysis_scores <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos_scores <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  if (!is.null(fao_col) && fao_col %in% colnames(seurat_obj@meta.data)) {
    fao_scores <- as.numeric(seurat_obj@meta.data[[fao_col]])
  } else {
    fao_scores <- rep(0, nrow(seurat_obj@meta.data))
  }
  
  glycolysis_scores[is.na(glycolysis_scores)] <- 0
  oxphos_scores[is.na(oxphos_scores)] <- 0
  fao_scores[is.na(fao_scores)] <- 0
  
  atp_scores <- switch(method,
    "weighted_sum" = {
      (weights["glycolysis"] * glycolysis_scores +
        weights["oxphos"] * oxphos_scores +
        weights["fao"] * fao_scores)
    },
    "flux_based" = {
      glycolysis_flux <- glycolysis_scores * atpYieldCoefficients["glycolysis"]
      oxphos_flux <- oxphos_scores * atpYieldCoefficients["oxidative_phosphorylation"]
      fao_flux <- fao_scores * atpYieldCoefficients["fatty_acid_oxidation"]
      glycolysis_flux + oxphos_flux + fao_flux
    },
    "relative" = {
      max_glycolysis <- max(glycolysis_scores, na.rm = TRUE)
      max_oxphos <- max(oxphos_scores, na.rm = TRUE)
      
      if (is.finite(max_glycolysis) && max_glycolysis > 0) {
        relative_glycolysis <- glycolysis_scores / max_glycolysis
      } else {
        relative_glycolysis <- glycolysis_scores
      }
      
      if (is.finite(max_oxphos) && max_oxphos > 0) {
        relative_oxphos <- oxphos_scores / max_oxphos
      } else {
        relative_oxphos <- oxphos_scores
      }
      
      (weights["glycolysis"] * relative_glycolysis +
        weights["oxphos"] * relative_oxphos)
    }
  )
  
  atp_scores <- as.numeric(atp_scores)
  
  if (normalize) {
    valid_scores <- atp_scores[is.finite(atp_scores)]
    if (length(valid_scores) > 0) {
      min_atp <- min(valid_scores)
      max_atp <- max(valid_scores)
      
      if (max_atp > min_atp) {
        atp_scores <- (atp_scores - min_atp) / (max_atp - min_atp)
      } else {
        atp_scores <- rep(0.5, length(atp_scores))
      }
    } else {
      atp_scores <- rep(0.5, length(atp_scores))
    }
  }
  
  seurat_obj@meta.data[[name]] <- atp_scores
  
  glycolytic_contrib <- (weights["glycolysis"] * glycolysis_scores) / 
    (atp_scores + 1e-10) * 100
  oxidative_contrib <- (weights["oxphos"] * oxphos_scores) / 
    (atp_scores + 1e-10) * 100
  
  seurat_obj@meta.data[[paste0(name, "_glycolytic_contribution")]] <- 
    pmax(0, pmin(100, glycolytic_contrib))
  
  seurat_obj@meta.data[[paste0(name, "_oxidative_contribution")]] <- 
    pmax(0, pmin(100, oxidative_contrib))
  
  if (verbose) {
    message("[scMetaboFlux] ATP estimation complete!")
    message("[scMetaboFlux] Mean ATP score: ", round(mean(atp_scores, na.rm = TRUE), 4))
    message("[scMetaboFlux] Median ATP score: ", round(median(atp_scores, na.rm = TRUE), 4))
    message("[scMetaboFlux] SD ATP score: ", round(sd(atp_scores, na.rm = TRUE), 4))
  }
  
  return(seurat_obj)
}

#' @title Calculate ATP Yield Ratio
#' 
#' @description Calculate ratio between glycolytic and oxidative ATP production.
#' Values > 1 indicate glycolytic dominance, values < 1 indicate oxidative dominance.
#' 
#' @param seurat_obj Seurat object with pathway scores
#' @param glycolysis_col Glycolysis score column
#' @param oxphos_col OXPHOS score column
#' @param name Column name for ratio
#' @param add_category Add categorical label (default TRUE)
#' 
#' @return Seurat object with ratio added to meta.data
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- calculateATPYieldRatio(seurat_obj)
#' }
calculateATPYieldRatio <- function(seurat_obj,
                                   glycolysis_col = "glycolysis",
                                   oxphos_col = "oxidative_phosphorylation",
                                   name = "ATP_yield_ratio",
                                   add_category = TRUE) {
  
  glycolysis_scores <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos_scores <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  ratio <- glycolysis_scores / (oxphos_scores + 1e-10)
  ratio[is.infinite(ratio) | is.nan(ratio)] <- NA
  
  seurat_obj@meta.data[[name]] <- ratio
  
  if (add_category) {
    seurat_obj@meta.data[[paste0(name, "_category")]] <- cut(
      ratio,
      breaks = c(-Inf, 0.5, 1.5, Inf),
      labels = c("Oxidative", "Balanced", "Glycolytic")
    )
  }
  
  return(seurat_obj)
}

#' @title Estimate Metabolic Efficiency
#' 
#' @description Calculate metabolic efficiency index as ATP per unit of 
#' total pathway activity.
#' 
#' @param seurat_obj Seurat object
#' @param atp_col ATP score column
#' @param glycolysis_col Glycolysis score column
#' @param oxphos_col OXPHOS score column
#' @param name Column name for efficiency
#' 
#' @return Seurat object with efficiency added
#' 
#' @export
estimateMetabolicEfficiency <- function(seurat_obj,
                                        atp_col = "ATP_score",
                                        glycolysis_col = "glycolysis",
                                        oxphos_col = "oxidative_phosphorylation",
                                        name = "metabolic_efficiency") {
  
  atp_scores <- as.numeric(seurat_obj@meta.data[[atp_col]])
  
  glycolysis <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  total_activity <- glycolysis + oxphos
  
  efficiency <- atp_scores / (total_activity + 1e-10)
  efficiency[is.infinite(efficiency) | is.nan(efficiency)] <- 0
  
  seurat_obj@meta.data[[name]] <- efficiency
  
  return(seurat_obj)
}

#' @title Classify ATP Levels
#' 
#' @description Classify cells into ATP level categories (Low, Medium, High)
#' based on quantile thresholds.
#' 
#' @param seurat_obj Seurat object with ATP scores
#' @param atp_col ATP score column
#' @param quantiles Numeric vector of two quantiles for classification
#' @param name Column name for classification
#' 
#' @return Seurat object with classification added
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- classifyATPLevels(seurat_obj)
#' table(seurat_obj@meta.data$ATP_level)
#' }
classifyATPLevels <- function(seurat_obj,
                             atp_col = "ATP_score",
                             quantiles = c(0.25, 0.75),
                             name = "ATP_level") {
  
  atp_scores <- as.numeric(seurat_obj@meta.data[[atp_col]])
  atp_scores[is.na(atp_scores)] <- median(atp_scores, na.rm = TRUE)
  
  q <- quantile(atp_scores, probs = quantiles, na.rm = TRUE)
  
  seurat_obj@meta.data[[name]] <- cut(
    atp_scores,
    breaks = c(-Inf, q[1], q[2], Inf),
    labels = c("Low", "Medium", "High")
  )
  
  return(seurat_obj)
}

#' @title Calculate Glycolytic-OXPHOS Index (GOX Index)
#' 
#' @description Calculate GOX index: positive values indicate glycolytic 
#' phenotype, negative values indicate oxidative phenotype.
#' 
#' @param seurat_obj Seurat object
#' @param glycolysis_col Glycolysis score column
#' @param oxphos_col OXPHOS score column
#' @param z_score Z-score normalize the index
#' @param name Column name
#' @param verbose Print progress messages
#' 
#' @return Seurat object with GOX index added
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- calculateGOXIndex(seurat_obj)
#' }
calculateGOXIndex <- function(seurat_obj,
                              glycolysis_col = "glycolysis",
                              oxphos_col = "oxidative_phosphorylation",
                              z_score = FALSE,
                              name = "GOX_index",
                              verbose = TRUE) {
  
  glycolysis_scores <- as.numeric(seurat_obj@meta.data[[glycolysis_col]])
  oxphos_scores <- as.numeric(seurat_obj@meta.data[[oxphos_col]])
  
  glycolysis_scores[is.na(glycolysis_scores)] <- 0
  oxphos_scores[is.na(oxphos_scores)] <- 0
  
  if (z_score) {
    gly_mean <- mean(glycolysis_scores, na.rm = TRUE)
    gly_sd <- sd(glycolysis_scores, na.rm = TRUE)
    gly_sd <- ifelse(gly_sd > 0, gly_sd, 1)
    
    oxphos_mean <- mean(oxphos_scores, na.rm = TRUE)
    oxphos_sd <- sd(oxphos_scores, na.rm = TRUE)
    oxphos_sd <- ifelse(oxphos_sd > 0, oxphos_sd, 1)
    
    glycolysis_z <- (glycolysis_scores - gly_mean) / gly_sd
    oxphos_z <- (oxphos_scores - oxphos_mean) / oxphos_sd
    
    gox_index <- glycolysis_z - oxphos_z
  } else {
    gox_range <- max(glycolysis_scores, na.rm = TRUE) - min(glycolysis_scores, na.rm = TRUE)
    oxphos_range <- max(oxphos_scores, na.rm = TRUE) - min(oxphos_scores, na.rm = TRUE)
    
    if (gox_range > 0) {
      glycolysis_norm <- (glycolysis_scores - min(glycolysis_scores, na.rm = TRUE)) / gox_range
    } else {
      glycolysis_norm <- glycolysis_scores
    }
    
    if (oxphos_range > 0) {
      oxphos_norm <- (oxphos_scores - min(oxphos_scores, na.rm = TRUE)) / oxphos_range
    } else {
      oxphos_norm <- oxphos_scores
    }
    
    gox_index <- glycolysis_norm - oxphos_norm
  }
  
  gox_index[is.na(gox_index)] <- 0
  
  seurat_obj@meta.data[[name]] <- as.numeric(gox_index)
  
  seurat_obj@meta.data[[paste0(name, "_category")]] <- cut(
    gox_index,
    breaks = c(-Inf, -0.5, 0.5, Inf),
    labels = c("Oxidative", "Balanced", "Glycolytic")
  )
  
  return(seurat_obj)
}

#' @title Estimate Pseudo-ATP Production Rate
#' 
#' @description Model ATP production using Michaelis-Menten-like kinetics
#' for more biologically realistic flux estimation.
#' 
#' @param seurat_obj Seurat object
#' @param pathway_scores Named list of pathway score columns
#' @param km Named vector of Michaelis constants for each pathway
#' @param vmax Named vector of maximum velocities for each pathway
#' @param name Column name
#' 
#' @return Seurat object with pseudo-ATP rate
#' 
#' @export
estimatePseudoATPRate <- function(seurat_obj,
                                  pathway_scores = NULL,
                                  km = c(glycolysis = 0.5, oxphos = 0.7),
                                  vmax = c(glycolysis = 2, oxphos = 32),
                                  name = "ATP_pseudo_rate") {
  
  if (is.null(pathway_scores)) {
    if ("glycolysis" %in% colnames(seurat_obj@meta.data) &&
        "oxidative_phosphorylation" %in% colnames(seurat_obj@meta.data)) {
      pathway_scores <- list(
        glycolysis = seurat_obj@meta.data$glycolysis,
        oxphos = seurat_obj@meta.data$oxidative_phosphorylation
      )
    } else {
      stop("Pathway scores not found and pathway_scores not provided")
    }
  }
  
  atp_rate <- sapply(names(pathway_scores), function(pathway) {
    s <- as.numeric(pathway_scores[[pathway]])
    s[is.na(s)] <- 0
    
    km_val <- km[pathway]
    vmax_val <- vmax[pathway]
    
    if (is.na(km_val)) km_val <- 0.5
    if (is.na(vmax_val)) vmax_val <- 2
    
    (vmax_val * s) / (km_val + s)
  })
  
  if (is.matrix(atp_rate)) {
    atp_rate <- rowSums(atp_rate, na.rm = TRUE)
  }
  
  atp_rate <- as.numeric(atp_rate)
  min_rate <- min(atp_rate, na.rm = TRUE)
  max_rate <- max(atp_rate, na.rm = TRUE)
  
  if (max_rate > min_rate) {
    atp_rate <- (atp_rate - min_rate) / (max_rate - min_rate)
  } else {
    atp_rate <- rep(0.5, length(atp_rate))
  }
  
  seurat_obj@meta.data[[name]] <- atp_rate
  
  return(seurat_obj)
}

#' @title Calculate Metabolic Power Output
#' 
#' @description Calculate metabolic power as ATP per unit time, accounting
#' for the number of expressed metabolic genes.
#' 
#' @param seurat_obj Seurat object
#' @param atp_col ATP score column
#' @param expression_rate Estimated expression rate
#' @param assay Assay to use for gene counts
#' @param name Column name
#' 
#' @return Seurat object with power output
#' 
#' @export
calculateMetabolicPowerOutput <- function(seurat_obj,
                                          atp_col = "ATP_score",
                                          expression_rate = 1,
                                          assay = NULL,
                                          name = "metabolic_power") {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  atp_scores <- as.numeric(seurat_obj@meta.data[[atp_col]])
  
  expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay)
  
  n_expressed_genes <- colSums(Matrix::rowSums(expr_data) > 0)
  
  power <- atp_scores * n_expressed_genes * expression_rate
  
  power <- as.numeric(power)
  min_power <- min(power, na.rm = TRUE)
  max_power <- max(power, na.rm = TRUE)
  
  if (max_power > min_power) {
    power <- (power - min_power) / (max_power - min_power)
  } else {
    power <- rep(0.5, length(power))
  }
  
  seurat_obj@meta.data[[name]] <- power
  
  return(seurat_obj)
}

#' @title Bootstrap ATP Confidence Intervals
#' 
#' @description Calculate bootstrap confidence intervals for ATP estimates
#' to assess uncertainty.
#' 
#' @param seurat_obj Seurat object
#' @param atp_col ATP score column
#' @param n_boot Number of bootstrap iterations (default 100)
#' @param ci Confidence interval level (default 0.95)
#' @param seed Random seed for reproducibility
#' 
#' @return Data frame with CI estimates
#' 
#' @export
#' @examples
#' \dontrun{
#' ci <- bootstrapATPConfidenceIntervals(seurat_obj, n_boot = 1000)
#' print(ci)
#' }
bootstrapATPConfidenceIntervals <- function(seurat_obj,
                                            atp_col = "ATP_score",
                                            n_boot = 100,
                                            ci = 0.95,
                                            seed = 42) {
  
  set.seed(seed)
  
  atp_scores <- as.numeric(seurat_obj@meta.data[[atp_col]])
  atp_scores <- atp_scores[!is.na(atp_scores)]
  
  if (length(atp_scores) == 0) {
    stop("No valid ATP scores found")
  }
  
  n_cells <- length(atp_scores)
  
  boot_means <- replicate(n_boot, {
    sample_idx <- sample(1:n_cells, n_cells, replace = TRUE)
    mean(atp_scores[sample_idx], na.rm = TRUE)
  })
  
  alpha <- (1 - ci) / 2
  ci_lower <- quantile(boot_means, alpha)
  ci_upper <- quantile(boot_means, 1 - alpha)
  
  results <- data.frame(
    metric = c("mean", "median", "CI_lower", "CI_upper", "se", "n_cells"),
    value = c(
      mean(atp_scores, na.rm = TRUE),
      median(atp_scores, na.rm = TRUE),
      as.numeric(ci_lower),
      as.numeric(ci_upper),
      sd(boot_means),
      n_cells
    ),
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' @title Get ATP Statistics by Group
#' 
#' @description Get ATP statistics stratified by a metadata column.
#' 
#' @param seurat_obj Seurat object
#' @param atp_col ATP score column
#' @param group_col Grouping column
#' 
#' @return Data frame with statistics per group
#' 
#' @export
#' @examples
#' \dontrun{
#' stats <- getATPStatsByGroup(seurat_obj, group_col = "cell_type")
#' print(stats)
#' }
getATPStatsByGroup <- function(seurat_obj,
                               atp_col = "ATP_score",
                               group_col = "seurat_clusters") {
  
  if (!group_col %in% colnames(seurat_obj@meta.data)) {
    stop("Group column not found in metadata")
  }
  
  if (!atp_col %in% colnames(seurat_obj@meta.data)) {
    stop("ATP column not found in metadata")
  }
  
  metadata <- seurat_obj@meta.data
  
  stats <- metadata %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarize(
      mean_atp = mean(.data[[atp_col]], na.rm = TRUE),
      median_atp = median(.data[[atp_col]], na.rm = TRUE),
      sd_atp = sd(.data[[atp_col]], na.rm = TRUE),
      min_atp = min(.data[[atp_col]], na.rm = TRUE),
      max_atp = max(.data[[atp_col]], na.rm = TRUE),
      n_cells = dplyr::n(),
      .groups = "drop"
    )
  
  stats$mean_atp <- round(stats$mean_atp, 4)
  stats$median_atp <- round(stats$median_atp, 4)
  stats$sd_atp <- round(stats$sd_atp, 4)
  
  return(stats)
}
