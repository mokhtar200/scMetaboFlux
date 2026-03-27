#' @title Main Wrapper Functions and Utilities
#' @description Primary workflow functions for scMetaboFlux
#' @name main_workflow
#' @rdname main-workflow
#' @keywords workflow analysis pipeline
#' @import methods stats grDevices utils
NULL

#' @title Run Complete Metabolic Analysis
#' 
#' @description Run the complete scMetaboFlux workflow on a Seurat object.
#' This is the main wrapper function that executes the full analysis pipeline:
#' 1. Pathway scoring
#' 2. ATP estimation
#' 3. Flux approximation
#' 4. Phenotype classification
#' 
#' @param seurat_obj Seurat object with expression data
#' @param organism Organism: "human" or "mouse"
#' @param pathways_to_analyze Pathways to include in analysis
#' @param scoring_method Method for pathway scoring: "AUCell", "GSVA", "mean", "weighted_mean"
#' @param classify_phenotype Include phenotype classification
#' @param cell_type_col Cell type annotation column
#' @param condition_col Condition column for differential analysis
#' @param compute_flux Compute flux approximation
#' @param parallel_cores Number of cores for parallel processing
#' @param verbose Print progress messages
#' 
#' @return Seurat object with all metabolic scores and classifications
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- runMetabolicAnalysis(seurat_obj,
#'   organism = "human",
#'   pathways_to_analyze = c("glycolysis", "tca_cycle", "oxidative_phosphorylation"),
#'   scoring_method = "AUCell",
#'   classify_phenotype = TRUE,
#'   cell_type_col = "cell_type",
#'   condition_col = "disease_status")
#' }
runMetabolicAnalysis <- function(seurat_obj,
                                organism = c("human", "mouse"),
                                pathways_to_analyze = NULL,
                                scoring_method = c("AUCell", "GSVA", "mean", "weighted_mean"),
                                classify_phenotype = TRUE,
                                cell_type_col = NULL,
                                condition_col = NULL,
                                compute_flux = TRUE,
                                parallel_cores = NULL,
                                verbose = TRUE) {
  
  organism <- match.arg(organism)
  scoring_method <- match.arg(scoring_method)
  
  if (is.null(parallel_cores)) {
    parallel_cores <- max(1, parallel::detectCores() - 1)
  }
  parallel_cores <- min(parallel_cores, parallel::detectCores())
  
  if (verbose && parallel_cores > 1) {
    message("[scMetaboFlux] Using ", parallel_cores, " cores for parallel processing")
  }
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (verbose) {
    message("========================================")
    message("[scMetaboFlux] Starting Metabolic Analysis")
    message("[scMetaboFlux] Version: ", as.character(packageVersion("scMetaboFlux")))
    message("========================================")
  }
  
  if (is.null(pathways_to_analyze)) {
    pathways_to_analyze <- c("glycolysis", "tca_cycle", "oxidative_phosphorylation",
                             "fatty_acid_oxidation")
  }
  
  gene_sets <- metabolicGeneSets[pathways_to_analyze]
  gene_sets <- gene_sets[!sapply(gene_sets, is.null)]
  
  if (verbose) {
    message("\n[1/5] Computing pathway activity scores...")
    message("      Pathways: ", paste(names(gene_sets), collapse = ", "))
    message("      Method: ", scoring_method)
  }
  
  seurat_obj <- computePathwayScores(seurat_obj,
                                    gene_sets = gene_sets,
                                    method = scoring_method,
                                    parallel_cores = parallel_cores,
                                    verbose = verbose)
  
  if (verbose) {
    message("\n[2/5] Estimating ATP production...")
  }
  
  seurat_obj <- estimateATPProduction(seurat_obj,
                                     normalize = TRUE,
                                     verbose = verbose)
  
  seurat_obj <- calculateATPYieldRatio(seurat_obj)
  seurat_obj <- calculateGOXIndex(seurat_obj)
  
  if (compute_flux) {
    if (verbose) {
      message("\n[3/5] Computing metabolic flux approximation...")
    }
    
    seurat_obj <- computeMetabolicFlux(seurat_obj,
                                      gene_sets = gene_sets,
                                      use_rate_limiting = TRUE,
                                      normalize = TRUE,
                                      verbose = verbose)
    seurat_obj <- calculateFluxContribution(seurat_obj)
  }
  
  if (classify_phenotype) {
    if (verbose) {
      message("\n[4/5] Classifying metabolic phenotypes...")
    }
    
    seurat_obj <- classifyMetabolicPhenotype(seurat_obj,
                                            method = "quantile",
                                            verbose = verbose)
  }
  
  if (!is.null(cell_type_col) && cell_type_col %in% colnames(seurat_obj@meta.data)) {
    if (verbose) {
      message("\n[5/5] Aggregating by cell type...")
    }
    
    cell_type_stats <- aggregateByCellType(seurat_obj, 
                                          cell_type_col = cell_type_col,
                                          include_sd = TRUE,
                                          include_sem = TRUE)
    
    phenotype_stats <- NULL
    if ("metabolic_phenotype" %in% colnames(seurat_obj@meta.data)) {
      phenotype_stats <- getPhenotypeStatistics(seurat_obj)
    }
    
    seurat_obj@misc$metabolic_summary <- list(
      cell_type_stats = cell_type_stats,
      phenotype_stats = phenotype_stats,
      analysis_params = list(
        pathways = names(gene_sets),
        scoring_method = scoring_method,
        organism = organism,
        timestamp = Sys.time()
      )
    )
  }
  
  if (verbose) {
    message("\n========================================")
    message("[scMetaboFlux] Analysis Complete!")
    message("========================================")
    message("\nAdded metadata columns:")
    new_cols <- grep("^(glycolysis|oxidative|ATP|GOX|flux|metabolic|contribution)",
                    colnames(seurat_obj@meta.data), value = TRUE)
    message(paste(head(new_cols, 10), collapse = "\n  "))
    if (length(new_cols) > 10) {
      message(paste0("  ... and ", length(new_cols) - 10, " more"))
    }
  }
  
  return(seurat_obj)
}

#' @title Create scMetaboFlux Object
#' 
#' @description Create an scMetaboFlux S4 object from Seurat for
#' structured storage and retrieval of results.
#' 
#' @param seurat_obj Seurat object
#' @param name Optional name for the object
#' 
#' @return scMetaboFluxObject
#' 
#' @export
#' @examples
#' \dontrun{
#' smf_obj <- createScMetaboFlux(seurat_obj, name = "my_analysis")
#' }
createScMetaboFlux <- function(seurat_obj, name = NULL) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object")
  }
  
  obj <- scMetaboFluxObject(
    seurat = seurat_obj,
    name = ifelse(is.null(name), "scMetaboFlux_Analysis", name),
    version = "1.0.0",
    results = list()
  )
  
  return(obj)
}

#' @title Set S4 Class Definition
#' @rdname scMetaboFluxObject-class
setClass("scMetaboFluxObject",
         representation(
           seurat = "ANY",
           name = "character",
           version = "character",
           results = "list"
         ))

#' @title Initialize scMetaboFlux Object
#' @rdname initialize-scMetaboFluxObject-method
#' @param .Object A scMetaboFlux object
#' @param ... Additional arguments
setMethod("initialize", "scMetaboFluxObject", function(.Object, ...) {
  .Object <- callNextMethod(.Object, ...)
  validObject(.Object)
  .Object
})

#' @title Validate scMetaboFlux Object
#' @name scMetaboFluxObject
setValidity("scMetaboFluxObject", function(object) {
  if (!inherits(object@seurat, "Seurat")) {
    return("seurat slot must contain a Seurat object")
  }
  if (length(object@name) == 0) {
    return("name slot must be non-empty")
  }
  return(TRUE)
})

#' @title Summary Method for scMetaboFlux
#' @rdname summary-scMetaboFluxObject-method
#' @param object A scMetaboFlux object
#' @export
setMethod("summary", "scMetaboFluxObject", function(object) {
  cat("scMetaboFlux Object:", object@name, "\n")
  cat("Version:", object@version, "\n")
  cat("Cells:", ncol(object@seurat), "\n")
  cat("Features:", nrow(object@seurat), "\n")
  
  new_cols <- grep("^(glycolysis|oxidative|ATP|GOX|flux|metabolic)",
                   colnames(object@seurat@meta.data), value = TRUE)
  cat("\nMetabolic metadata columns (", length(new_cols), "):\n")
  cat(paste("  -", head(new_cols, 15)), sep = "\n")
  if (length(new_cols) > 15) {
    cat("  ... and", length(new_cols) - 15, "more\n")
  }
  
  if (!is.null(object@seurat@misc$metabolic_summary)) {
    cat("\nAnalysis parameters:\n")
    params <- object@seurat@misc$metabolic_summary$analysis_params
    cat("  Pathways:", paste(params$pathways, collapse = ", "), "\n")
    cat("  Scoring method:", params$scoring_method, "\n")
  }
  
  invisible(NULL)
})

#' @title Print Method for scMetaboFlux
#' @rdname print-scMetaboFluxObject-method
#' @param x A scMetaboFlux object
#' @export
setMethod("print", "scMetaboFluxObject", function(x) {
  summary(x)
})

#' @title Extract Results from scMetaboFlux
#' 
#' @description Extract specific results from scMetaboFlux object.
#' 
#' @param object scMetaboFlux object
#' @param what What to extract: "scores", "phenotypes", "flux", "all"
#' 
#' @return Requested data
#' 
#' @export
#' @examples
#' \dontrun{
#' smf_obj <- createScMetaboFlux(seurat_obj)
#' scores <- extractResults(smf_obj, what = "scores")
#' }
#' @export
setGeneric("extractResults", function(object, what) standardGeneric("extractResults"))

#' @rdname extractResults
setMethod("extractResults", "scMetaboFluxObject", function(object, what = c("scores", "phenotypes", "flux", "all")) {
  
  what <- match.arg(what)
  
  seurat_obj <- object@seurat
  
  if (what == "scores") {
    score_cols <- grep("^(glycolysis|oxidative|ATP|GOX)",
                      colnames(seurat_obj@meta.data), value = TRUE)
    return(seurat_obj@meta.data[, score_cols, drop = FALSE])
  } else if (what == "phenotypes") {
    if ("metabolic_phenotype" %in% colnames(seurat_obj@meta.data)) {
      return(seurat_obj@meta.data$metabolic_phenotype)
    }
    return(NULL)
  } else if (what == "flux") {
    flux_cols <- grep("_flux$", colnames(seurat_obj@meta.data), value = TRUE)
    return(seurat_obj@meta.data[, flux_cols, drop = FALSE])
  } else {
    return(seurat_obj@meta.data)
  }
})

#' @title Check if Scores Exist
#' 
#' @description Check if metabolic scores have been computed.
#' 
#' @param seurat_obj Seurat object
#' 
#' @return Logical indicating if scores exist
#' 
#' @export
#' @examples
#' \dontrun{
#' hasMetabolicScores(seurat_obj)
#' }
hasMetabolicScores <- function(seurat_obj) {
  score_cols <- c("glycolysis", "oxidative_phosphorylation", "ATP_score")
  return(all(score_cols %in% colnames(seurat_obj@meta.data)))
}

#' @title Get Analysis Parameters
#' 
#' @description Get parameters used for analysis.
#' 
#' @param seurat_obj Seurat object
#' 
#' @return List of analysis parameters or NULL
#' 
#' @export
#' @examples
#' \dontrun{
#' params <- getAnalysisParameters(seurat_obj)
#' }
getAnalysisParameters <- function(seurat_obj) {
  if (!is.null(seurat_obj@misc$metabolic_summary$analysis_params)) {
    return(seurat_obj@misc$metabolic_summary$analysis_params)
  }
  return(NULL)
}

#' @title Reset Metabolic Results
#' 
#' @description Remove all scMetaboFlux results from Seurat object.
#' 
#' @param seurat_obj Seurat object
#' @param keep_original Keep original scores if present
#' 
#' @return Seurat object with results removed
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- resetMetabolicResults(seurat_obj)
#' }
resetMetabolicResults <- function(seurat_obj, keep_original = FALSE) {
  
  meta_cols <- colnames(seurat_obj@meta.data)
  
  metabolic_cols <- c(
    grep("^(glycolysis|oxidative_phosphorylation|tca_cycle|ATP|GOX|flux|metabolic|contribution|dysregulation)",
         meta_cols, value = TRUE)
  )
  
  if (keep_original && "glycolysis" %in% meta_cols) {
    metabolic_cols <- setdiff(metabolic_cols, 
                            c("glycolysis", "oxidative_phosphorylation", "ATP_score"))
  }
  
  if (length(metabolic_cols) > 0) {
    seurat_obj@meta.data[, metabolic_cols] <- NULL
  }
  
  seurat_obj@misc$metabolic_summary <- NULL
  
  return(seurat_obj)
}

#' @title Export Results to File
#' 
#' @description Export metabolic analysis results to CSV files.
#' 
#' @param seurat_obj Seurat object
#' @param output_dir Output directory
#' @param include_embeddings Include UMAP/tSNE embeddings
#' 
#' @export
#' @examples
#' \dontrun{
#' exportMetabolicResults(seurat_obj, 
#'   output_dir = "scMetaboFlux_results",
#'   include_embeddings = TRUE)
#' }
exportMetabolicResults <- function(seurat_obj,
                                  output_dir = "scMetaboFlux_results",
                                  include_embeddings = TRUE) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  score_cols <- grep("^(glycolysis|oxidative|ATP|GOX|flux|metabolic|contribution)",
                    colnames(seurat_obj@meta.data), value = TRUE)
  
  scores_df <- seurat_obj@meta.data[, c(score_cols), drop = FALSE]
  scores_df <- cbind(cell_id = rownames(scores_df), scores_df)
  rownames(scores_df) <- NULL
  
  write.csv(scores_df, 
           file.path(output_dir, "metabolic_scores.csv"),
           row.names = FALSE)
  
  if ("metabolic_phenotype" %in% colnames(seurat_obj@meta.data)) {
    phenotype_summary <- getPhenotypeDistribution(seurat_obj)
    write.csv(phenotype_summary, 
             file.path(output_dir, "phenotype_distribution.csv"),
             row.names = FALSE)
  }
  
  if (include_embeddings && !is.null(seurat_obj@reductions)) {
    for (reduction in names(seurat_obj@reductions)) {
      emb <- Seurat::Embeddings(seurat_obj, reduction = reduction)
      emb_df <- as.data.frame(emb)
      emb_df <- cbind(cell_id = rownames(emb_df), emb_df)
      rownames(emb_df) <- NULL
      
      write.csv(emb_df,
               file.path(output_dir, paste0("embeddings_", reduction, ".csv")),
               row.names = FALSE)
    }
  }
  
  message("[scMetaboFlux] Results exported to: ", output_dir)
  
  return(invisible(output_dir))
}

#' @title Validate Input Data
#' 
#' @description Validate Seurat object for scMetaboFlux analysis.
#' 
#' @param seurat_obj Seurat object
#' @param check_genes Check for metabolic genes
#' 
#' @return List with validation results containing:
#'   \item{valid}{Logical, TRUE if valid}
#'   \item{errors}{Character vector of errors}
#'   \item{warnings}{Character vector of warnings}
#' 
#' @export
#' @examples
#' \dontrun{
#' validation <- validateInputData(seurat_obj)
#' if (!validation$valid) {
#'   print(validation$errors)
#' }
#' }
validateInputData <- function(seurat_obj, check_genes = TRUE) {
  
  errors <- c()
  warnings <- c()
  
  if (!inherits(seurat_obj, "Seurat")) {
    errors <- c(errors, "Input must be a Seurat object")
    return(list(valid = FALSE, errors = errors, warnings = warnings))
  }
  
  if (ncol(seurat_obj) == 0) {
    errors <- c(errors, "Seurat object contains no cells")
  }
  
  if (nrow(seurat_obj) == 0) {
    errors <- c(errors, "Seurat object contains no features")
  }
  
  if (is.null(Seurat::DefaultAssay(seurat_obj))) {
    errors <- c(errors, "No default assay set")
  }
  
  expr_data <- Seurat::GetAssayData(seurat_obj)
  if (nrow(expr_data) > 0 && ncol(expr_data) > 0) {
    if (all(Matrix::rowSums(expr_data) == 0)) {
      errors <- c(errors, "Expression matrix contains only zeros")
    }
  } else {
    errors <- c(errors, "Expression matrix is empty")
  }
  
  if (check_genes && nrow(seurat_obj) > 0) {
    genes <- rownames(seurat_obj)
    genes_upper <- toupper(genes)
    
    mapped_count <- sum(sapply(metabolicGeneSets, function(gs) {
      sum(toupper(gs) %in% genes_upper)
    }))
    
    unique_metabolic_genes <- length(unique(unlist(lapply(metabolicGeneSets, toupper))))
    
    coverage_pct <- mapped_count / unique_metabolic_genes * 100
    
    if (coverage_pct < 5) {
      warnings <- c(warnings, 
                   sprintf("Only %.1f%% of metabolic genes found in dataset", coverage_pct))
    }
    
    if (mapped_count < 10) {
      warnings <- c(warnings, 
                   sprintf("Only %d metabolic genes found - may affect accuracy", mapped_count))
    }
  }
  
  return(list(
    valid = length(errors) == 0,
    errors = errors,
    warnings = warnings
  ))
}
