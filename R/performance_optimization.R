# ============================================================================
# PERFORMANCE OPTIMIZATION MODULE FOR scMetaboFlux
# ============================================================================
# High-performance utilities for large-scale single-cell analysis
# ============================================================================

#' @title Initialize Parallel Backend
#' 
#' @description Initialize parallel processing backend for scMetaboFlux.
#' 
#' @param n_cores Number of cores to use. If NULL, uses all available - 1.
#' @param backend Parallel backend: "multicore", "snow", or "future"
#' 
#' @return Invisible NULL
#' 
#' @export
#' @examples
#' \dontrun{
#' initParallelBackend(n_cores = 4)
#' }
initParallelBackend <- function(n_cores = NULL, backend = c("future", "multicore", "snow")) {
  
  backend <- match.arg(backend)
  
  if (is.null(n_cores)) {
    n_cores <- max(1, parallel::detectCores() - 1)
  }
  
  n_cores <- min(n_cores, parallel::detectCores())
  
  if (backend == "future") {
    if (!requireNamespace("future", quietly = TRUE)) {
      warning("future package not available, using sequential processing")
      return(invisible(NULL))
    }
    plan("multisession", workers = n_cores)
    message("[scMetaboFlux] Parallel backend initialized with ", n_cores, " cores (future)")
  } else if (backend == "multicore") {
    options(mc.cores = n_cores)
    message("[scMetaboFlux] Parallel backend initialized with ", n_cores, " cores (multicore)")
  } else {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      message("[scMetaboFlux] doParallel not available, using sequential processing")
      return(NULL)
    }
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    message("[scMetaboFlux] Parallel backend initialized with ", n_cores, " cores (snow)")
    return(cl)
  }
  
  invisible(NULL)
}

#' @title Scale Data with Progress
#' 
#' @description Memory-efficient scaling for large datasets.
#' 
#' @param expr_matrix Expression matrix (dgCMatrix)
#' @param block_size Number of genes to process at once
#' @param verbose Print progress
#' 
#' @return Scaled matrix
#' 
#' @export
scaleMatrixBlocks <- function(expr_matrix, block_size = 1000, verbose = TRUE) {
  
  n_genes <- nrow(expr_matrix)
  n_blocks <- ceiling(n_genes / block_size)
  
  scaled_matrix <- matrix(0, nrow = n_genes, ncol = ncol(expr_matrix))
  
  for (i in seq_len(n_blocks)) {
    start_idx <- (i - 1) * block_size + 1
    end_idx <- min(i * block_size, n_genes)
    
    block <- expr_matrix[start_idx:end_idx, , drop = FALSE]
    
    block_means <- Matrix::rowMeans(block)
    block_sds <- apply(block, 1, sd)
    block_sds[block_sds == 0] <- 1
    
    scaled_block <- (block - block_means) / block_sds
    scaled_matrix[start_idx:end_idx, ] <- as.matrix(scaled_block)
    
    if (verbose && i %% 10 == 0) {
      message("[scMetaboFlux] Scaled block ", i, "/", n_blocks)
    }
  }
  
  return(scaled_matrix)
}

#' @title Batch Process Large Dataset
#' 
#' @description Process large datasets in batches to manage memory.
#' 
#' @param seurat_obj Seurat object
#' @param batch_size Number of cells per batch
#' @param FUN Function to apply to each batch
#' @param ... Additional arguments to FUN
#' 
#' @return Seurat object with results
#' 
#' @export
#' @examples
#' \dontrun{
#' obj <- batchProcess(seurat_obj, batch_size = 5000, myFunction, arg1 = val)
#' }
batchProcess <- function(seurat_obj, batch_size = 5000, FUN, ...) {
  
  n_cells <- ncol(seurat_obj)
  n_batches <- ceiling(n_cells / batch_size)
  
  results_list <- vector("list", n_batches)
  
  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_cells)
    
    batch_indices <- start_idx:end_idx
    
    batch_obj <- seurat_obj[, batch_indices]
    
    batch_result <- FUN(batch_obj, ...)
    
    results_list[[i]] <- batch_result
    
    gc(verbose = FALSE)
  }
  
  return(results_list)
}

#' @title Memory-Efficient Pathway Scoring
#' 
#' @description Pathway scoring optimized for large datasets.
#' 
#' @param seurat_obj Seurat object
#' @param gene_sets Gene sets
#' @param batch_size Cells per batch
#' @param verbose Print progress
#' 
#' @return Matrix of pathway scores
#' 
#' @export
batchPathwayScores <- function(seurat_obj, gene_sets, batch_size = 5000, verbose = TRUE) {
  
  n_cells <- ncol(seurat_obj)
  n_batches <- ceiling(n_cells / batch_size)
  
  assay <- Seurat::DefaultAssay(seurat_obj)
  expr_matrix <- as(Seurat::GetAssayData(seurat_obj, assay = assay), "dgCMatrix")
  
  genes <- rownames(expr_matrix)
  genes_upper <- toupper(genes)
  
  pathway_names <- names(gene_sets)
  all_scores <- matrix(0, nrow = length(pathway_names), ncol = n_cells)
  
  for (batch_idx in seq_len(n_batches)) {
    start_cell <- (batch_idx - 1) * batch_size + 1
    end_cell <- min(batch_idx * batch_size, n_cells)
    
    batch_matrix <- expr_matrix[, start_cell:end_cell, drop = FALSE]
    
    for (pw_idx in seq_along(pathway_names)) {
      pathway <- pathway_names[pw_idx]
      pathway_genes <- toupper(gene_sets[[pathway]])
      
      matched_idx <- match(pathway_genes, genes_upper)
      matched_idx <- matched_idx[!is.na(matched_idx)]
      
      if (length(matched_idx) == 0) {
        all_scores[pw_idx, start_cell:end_cell] <- 0
        next
      }
      
      gene_expr <- batch_matrix[matched_idx, , drop = FALSE]
      scores <- Matrix::colMeans(gene_expr, na.rm = TRUE)
      all_scores[pw_idx, start_cell:end_cell] <- scores
    }
    
    if (verbose && batch_idx %% 5 == 0) {
      message("[scMetaboFlux] Processed batch ", batch_idx, "/", n_batches)
    }
    
    rm(batch_matrix)
    gc(verbose = FALSE)
  }
  
  rownames(all_scores) <- pathway_names
  colnames(all_scores) <- colnames(expr_matrix)
  
  return(all_scores)
}

#' @title Get Memory Usage
#' 
#' @description Get current memory usage statistics.
#' 
#' @return List with memory statistics
#' 
#' @export
getMemoryUsage <- function() {
  
  gc_info <- gc()
  
  memory_used_mb <- sum(gc_info[, 2]) / 1024
  memory_max_mb <- memory.limit() / 1024
  
  list(
    memory_used_mb = round(memory_used_mb, 2),
    memory_limit_mb = memory_max_mb,
    memory_usage_pct = round(memory_used_mb / memory_max_mb * 100, 2),
    n_garbage_collections = sum(gc_info[, 6])
  )
}

#' @title Clear Memory
#' 
#' @description Aggressively clear memory for large operations.
#' 
#' @export
clearMemory <- function() {
  gc(verbose = FALSE)
  invisible(NULL)
}

#' @title Efficient Correlation Matrix
#' 
#' @description Compute correlation matrix in a memory-efficient way.
#' 
#' @param data_matrix Data matrix
#' @param method Correlation method
#' @param batch_size Batch size for processing
#' 
#' @return Correlation matrix
#' 
#' @export
efficientCorMatrix <- function(data_matrix, method = c("pearson", "spearman"),
                               batch_size = 1000) {
  
  method <- match.arg(method)
  
  n_cols <- ncol(data_matrix)
  cor_matrix <- matrix(0, n_cols, n_cols)
  
  for (i in seq_len(ceiling(n_cols / batch_size))) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_cols)
    
    batch_data <- data_matrix[, start_idx:end_idx, drop = FALSE]
    cor_matrix[start_idx:end_idx, ] <- cor(batch_data, data_matrix, 
                                           method = method, use = "pairwise")
  }
  
  rownames(cor_matrix) <- colnames(data_matrix)
  colnames(cor_matrix) <- colnames(data_matrix)
  
  return(cor_matrix)
}

#' @title Sparse Matrix Statistics
#' 
#' @description Compute statistics on sparse matrices efficiently.
#' 
#' @param sparse_mat Sparse matrix
#' @param stats Stats to compute: "mean", "var", "nz_mean", "sparsity"
#' 
#' @return Named vector of statistics
#' 
#' @export
sparseMatrixStats <- function(sparse_mat, stats = c("mean", "nz_mean", "sparsity", "n_nonzero")) {
  
  results <- c()
  
  if ("mean" %in% stats) {
    results["mean"] <- mean(sparse_mat)
  }
  
  if ("nz_mean" %in% stats) {
    nz_vals <- sparse_mat[sparse_mat != 0]
    results["nz_mean"] <- if (length(nz_vals) > 0) mean(nz_vals) else 0
  }
  
  if ("sparsity" %in% stats) {
    results["sparsity"] <- sum(sparse_mat == 0) / length(sparse_mat)
  }
  
  if ("n_nonzero" %in% stats) {
    results["n_nonzero"] <- sum(sparse_mat != 0)
  }
  
  return(results)
}

#' @title Progress Bar Wrapper
#' 
#' @description Create progress bar wrapper with timing.
#' 
#' @param n Total iterations
#' @param style Progress bar style
#' 
#' @return Progress bar object
#' 
#' @export
createProgressBar <- function(n, style = 3) {
  pb <- txtProgressBar(min = 0, max = n, style = style)
  start_time <- Sys.time()
  
  list(
    pb = pb,
    start_time = start_time,
    update = function(i) {
      setTxtProgressBar(pb, i)
    },
    close = function() {
      close(pb)
      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      message(sprintf("Completed in %.1f seconds", elapsed))
    }
  )
}

#' @title Chunk Apply
#' 
#' @description Apply function to chunks of a vector/matrix in parallel.
#' 
#' @param data Data to chunk
#' @param chunk_size Size of each chunk
#' @param FUN Function to apply
#' @param n_cores Number of cores
#' @param ... Additional arguments to FUN
#' 
#' @return List of results
#' 
#' @export
chunkApply <- function(data, chunk_size, FUN, n_cores = 1, ...) {
  
  n <- length(data)
  n_chunks <- ceiling(n / chunk_size)
  
  chunks <- lapply(seq_len(n_chunks), function(i) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n)
    data[start_idx:end_idx]
  })
  
  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    
    results <- parallel::parLapply(cl, chunks, FUN, ...)
  } else {
    results <- lapply(chunks, FUN, ...)
  }
  
  return(results)
}

#' @title Optimize Seurat Object
#' 
#' @description Optimize Seurat object for memory efficiency.
#' 
#' @param seurat_obj Seurat object
#' @param remove_scale_data Remove scale.data to save memory
#' @param compress Compress assay data
#' 
#' @return Optimized Seurat object
#' 
#' @export
optimizeSeuratObject <- function(seurat_obj, 
                                  remove_scale_data = TRUE,
                                  compress = TRUE) {
  
  if (remove_scale_data) {
    assays <- Seurat::Assays(seurat_obj)
    for (assay_name in assays) {
      if (!is.null(seurat_obj[[assay_name]]@scale.data)) {
        seurat_obj[[assay_name]]@scale.data <- matrix()
      }
    }
  }
  
  gc(verbose = FALSE)
  
  return(seurat_obj)
}

#' @title Benchmark Function
#' 
#' @description Benchmark a function with timing and memory tracking.
#' 
#' @param FUN Function to benchmark
#' @param n_reps Number of repetitions
#' @param ... Arguments to FUN
#' 
#' @return List with benchmark results
#' 
#' @export
benchmarkFunction <- function(FUN, n_reps = 3, ...) {
  
  mem_before <- getMemoryUsage()
  
  times <- numeric(n_reps)
  results <- list()
  
  for (i in seq_len(n_reps)) {
    gc()
    start_time <- Sys.time()
    
    results[[i]] <- FUN(...)
    
    end_time <- Sys.time()
    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  mem_after <- getMemoryUsage()
  
  list(
    mean_time = mean(times),
    min_time = min(times),
    max_time = max(times),
    sd_time = sd(times),
    memory_delta_mb = mem_after$memory_used_mb - mem_before$memory_used_mb,
    results = results[[which.min(times)]]
  )
}

#' @title Cached Pathway Scoring
#' 
#' @description Pathway scoring with caching for repeated analyses.
#' 
#' @param seurat_obj Seurat object
#' @param cache_file Path to cache file
#' @param gene_sets Gene sets
#' @param method Scoring method
#' 
#' @return Matrix of scores
#' 
#' @export
cachedPathwayScores <- function(seurat_obj, cache_file, gene_sets, method = "mean") {
  
  if (file.exists(cache_file)) {
    message("[scMetaboFlux] Loading cached scores from ", cache_file)
    return(readRDS(cache_file))
  }
  
  message("[scMetaboFlux] Computing scores...")
  scores <- computePathwayScores(seurat_obj, gene_sets = gene_sets, 
                                 method = method, verbose = FALSE)
  
  dir.create(dirname(cache_file), showWarnings = FALSE, recursive = TRUE)
  saveRDS(scores, cache_file)
  message("[scMetaboFlux] Scores cached to ", cache_file)
  
  return(scores)
}

#' @title Get Optimal Batch Size
#' 
#' @description Calculate optimal batch size based on available memory.
#' 
#' @param seurat_obj Seurat object
#' @param target_memory_pct Target memory usage percentage
#' 
#' @return Recommended batch size
#' 
#' @export
getOptimalBatchSize <- function(seurat_obj, target_memory_pct = 50) {
  
  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj)
  
  available_memory_mb <- memory.limit() / 1024
  current_usage_mb <- sum(gc()[, 2]) / 1024
  
  memory_per_cell_mb <- (current_usage_mb * 1024) / n_cells
  
  target_memory_mb <- available_memory_mb * (target_memory_pct / 100)
  usable_memory_mb <- target_memory_mb - current_usage_mb
  
  optimal_batch_size <- floor(usable_memory_mb / memory_per_cell_mb)
  
  optimal_batch_size <- max(100, min(optimal_batch_size, n_cells))
  
  return(optimal_batch_size)
}
