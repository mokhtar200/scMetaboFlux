# ============================================================================
# CORE TESTS FOR scMetaboFlux
# ============================================================================
# Minimal test suite that validates core functionality

library(testthat)
library(scMetaboFlux)

test_that("metabolicGeneSets is properly defined", {
  expect_type(metabolicGeneSets, "list")
  expect_true(length(metabolicGeneSets) > 0)
  expect_true("glycolysis" %in% names(metabolicGeneSets))
  expect_true("oxidative_phosphorylation" %in% names(metabolicGeneSets))
})

test_that("rateLimitingWeights is properly defined", {
  expect_type(rateLimitingWeights, "double")
  expect_true(length(rateLimitingWeights) > 0)
})

test_that("metabolicGeneSets contain valid human genes", {
  expect_true("HK2" %in% metabolicGeneSets$glycolysis)
  expect_true("LDHA" %in% metabolicGeneSets$glycolysis)
  expect_true("NDUFA9" %in% metabolicGeneSets$oxidative_phosphorylation)
  expect_true("CS" %in% metabolicGeneSets$tca_cycle)
})

test_that("ATP yield coefficients are defined", {
  expect_true(is.atomic(atpYieldCoefficients) || is.list(atpYieldCoefficients))
  expect_true("glycolysis" %in% names(atpYieldCoefficients))
  expect_true("oxidative_phosphorylation" %in% names(atpYieldCoefficients))
})

test_that("validateInputData works", {
  # Create minimal test data
  set.seed(42)
  n_cells <- 50
  n_genes <- 100
  
  # Include metabolic genes in the matrix
  metabolic_genes <- c("HK2", "LDHA", "NDUFA9", "SDHA", "CS", "IDH1", "PKM2")
  other_genes <- paste0("GENE", 1:(n_genes - length(metabolic_genes)))
  all_genes <- c(metabolic_genes, other_genes)
  
  # Create expression matrix
  expr_matrix <- Matrix::Matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 5),
    nrow = n_genes,
    ncol = n_cells,
    sparse = TRUE,
    dimnames = list(all_genes, paste0("Cell_", 1:n_cells))
  )
  
  # Add higher expression for metabolic genes
  for (gene in metabolic_genes) {
    expr_matrix[gene, ] <- rpois(n_cells, lambda = 10)
  }
  
  # Create Seurat object
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  
  # Test validation
  result <- validateInputData(seurat_obj)
  expect_type(result, "list")
  expect_true("valid" %in% names(result))
})

test_that("computePathwayScores with mean method works", {
  set.seed(42)
  n_cells <- 50
  n_genes <- 100
  
  # Include actual metabolic genes
  glycolysis_genes <- c("HK2", "LDHA", "PKM2")
  oxphos_genes <- c("NDUFA9", "SDHA", "COX5A")
  other_genes <- paste0("GENE", 1:(n_genes - length(c(glycolysis_genes, oxphos_genes))))
  all_genes <- c(glycolysis_genes, oxphos_genes, other_genes)
  
  expr_matrix <- Matrix::Matrix(
    rnbinom(n_genes * n_cells, size = 2, mu = 5),
    nrow = n_genes,
    ncol = n_cells,
    sparse = TRUE,
    dimnames = list(all_genes, paste0("Cell_", 1:n_cells))
  )
  
  # Add higher expression for metabolic genes
  for (gene in glycolysis_genes) {
    expr_matrix[gene, ] <- rpois(n_cells, lambda = 15)
  }
  for (gene in oxphos_genes) {
    expr_matrix[gene, ] <- rpois(n_cells, lambda = 15)
  }
  
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  
  # Test pathway scoring
  result <- computePathwayScores(seurat_obj,
                                gene_sets = list(
                                  glycolysis = glycolysis_genes,
                                  oxphos = oxphos_genes
                                ),
                                method = "mean",
                                verbose = FALSE)
  
  expect_true("glycolysis" %in% colnames(result@meta.data))
  expect_true("oxphos" %in% colnames(result@meta.data))
  expect_true(all(!is.na(result@meta.data$glycolysis)))
})

test_that("estimateATPProduction works", {
  set.seed(42)
  n_cells <- 30
  
  # Create mock Seurat with required columns
  expr_matrix <- Matrix::Matrix(
    rnbinom(100 * n_cells, size = 2, mu = 5),
    nrow = 100,
    ncol = n_cells,
    sparse = TRUE
  )
  rownames(expr_matrix) <- paste0("GENE", 1:100)
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  seurat_obj@meta.data$glycolysis <- runif(n_cells, 0.3, 0.8)
  seurat_obj@meta.data$oxidative_phosphorylation <- runif(n_cells, 0.3, 0.8)
  
  result <- estimateATPProduction(seurat_obj, verbose = FALSE)
  
  expect_true("ATP_score" %in% colnames(result@meta.data))
  expect_true(all(is.finite(result@meta.data$ATP_score)))
})

test_that("calculateGOXIndex works", {
  set.seed(42)
  n_cells <- 30
  
  expr_matrix <- Matrix::Matrix(
    rnbinom(100 * n_cells, size = 2, mu = 5),
    nrow = 100,
    ncol = n_cells,
    sparse = TRUE
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  
  # Add required columns
  seurat_obj@meta.data$glycolysis <- runif(n_cells, 0.3, 0.9)
  seurat_obj@meta.data$oxidative_phosphorylation <- runif(n_cells, 0.1, 0.7)
  
  result <- calculateGOXIndex(seurat_obj, verbose = FALSE)
  
  expect_true("GOX_index" %in% colnames(result@meta.data))
  expect_true(is.numeric(result@meta.data$GOX_index))
})

test_that("classifyMetabolicPhenotype works", {
  set.seed(42)
  n_cells <- 50
  
  expr_matrix <- Matrix::Matrix(
    rnbinom(100 * n_cells, size = 2, mu = 5),
    nrow = 100,
    ncol = n_cells,
    sparse = TRUE
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  
  # Add required columns
  seurat_obj@meta.data$glycolysis <- runif(n_cells, 0.2, 0.9)
  seurat_obj@meta.data$oxidative_phosphorylation <- runif(n_cells, 0.2, 0.9)
  seurat_obj@meta.data$ATP_score <- runif(n_cells, 0.3, 0.8)
  
  result <- classifyMetabolicPhenotype(seurat_obj, verbose = FALSE)
  
  expect_true("metabolic_phenotype" %in% colnames(result@meta.data))
  expect_true(all(!is.na(result@meta.data$metabolic_phenotype)))
})

test_that("aggregateByCellType works", {
  set.seed(42)
  n_cells <- 100
  
  expr_matrix <- Matrix::Matrix(
    rnbinom(100 * n_cells, size = 2, mu = 5),
    nrow = 100,
    ncol = n_cells,
    sparse = TRUE
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  
  seurat_obj@meta.data$cell_type <- sample(c("T_cell", "B_cell", "Macrophage"), 
                                           n_cells, replace = TRUE)
  seurat_obj@meta.data$glycolysis <- runif(n_cells, 0.3, 0.8)
  seurat_obj@meta.data$ATP_score <- runif(n_cells, 0.3, 0.8)
  
  result <- aggregateByCellType(seurat_obj, cell_type_col = "cell_type", verbose = FALSE)
  
  expect_type(result, "list")
  expect_true("cell_type" %in% colnames(result))
  expect_true(nrow(result) >= 1)
})

test_that("validatePackageIntegrity works", {
  result <- validatePackageIntegrity()
  expect_type(result, "list")
  expect_true("status" %in% names(result))
})

test_that("zscorePathwayScores works", {
  set.seed(42)
  n_cells <- 30
  
  expr_matrix <- Matrix::Matrix(
    rnbinom(100 * n_cells, size = 2, mu = 5),
    nrow = 100,
    ncol = n_cells,
    sparse = TRUE
  )
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = expr_matrix, min.features = 10)
  seurat_obj@meta.data$glycolysis <- rnorm(n_cells)
  seurat_obj@meta.data$oxidative_phosphorylation <- rnorm(n_cells)
  
  result <- zscorePathwayScores(seurat_obj, score_cols = c("glycolysis", "oxidative_phosphorylation"))
  
  expect_true("zscore_glycolysis" %in% colnames(result@meta.data))
  expect_true("zscore_oxidative_phosphorylation" %in% colnames(result@meta.data))
})

test_that("mapGenesToPathways works", {
  result <- mapGenesToPathways(c("HK2", "LDHA", "NDUFA9", "CS"))
  expect_type(result, "list")
  expect_true("mapped_genes" %in% names(result))
  expect_true("coverage" %in% names(result))
})

cat("\n[scMetaboFlux] Core tests completed successfully!\n")
