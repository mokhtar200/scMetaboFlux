#' scMetaboFlux Example Analysis Script
#' 
#' This script demonstrates a complete analysis workflow using scMetaboFlux.
#' 
#' To run this example:
#' 1. Load scMetaboFlux: library(scMetaboFlux)
#' 2. Have a Seurat object with scRNA-seq data ready
#' 3. Replace "your_seurat_object" with your actual object

# Load the package
library(scMetaboFlux)
library(Seurat)

# ============================================================================
# EXAMPLE 1: QUICK ANALYSIS WITH DEFAULTS
# ============================================================================

# If you have a Seurat object named 'seurat_obj':
# seurat_obj <- runMetabolicAnalysis(seurat_obj,
#                                   pathways_to_analyze = c("glycolysis", 
#                                                           "tca_cycle",
#                                                           "oxidative_phosphorylation"),
#                                   scoring_method = "AUCell",
#                                   classify_phenotype = TRUE)

# ============================================================================
# EXAMPLE 2: STEP-BY-STEP ANALYSIS
# ============================================================================

# Step 1: Load and preprocess your data
# ------------------------------------
# seurat_obj <- readRDS("path/to/your/data.rds")
# seurat_obj <- NormalizeData(seurat_obj)
# seurat_obj <- FindVariableFeatures(seurat_obj, nfeatures = 2000)
# seurat_obj <- ScaleData(seurat_obj)
# seurat_obj <- RunPCA(seurat_obj, npcs = 50)
# seurat_obj <- RunUMAP(seurat_obj, dims = 1:50)
# seurat_obj <- FindNeighbors(seurat_obj, dims = 1:50)
# seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# Step 2: Map genes to pathways
# -----------------------------
# Get statistics about metabolic gene coverage
# stats <- getMetabolicGeneSetStats(seurat_obj)
# print(stats)

# Step 3: Compute pathway scores
# ------------------------------
# seurat_obj <- computePathwayScores(seurat_obj,
#                                   gene_sets = metabolicGeneSets[c("glycolysis", 
#                                                                  "tca_cycle",
#                                                                  "oxidative_phosphorylation",
#                                                                  "fatty_acid_oxidation")],
#                                   method = "AUCell",
#                                   parallel_cores = 4)

# Step 4: Estimate ATP production
# ------------------------------
# seurat_obj <- estimateATPProduction(seurat_obj,
#                                    normalize = TRUE,
#                                    verbose = TRUE)

# Step 5: Calculate metabolic ratios
# -----------------------------------
# seurat_obj <- calculateATPYieldRatio(seurat_obj)
# seurat_obj <- calculateGOXIndex(seurat_obj, z_score = TRUE)

# Step 6: Compute metabolic flux
# ------------------------------
# seurat_obj <- computeMetabolicFlux(seurat_obj,
#                                   gene_sets = metabolicGeneSets,
#                                   use_rate_limiting = TRUE,
#                                   normalize = TRUE)
# seurat_obj <- calculateFluxContribution(seurat_obj)

# Step 7: Classify metabolic phenotypes
# -------------------------------------
# seurat_obj <- classifyMetabolicPhenotype(seurat_obj,
#                                         method = "quantile",
#                                         verbose = TRUE)

# Step 8: Aggregate by cell type
# ------------------------------
# if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
#   cell_type_stats <- aggregateByCellType(seurat_obj, 
#                                         cell_type_col = "cell_type")
#   print(cell_type_stats)
# }

# ============================================================================
# EXAMPLE 3: DIFFERENTIAL METABOLIC ANALYSIS
# ============================================================================

# Compare metabolic states between conditions
# diff_results <- differentialMetabolicAnalysis(seurat_obj,
#                                             condition_col = "condition",
#                                             control_group = "Control",
#                                             case_group = "Disease",
#                                             method = "wilcox.test")
# print(diff_results$results)

# Pairwise comparison across multiple groups
# pairwise_results <- pairwiseMetabolicComparison(seurat_obj,
#                                               group_col = "treatment",
#                                               method = "t.test")
# print(pairwise_results)

# ============================================================================
# EXAMPLE 4: VISUALIZATION
# ============================================================================

# Plot metabolic scores on UMAP
# p1 <- plotMetabolicUMAP(seurat_obj, 
#                        score_cols = c("glycolysis", "oxidative_phosphorylation", "ATP_score"),
#                        ncol = 3)
# print(p1)

# Plot phenotype distribution
# p2 <- plotPhenotypeDistribution(seurat_obj, 
#                                group_col = "cell_type")
# print(p2)

# Plot glycolysis vs OXPHOS
# p3 <- plotGlycolysisVsOxphos(seurat_obj, 
#                             color_by = "metabolic_phenotype",
#                             density = TRUE)
# print(p3)

# Plot ATP distribution
# p4 <- plotATPDistribution(seurat_obj, 
#                          group_col = "cell_type",
#                          type = "violin")
# print(p4)

# Plot pathway correlation
# p5 <- plotPathwayCorrelation(seurat_obj)
# print(p5)

# Create comprehensive dashboard
# dashboard <- createMetabolicDashboard(seurat_obj, cell_type_col = "cell_type")
# print(dashboard)

# ============================================================================
# EXAMPLE 5: EXPORT RESULTS
# ============================================================================

# Export all results
# exportMetabolicResults(seurat_obj, 
#                      output_dir = "scMetaboFlux_results",
#                      include_embeddings = TRUE)

# Export cell type summary
# if ("cell_type" %in% colnames(seurat_obj@meta.data)) {
#   summary <- exportCellTypeSummary(seurat_obj,
#                                   cell_type_col = "cell_type",
#                                   output_dir = "cell_type_summary")
# }

# ============================================================================
# EXAMPLE 6: CUSTOM ANALYSIS
# ============================================================================

# Use custom gene sets
# custom_geneset <- createCustomGeneSet(
#   genes = c("GENE1", "GENE2", "GENE3"),
#   name = "my_pathway"
# )

# Create custom gene set collection
# custom_gsc <- createMetabolicGeneSetCollection(
#   gene_sets = list(
#     my_custom_pathway = c("GENE1", "GENE2"),
#     glycolysis = metabolicGeneSets$glycolysis
#   )
# )

# Compare scoring methods
# method_comparison <- comparePathwayMethods(seurat_obj,
#                                          gene_sets = metabolicGeneSets[1:3],
#                                          methods = c("AUCell", "GSVA", "mean"))
# print(method_comparison)

# ============================================================================
# EXAMPLE 7: VALIDATION AND QUALITY CONTROL
# ============================================================================

# Validate input data
# validation <- validateInputData(seurat_obj, check_genes = TRUE)
# if (!validation$valid) {
#   print(validation$errors)
# }
# if (length(validation$warnings) > 0) {
#   print(validation$warnings)
# }

# Check if analysis has been run
# hasMetabolicScores(seurat_obj)  # Returns TRUE if scores exist

# Get analysis parameters
# params <- getAnalysisParameters(seurat_obj)

# ============================================================================
# EXAMPLE 8: MULTI-OMICS INTEGRATION
# ============================================================================

# Integrate spatial data (requires spatial coordinates)
# spatial_coords <- read.csv("spatial_coordinates.csv")
# seurat_obj <- integrateSpatialMetabolism(seurat_obj, 
#                                         spatial_coords = spatial_coords)

# Create multi-omics landscape
# landscape <- createMultiOmicsMetabolicLandscape(seurat_obj,
#                                               include_spatial = FALSE,
#                                               include_proteomics = FALSE)

# ============================================================================
# NOTES AND TIPS
# ============================================================================

# 1. Minimum requirements:
#    - Seurat object with normalized expression data
#    - At least 50 cells recommended
#    - Metabolic genes detected (automatically checked)

# 2. For best results:
#    - Use SCTransform normalization for better detection
#    - Run UMAP/PCA before analysis for visualizations
#    - Include cell type annotations for cell-type level analysis

# 3. Troubleshooting:
#    - If AUCell/GSVA fails, use method = "mean"
#    - If memory issues, reduce parallel_cores
#    - Check gene coverage with getMetabolicGeneSetStats()

# 4. Interpreting GOX Index:
#    - GOX > 0: Glycolytic phenotype
#    - GOX < 0: Oxidative phenotype
#    - GOX ≈ 0: Balanced metabolism

# 5. Phenotype meanings:
#    - Glycolytic: Cancer-like, high glucose fermentation
#    - Oxidative: Mitochondrial-dependent, high ATP efficiency
#    - Energetically Balanced: Normal metabolic state
#    - Energy-Stressed: Low ATP, potential dysfunction
#    - Hypermetabolic: High energy demand, active metabolism

# ============================================================================
# CITATION
# ============================================================================

# citation()  # Print citation information
