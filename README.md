# scMetaboFlux: Single-Cell Energy Metabolism Inference Engine

[![Release](https://www.bioconductor.org/images/release_badge_2024_short.svg)](https://www.bioconductor.org/packages/release/bioc/html/scMetaboFlux.html)
[![Build Status](https://travis-ci.com/scMetaboFlux/scMetaboFlux.svg?branch=master)](https://travis-ci.com/scMetaboFlux/scMetaboFlux)
[![Coverage](https://codecov.io/gh/scMetaboFlux/scMetaboFlux/branch/master/graph/badge.svg)](https://codecov.io/gh/scMetaboFlux/scMetaboFlux)

## Overview

scMetaboFlux is a Bioconductor package for inferring cellular energy states from single-cell transcriptomics data via metabolic flux approximation. It converts gene expression data into quantitative metabolic estimates including ATP production, glycolysis activity, OXPHOS activity, and metabolic phenotype classification at single-cell resolution.

## Installation

### From Bioconductor (recommended)

```r
if (!require("BiocManager"))
    install.packages("BiocManager")

BiocManager::install("scMetaboFlux")
```

### Development version from GitHub

```r
# Install from GitHub
if (!require("devtools"))
    install.packages("devtools")
devtools::install_github("scMetaboFlux/scMetaboFlux")
```

## Quick Start

```r
library(scMetaboFlux)

# Load your Seurat object
seurat_obj <- readRDS("your_data.rds")

# Run complete metabolic analysis
seurat_obj <- runMetabolicAnalysis(seurat_obj,
                                  pathways_to_analyze = c("glycolysis", "tca_cycle", 
                                                          "oxidative_phosphorylation"),
                                  scoring_method = "AUCell",
                                  classify_phenotype = TRUE,
                                  cell_type_col = "cell_type")
```

## Key Features

### 1. ATP Production Estimation (Core Novelty)

```
ATP_total = w1 × Glycolysis + w2 × OXPHOS + w3 × FAO
```

Where weights reflect relative ATP yield per pathway.

### 2. Metabolic Phenotype Classification

Classifies cells into 5 metabolic phenotypes:

| Phenotype | Description |
|-----------|-------------|
| Glycolytic | High glycolysis, low OXPHOS (Warburg-like) |
| Oxidative | High OXPHOS, low glycolysis |
| Energetically Balanced | Both pathways moderate |
| Energy-Stressed | Low overall ATP |
| Hypermetabolic | High ATP from both pathways |

### 3. Multiple Scoring Methods

- **AUCell**: AUC-based gene set enrichment
- **GSVA**: Gene Set Variation Analysis
- **Mean Expression**: Simple average
- **Weighted Mean**: Weighted by rate-limiting enzymes

### 4. Pre-built Gene Sets

| Pathway | Genes | Description |
|---------|-------|-------------|
| Glycolysis | 45 | Glucose breakdown to pyruvate |
| TCA Cycle | 52 | Citric acid cycle enzymes |
| OXPHOS | 150+ | Electron transport chain complexes |
| Fatty Acid Oxidation | 53 | Beta-oxidation pathway |
| Pentose Phosphate | 21 | PPP enzymes |
| Glutaminolysis | 14 | Glutamine metabolism |
| Hypoxia Response | 25 | HIF1A pathway genes |

## Usage Examples

### Complete Workflow

```r
library(scMetaboFlux)

# Load data
seurat_obj <- readRDS("pbmc_tumor.rds")

# Run analysis
seurat_obj <- runMetabolicAnalysis(seurat_obj,
                                  pathways_to_analyze = c("glycolysis", 
                                                          "tca_cycle",
                                                          "oxidative_phosphorylation",
                                                          "fatty_acid_oxidation"),
                                  scoring_method = "AUCell",
                                  classify_phenotype = TRUE,
                                  cell_type_col = "cell_type",
                                  condition_col = "condition")

# Visualize
plotMetabolicUMAP(seurat_obj, score_cols = "ATP_score")
plotPhenotypeDistribution(seurat_obj, group_col = "cell_type")
plotGlycolysisVsOxphos(seurat_obj, color_by = "metabolic_phenotype")
```

### Step-by-Step Analysis

```r
# 1. Compute pathway scores
seurat_obj <- computePathwayScores(seurat_obj,
                                   gene_sets = metabolicGeneSets[c("glycolysis", "oxphos")],
                                   method = "AUCell")

# 2. Estimate ATP production
seurat_obj <- estimateATPProduction(seurat_obj,
                                   normalize = TRUE)

# 3. Calculate GOX index
seurat_obj <- calculateGOXIndex(seurat_obj)

# 4. Classify phenotypes
seurat_obj <- classifyMetabolicPhenotype(seurat_obj)

# 5. Compare conditions
diff_results <- differentialMetabolicAnalysis(seurat_obj,
                                             condition_col = "condition",
                                             control_group = "Normal",
                                             case_group = "Tumor")
```

### Visualization

```r
# ATP score on UMAP
plotMetabolicUMAP(seurat_obj, score_cols = "ATP_score")

# Phenotype distribution
plotPhenotypeDistribution(seurat_obj, group_col = "cell_type")

# Glycolysis vs OXPHOS scatter
plotGlycolysisVsOxphos(seurat_obj, color_by = "cell_type")

# Pathway correlation heatmap
plotPathwayCorrelation(seurat_obj)

# Create dashboard
dashboard <- createMetabolicDashboard(seurat_obj, cell_type_col = "cell_type")
print(dashboard)
```

## Biological Interpretation

### GOX Index

- **Positive values**: Glycolytic phenotype
- **Negative values**: Oxidative phenotype
- **Near zero**: Balanced metabolism

### Validation

The package has been validated against:

1. **Known biology**: Cancer cells show glycolytic phenotype, neurons show oxidative
2. **External datasets**: Metabolomics data correlations
3. **Perturbation experiments**: Hypoxia-induced metabolic switching
4. **Literature**: Known metabolic gene expression patterns

## Functions

### Core Functions

| Function | Description |
|----------|-------------|
| `runMetabolicAnalysis()` | Complete workflow wrapper |
| `computePathwayScores()` | Pathway activity scoring |
| `estimateATPProduction()` | ATP estimation |
| `computeMetabolicFlux()` | Flux approximation |
| `classifyMetabolicPhenotype()` | Phenotype classification |

### Analysis Functions

| Function | Description |
|----------|-------------|
| `differentialMetabolicAnalysis()` | Disease vs control comparison |
| `compareCellTypesMetabolism()` | Cell type comparison |
| `aggregateByCellType()` | Cell-type level statistics |
| `calculateGOXIndex()` | Glycolytic-OXPHOS index |

### Visualization Functions

| Function | Description |
|----------|-------------|
| `plotMetabolicUMAP()` | UMAP overlay |
| `plotPhenotypeDistribution()` | Phenotype bar plots |
| `plotGlycolysisVsOxphos()` | Scatter plot |
| `plotATPDistribution()` | ATP violin/boxplot |
| `plotPathwayCorrelation()` | Correlation heatmap |
| `createMetabolicDashboard()` | Summary dashboard |

## Citation

```r
# Get citation
citation("scMetaboFlux")
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## Issues

Please report issues at: https://github.com/scMetaboFlux/scMetaboFlux/issues

## Related Packages

- [Seurat](https://satijalab.org/seurat/) - Single-cell RNA-seq analysis
- [AUCell](https://bioconductor.org/packages/AUCell) - AUC-based gene set enrichment
- [GSVA](https://bioconductor.org/packages/GSVA) - Gene Set Variation Analysis
- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) - S4 class for single-cell data


---

**scMetaboFlux**: Decoding cellular bioenergetics at single-cell resolution.
