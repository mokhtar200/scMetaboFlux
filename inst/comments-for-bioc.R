# scMetaboFlux Submission Comments for Bioconductor

## Bioconductor Package Submission

### Package: scMetaboFlux
### Version: 1.0.0
### Submission Date: 2024-01-15

## General Comments

This is an initial submission of scMetaboFlux to Bioconductor. The package provides
comprehensive tools for inferring cellular energy metabolism from single-cell transcriptomics
data.

## Key Features

1. **ATP Production Estimation**: Novel algorithm to estimate ATP production from gene
   expression data using metabolic flux approximation.

2. **Metabolic Phenotype Classification**: Classifies cells into 5 metabolic phenotypes
   based on glycolysis/OXPHOS activity patterns.

3. **Multiple Scoring Methods**: Supports AUCell, GSVA, mean expression, and weighted
   mean methods for pathway scoring.

4. **Comprehensive Gene Sets**: Pre-built metabolic gene sets from KEGG/Reactome for
   12+ pathways including glycolysis, TCA, OXPHOS, FAO, etc.

5. **Full Integration**: Works seamlessly with Seurat and SingleCellExperiment objects.

## Testing

- The package has been tested on multiple simulated datasets representing:
  * PBMC datasets (standard 10X format)
  * Cancer tumor microenvironment
  * Differentiation time courses
  * Hypoxia treatment experiments
  * Large-scale datasets (2000+ cells)

- All tests pass locally with the included testthat suite.

## Dependencies

All dependencies are from Bioconductor or standard R packages available on CRAN.
The package has been checked for:
- R CMD build
- R CMD check --no-manual --no-vignettes
- Windows compatibility
- macOS compatibility (via Bioconductor)

## Future Plans

- Integration with spatial transcriptomics (Visium, MERFISH)
- Support for multi-modal integration (CITE-seq, SHARE-seq)
- Web-based visualization application
- Shiny app for interactive analysis

## Contact

Maintainer: Metabolic Analysis Team <dev@scmetaboFlux.org>
