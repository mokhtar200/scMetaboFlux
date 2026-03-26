# scMetaboFlux NEWS

## scMetaboFlux 1.0.0 (2024-01-15)

### Major Features

* Initial release to Bioconductor
* Comprehensive single-cell energy metabolism inference from transcriptomics
* Novel ATP production estimation algorithm
* Metabolic flux approximation with rate-limiting enzyme weighting
* Classification of cells into 5 metabolic phenotypes (Glycolytic, Oxidative, Energetically Balanced, Energy-Stressed, Hypermetabolic)
* Support for multiple pathway scoring methods (AUCell, GSVA, mean expression, weighted mean)
* Integration with Seurat and SingleCellExperiment workflows

### Metabolic Pathways Supported

* Glycolysis and gluconeogenesis
* TCA cycle
* Oxidative phosphorylation (Complex I-V)
* Fatty acid oxidation
* Pentose phosphate pathway
* Glutaminolysis
* Pyruvate metabolism
* Electron transport chain
* Mitochondrial biogenesis
* Hypoxia response

### Key Functions

* `runMetabolicAnalysis()` - Complete workflow in one function
* `computePathwayScores()` - Pathway activity scoring
* `estimateATPProduction()` - ATP estimation from transcriptomics
* `computeMetabolicFlux()` - Pseudo-flux calculation
* `classifyMetabolicPhenotype()` - Metabolic phenotype classification
* `differentialMetabolicAnalysis()` - Disease vs control comparison
* `aggregateByCellType()` - Cell-type level aggregation
* `plotMetabolicUMAP()` - Visualization functions

### Documentation

* Full vignette with step-by-step analysis workflow
* 70+ exported functions with documentation
* Comprehensive test suite (150+ tests)
* Example analysis scripts

### Performance

* Parallel processing support for pathway scoring
* Memory-efficient sparse matrix handling
* Tested on datasets up to 500,000 cells
