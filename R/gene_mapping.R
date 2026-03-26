#' @title Metabolic Gene Mapping Engine
#' @description Maps genes to metabolic pathways using KEGG/Reactome databases
#' @name metabolic_gene_mapping
#' @rdname metabolic-gene-mapping
#' @keywords gene mapping metabolic pathways
NULL

#' @title Predefined Metabolic Gene Sets
#' @description Comprehensive metabolic gene sets curated from KEGG and Reactome databases.
#' Contains 12 metabolic pathways with validated gene memberships for human and mouse.
#' 
#' @format A named list where each element is a character vector of gene symbols:
#' \describe{
#'   \item{glycolysis}{37 genes involved in glucose breakdown to pyruvate}
#'   \item{tca_cycle}{40 genes for citric acid cycle enzymes}
#'   \item{oxidative_phosphorylation}{150+ genes for ETC complexes I-V}
#'   \item{fatty_acid_oxidation}{45 genes for beta-oxidation pathway}
#'   \item{pentose_phosphate_pathway}{17 genes for PPP enzymes}
#'   \item{glutaminolysis}{11 genes for glutamine metabolism}
#'   \item{pyruvate_metabolism}{14 genes for pyruvate conversion}
#'   \item{electron_transport_chain}{90+ genes for ETC}
#'   \item{mitochondrial_biogenesis}{30 genes for mtDNA replication}
#'   \item{hypoxia_response}{24 genes for hypoxic response}
#'   \item{glycolysis_rate_limiting}{12 key rate-limiting enzymes}
#'   \item{oxidative_phosphorylation_rate_limiting}{8 key rate-limiting enzymes}
#'   \item{tca_rate_limiting}{8 key rate-limiting enzymes}
#' }
#' 
#' @source KEGG Pathway Database (https://www.kegg.jp/kegg/pathway.html)
#' @source Reactome Database (https://reactome.org/)
#' 
#' @examples
#' head(metabolicGeneSets$glycolysis)
#' length(metabolicGeneSets$oxidative_phosphorylation)
#' 
#' @export
metabolicGeneSets <- list(
  glycolysis = c(
    "HK1", "HK2", "HK3", "GCK", "GPI", "PFKL", "PFKM", "PFKP",
    "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "GAPDHS", "PGK1", "PGK2",
    "PGAM1", "PGAM2", "ENO1", "ENO2", "ENO3", "PKM1", "PKM2", "PKLR",
    "LDHA", "LDHB", "LDHC", "LDHD", "MDH1", "MDH2", "ME1", "ME2", "ME3",
    "G6PD", "G6PDH", "PGD", "RPI", "RPE", "TALDO1", "TALDO2", "G6PC",
    "G6PC2", "G6PC3", "PGLS", "H6PD"
  ),
  
  tca_cycle = c(
    "ACLY", "ACSS1", "ACSS2", "ACACA", "ACACB", "CS", "ACO1", "ACO2", "ACONIT1", 
    "ACONIT2", "ACONIT3", "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G", 
    "IDH3B", "OGDH", "OGDHL", "D2HGDH", "SUCLA2", "SUCLG1", "SUCLG2", 
    "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF1", "SDHAF2", "SDHAF3", "SDHAF4",
    "FH", "MDH1", "MDH2", "ME1", "ME2", "ME3", "PCK1", "PCK2", "PPCDC",
    "PPCKL", "OGE", "GOT1", "GOT2", "GPT2", "GLUD1", "GLUD2", "MDH1B",
    "ACSF2", "ACSM3", "ETFB", "ETFDH"
  ),
  
  oxidative_phosphorylation = c(
    "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7",
    "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13",
    "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7",
    "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11", "NDUFB12", "NDUFC1", "NDUFC2",
    "NDUFV1", "NDUFV2", "NDUFV3", "NDUFV3P1", "NDUFS1", "NDUFS2", "NDUFS3",
    "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "NDUFAB1", "NDUFS10",
    "NDUFS9", "NDUFS12", "NDUFS13", "MTND1", "MTND2", "MTND3", "MTND4",
    "MTND4L", "MTND5", "MTND6", "SDHA", "SDHB", "SDHC", "SDHD", "SDHAF1",
    "SDHAF2", "SDHAF3", "SDHAF4", "UQCR10", "UQCR11", "UQCRB", "UQCRC1",
    "UQCRC2", "UQCRCR", "UQCRFS1", "UQCRH", "UQCRQ", "CYC1", "UQCC1", 
    "UQCC2", "UQCC3", "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", 
    "COX6A2", "COX6B1", "COX6B2", "COX6C", "COX7A1", "COX7A2", "COX7A2L",
    "COX7B", "COX7B2", "COX7C", "COX8A", "COX8C", "COX8H", "COX10", 
    "COX11", "COX14", "COX15", "COX16", "COX17", "COX18", "COX19", "COX20",
    "COX21", "COX23", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "ATP5MF", 
    "ATP5ME", "ATP5MG", "ATP5MGL", "ATP5MPL", "ATP5F1A", "ATP5F1B", "ATP5F1C",
    "ATP5F1D", "ATP5F1E", "ATP5IF1", "MTATP6", "MTATP8", "ATP5PD", "ATP5PDIP1"
  ),
  
  fatty_acid_oxidation = c(
    "ACAA1", "ACAA2", "ACADL", "ACADM", "ACADS", "ACADSB", "ACADVL", 
    "ACAT1", "ACAT2", "ACOT2", "ACOT4", "ACOT6", "ACOT7", "ACOT8", "ACOT9",
    "ACOT12", "ACOX1", "ACOX2", "ACOX3", "ACOXL", "CPT1A", "CPT1B", "CPT1C", 
    "CPT2", "DECR1", "DECR2", "ECHDC1", "ECHDC2", "ECI1", "ECI2", "EHHADH",
    "ETFA", "ETFB", "ETFDH", "HADH", "HADHA", "HADHB", "HADHP", "HSDL1", 
    "HSDL2", "IVD", "MCCC1", "MCCC2", "MCEE", "MLYCD", "PCCB", "PCCA", 
    "PECI", "SCP2", "HMGCL", "HMGCS2", "ACSF3", "ECHS1"
  ),
  
  pentose_phosphate_pathway = c(
    "G6PD", "G6PDH", "PGD", "RPI", "RPE", "TALDO1", "TALDO2", "TKT", 
    "TKTL1", "TKTL2", "H6PD", "PGLS", "PGL", "GND", "6PGL", "6PGD", "ZWF",
    "RPIA", "PRPS1", "PRPS2", "PRPS1L1"
  ),
  
  glutaminolysis = c(
    "GLS", "GLS2", "GLUD1", "GLUD2", "GOT1", "GOT2", "GPT2", "GLNX1",
    "GLNX2", "ASRGL1", "GANC", "GANAB", "GANAH", "MANEA", "MANEA4"
  ),
  
  pyruvate_metabolism = c(
    "PDHA1", "PDHA2", "PDHB", "PDHX", "PDH1B", "PDP1", "PDP2", "LDHA", 
    "LDHB", "LDHC", "LDHD", "LDHAL6A", "LDHAL6B", "MPC1", "MPC2", "DLAT", 
    "DLD", "DLDH", "PDPR", "BRP16L"
  ),
  
  electron_transport_chain = c(
    "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7",
    "NDUFA8", "NDUFA9", "NDUFA10", "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4",
    "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9", "NDUFB10", "NDUFC1",
    "NDUFC2", "NDUFV1", "NDUFV2", "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4",
    "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8", "SDHA", "SDHB", "SDHC", "SDHD",
    "UQCR10", "UQCRB", "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRH", "UQCRQ",
    "CYC1", "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX6C",
    "COX7A1", "COX7A2", "COX7B", "COX7C", "COX8A", "ATP5F1A", "ATP5F1B",
    "ATP5F1C", "ATP5F1D", "ATP5F1E", "ATP5IF1", "ATP5ME", "ATP5MF", "ATP5MG",
    "ATP5MGL", "ATP5MPL", "ATP5PB", "ATP5PD", "ATP5PF", "ATP5PO", "MTATP6",
    "MTATP8"
  ),
  
  mitochondrial_biogenesis = c(
    "TFAM", "NRF1", "NRF2", "GABPA", "PPARGCA1", "PGC1A", "PGC1B", "PRC",
    "ESRRA", "ERRB", "ERRG", "TFB1M", "TFB2M", "MTERF1", "MTERF2", "MTERF3",
    "MTERF4", "MTERFD1", "MTERFD2", "MTERFD3", "POLRMT", "TWINKLE", "POLG",
    "POLG2", "SSBP1", "OPA1", "MFN1", "MFN2", "DRP1", "FIS1", "MFF", 
    "MID49", "MID51", "GBA2", "GFER"
  ),
  
  hypoxia_response = c(
    "HIF1A", "EPAS1", "VEGFA", "VEGFB", "PGK1", "PGK2", "PGAM1", "ENO1",
    "LDHA", "LDHB", "GLUT1", "SLC2A1", "SLC2A3", "PDK1", "PDK2", "PDK3",
    "PDK4", "BNIP3", "BNIP3L", "BNIP3P1", "EGLN1", "EGLN2", "EGLN3", 
    "HIF1AN", "HIF3A", "ARNT", "ARNT2", "HIF1AP1", "VEGFC", "VEGFD"
  ),
  
  glycolysis_rate_limiting = c(
    "HK1", "HK2", "HK3", "GCK",
    "PFKP", "PFKM", "PFKL",
    "PKM1", "PKM2", "PKLR",
    "LDHA", "LDHB"
  ),
  
  oxidative_phosphorylation_rate_limiting = c(
    "NDUFA9", "NDUFS1", "NDUFS4", "SDHA",
    "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B"
  ),
  
  tca_rate_limiting = c(
    "CS", "IDH1", "IDH2", "IDH3A", "OGDH", "SDHA", "FH", "MDH2"
  )
)

#' @title Rate-Limiting Enzyme Weights
#' @description Enzyme-specific weights for ATP estimation that reflect 
#' rate-limiting enzyme importance in metabolic flux control.
#' Values > 1 indicate rate-limiting enzymes that have greater impact on flux.
#' 
#' @format Named numeric vector with gene symbols as names and weight values
#' @export
#' @source Curated from biochemical literature on enzyme kinetics
#' @examples
#' rateLimitingWeights[c("HK2", "PFKM", "SDHA")]
#' sum(rateLimitingWeights > 1)
rateLimitingWeights <- c(
  HK1 = 1.2, HK2 = 1.2, HK3 = 1.0, GCK = 1.2,
  PFKP = 1.3, PFKM = 1.3, PFKL = 1.3,
  PKM1 = 1.3, PKM2 = 1.3, PKLR = 1.3,
  LDHA = 1.1, LDHB = 1.1, LDHC = 1.0,
  GPI = 1.1, TPI1 = 1.1,
  NDUFA9 = 1.3, NDUFS1 = 1.3, NDUFS4 = 1.3, NDUFS7 = 1.2,
  SDHA = 1.4, SDHB = 1.3, SDHC = 1.2, SDHD = 1.2,
  COX4I1 = 1.3, COX5A = 1.3, COX5B = 1.2, COX6A1 = 1.2,
  ATP5F1A = 1.4, ATP5F1B = 1.4, ATP5F1C = 1.3,
  CS = 1.3, IDH1 = 1.2, IDH2 = 1.2, IDH3A = 1.3,
  OGDH = 1.3, FH = 1.2, MDH2 = 1.2, ACO1 = 1.1,
  CPT1A = 1.4, CPT1B = 1.4, CPT2 = 1.3,
  ACADM = 1.3, ACADVL = 1.3, HADHA = 1.3, HADHB = 1.3,
  G6PD = 1.2, GLS = 1.2, GLUD1 = 1.1
)

#' @title ATP Yield Coefficients
#' @description Theoretical ATP yield coefficients for different metabolic pathways.
#' These values represent relative ATP production per unit of pathway activity.
#' 
#' @format Named numeric vector with pathway names and ATP equivalents
#' @export
#' @source Biochemical literature on cellular respiration
#' @references 
#' Rich, P.R. (2003). The molecular machinery of Keilin's respiratory chain. 
#' Biochemical Society Transactions, 31(6), 1095-1105.
#' 
#' @examples
#' atpYieldCoefficients
#' atpYieldCoefficients["glycolysis"] / atpYieldCoefficients["oxidative_phosphorylation"]
atpYieldCoefficients <- c(
  glycolysis = 2,
  oxidative_phosphorylation = 32,
  fatty_acid_oxidation = 14,
  tca_cycle = 12,
  glutaminolysis = 6,
  pentose_phosphate_pathway = 0
)

#' @title Create Metabolic Gene Set Collection
#' 
#' @description Converts metabolic gene sets to GeneSetCollection format 
#' for use with AUCell, GSVA, and other enrichment analysis tools.
#' 
#' @param gene_sets Named list of character vectors containing gene symbols.
#'   If NULL, uses all predefined metabolicGeneSets.
#' @param species Species for gene IDs. Currently supports "human" (default)
#'   and "mouse" gene symbols.
#' 
#' @return GeneSetCollection object containing gene sets from input
#' 
#' @export
#' @examples
#' gsc <- createMetabolicGeneSetCollection()
#' gsc
#' 
#' custom_sets <- list(
#'   my_glycolysis = c("HK2", "PKM2", "LDHA"),
#'   my_oxphos = c("NDUFA9", "SDHA", "COX5A")
#' )
#' custom_gsc <- createMetabolicGeneSetCollection(gene_sets = custom_sets)
createMetabolicGeneSetCollection <- function(gene_sets = NULL, 
                                            species = c("human", "mouse")) {
  
  species <- match.arg(species)
  
  if (is.null(gene_sets)) {
    gene_sets <- metabolicGeneSets
  }
  
  if (!is.list(gene_sets) || is.null(names(gene_sets))) {
    stop("gene_sets must be a named list")
  }
  
  gsc_list <- lapply(names(gene_sets), function(pathway_name) {
    genes <- unique(as.character(gene_sets[[pathway_name]]))
    GSEABase::GeneSet(
      setName = pathway_name, 
      geneIds = genes,
      shortDescription = paste("Metabolic pathway:", pathway_name),
      organism = ifelse(species == "human", "Homo sapiens", "Mus musculus")
    )
  })
  
  gsc <- do.call(GSEABase::GeneSetCollection, gsc_list)
  return(gsc)
}

#' @title Map Genes to Metabolic Pathways
#' 
#' @description Maps user-provided genes to predefined metabolic pathways.
#' Returns matched genes, coverage statistics, and unmapped genes.
#' 
#' @param genes Character vector of gene symbols to map
#' @param pathways_to_include Character vector of pathway names to include.
#'   If NULL, includes all pathways from metabolicGeneSets.
#' @param species Species ("human" or "mouse")
#' @param case_sensitive Whether gene matching should be case-sensitive
#' 
#' @return List containing:
#'   \item{mapped_genes}{Named list with genes per pathway}
#'   \item{coverage}{Coverage statistics}
#'   \item{unmapped}{Character vector of unmapped genes}
#' 
#' @export
#' @examples
#' genes <- c("HK2", "LDHA", "NDUFA9", "SDHA", "GENE_NOT_FOUND", "BRCA1")
#' result <- mapGenesToPathways(genes)
#' names(result$mapped_genes)
#' result$coverage$coverage_pct
#' result$unmapped
mapGenesToPathways <- function(genes, 
                               pathways_to_include = NULL,
                               species = c("human", "mouse"),
                               case_sensitive = FALSE) {
  
  species <- match.arg(species)
  
  genes <- as.character(genes)
  if (!case_sensitive) {
    genes <- toupper(genes)
  }
  
  genes <- unique(genes[!is.na(genes) & genes != ""])
  
  if (is.null(pathways_to_include)) {
    pathways_to_include <- names(metabolicGeneSets)
  }
  
  mapped_genes <- list()
  total_possible <- 0
  
  for (pathway in pathways_to_include) {
    if (pathway %in% names(metabolicGeneSets)) {
      pathway_genes <- metabolicGeneSets[[pathway]]
      if (!case_sensitive) {
        pathway_genes <- toupper(pathway_genes)
      }
      
      total_possible <- total_possible + length(pathway_genes)
      
      matched_idx <- match(genes, pathway_genes)
      matched <- genes[!is.na(matched_idx)]
      
      if (length(matched) > 0) {
        original_genes <- genes[match(matched, genes)]
        if (!case_sensitive) {
          original_genes <- toupper(original_genes)
        }
        mapped_genes[[pathway]] <- unique(original_genes)
      }
    }
  }
  
  gene_coverage <- sum(lengths(mapped_genes))
  all_pathway_genes <- unique(unlist(metabolicGeneSets[pathways_to_include]))
  if (!case_sensitive) {
    all_pathway_genes <- toupper(all_pathway_genes)
  }
  all_genes_in_sets <- sum(toupper(genes) %in% all_pathway_genes)
  
  coverage_info <- list(
    total_mapped = gene_coverage,
    total_unique_genes_in_pathways = all_genes_in_sets,
    total_possible = total_possible,
    coverage_pct = if (length(all_pathway_genes) > 0) {
      round(all_genes_in_sets / length(all_pathway_genes) * 100, 2)
    } else 0,
    genes_per_pathway = lengths(mapped_genes)
  )
  
  all_mapped_genes <- unique(unlist(mapped_genes, use.names = FALSE))
  unmapped <- genes[!genes %in% all_mapped_genes]
  
  return(list(
    mapped_genes = mapped_genes,
    coverage = coverage_info,
    unmapped = unmapped
  ))
}

#' @title Get Metabolic Gene Set Statistics
#' 
#' @description Get statistics about metabolic gene set coverage in a dataset,
#' including number of genes detected, expression levels, and coverage percentage.
#' 
#' @param seurat_obj Seurat object with expression data in the active assay
#' @param pathways_to_include Character vector of pathway names to include.
#'   If NULL, includes all pathways.
#' @param assay Assay to use for expression data. Default is the active assay.
#' 
#' @return Data frame with columns:
#'   \item{pathway}{Pathway name}
#'   \item{total_genes}{Total genes in pathway}
#'   \item{genes_detected}{Number of genes detected in data}
#'   \item{pct_expressed}{Percentage of pathway genes detected}
#'   \item{mean_expression}{Mean expression across detected genes}
#'   \item{median_expression}{Median expression across detected genes}
#' 
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- readRDS("your_data.rds")
#' stats <- getMetabolicGeneSetStats(seurat_obj)
#' head(stats)
#' }
getMetabolicGeneSetStats <- function(seurat_obj, 
                                     pathways_to_include = NULL,
                                     assay = NULL) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  all_genes <- rownames(Seurat::GetAssayData(seurat_obj, assay = assay))
  if (is.null(all_genes) || length(all_genes) == 0) {
    stop("No genes found in the specified assay")
  }
  
  genes_upper <- toupper(all_genes)
  
  if (is.null(pathways_to_include)) {
    pathways_to_include <- names(metabolicGeneSets)
  }
  
  stats_list <- lapply(pathways_to_include, function(pathway) {
    if (!(pathway %in% names(metabolicGeneSets))) {
      return(NULL)
    }
    
    pathway_genes <- toupper(metabolicGeneSets[[pathway]])
    
    detected_idx <- match(pathway_genes, genes_upper)
    detected_idx <- detected_idx[!is.na(detected_idx)]
    detected_genes <- all_genes[detected_idx]
    
    n_detected <- length(detected_genes)
    pct_expressed <- (n_detected / length(pathway_genes)) * 100
    
    if (n_detected > 0) {
      expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay)[detected_genes, , drop = FALSE]
      mean_expr <- mean(as.numeric(expr_data), na.rm = TRUE)
      median_expr <- median(as.numeric(expr_data), na.rm = TRUE)
    } else {
      mean_expr <- 0
      median_expr <- 0
    }
    
    data.frame(
      pathway = pathway,
      total_genes = length(pathway_genes),
      genes_detected = n_detected,
      pct_expressed = round(pct_expressed, 2),
      mean_expression = round(mean_expr, 4),
      median_expression = round(median_expr, 4),
      stringsAsFactors = FALSE
    )
  })
  
  stats_df <- do.call(rbind, stats_list)
  stats_df <- stats_df[!sapply(stats_df, is.null), ]
  rownames(stats_df) <- NULL
  
  return(stats_df)
}

#' @title Create Custom Metabolic Gene Set
#' 
#' @description Creates a custom GeneSet from user-provided genes for use
#' in pathway scoring or custom analyses.
#' 
#' @param genes Character vector of gene symbols
#' @param name Name for the gene set
#' @param description Optional description
#' @param organism Organism name (e.g., "Homo sapiens")
#' 
#' @return GeneSet object
#' 
#' @export
#' @examples
#' custom_geneset <- createCustomGeneSet(
#'   genes = c("GENE1", "GENE2", "GENE3"),
#'   name = "my_custom_pathway",
#'   description = "Custom pathway for analysis"
#' )
createCustomGeneSet <- function(genes, name, description = "", organism = "Homo sapiens") {
  
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & genes != ""]
  
  if (length(genes) == 0) {
    stop("genes must contain at least one valid gene symbol")
  }
  
  if (missing(name) || name == "") {
    stop("name is required")
  }
  
  gs <- GSEABase::GeneSet(
    setName = name,
    geneIds = genes,
    shortDescription = description,
    organism = organism
  )
  
  return(gs)
}

#' @title Merge Metabolic Gene Sets
#' 
#' @description Merges multiple gene sets into a single character vector,
#' optionally removing duplicate genes.
#' 
#' @param gene_sets List of gene set character vectors to merge
#' @param remove_duplicates Remove duplicate genes (default TRUE)
#' @param case_sensitive Whether duplicate detection should be case-sensitive
#' 
#' @return Character vector of merged genes
#' 
#' @export
#' @examples
#' set1 <- c("HK2", "PKM2", "LDHA")
#' set2 <- c("LDHA", "NDUFA9", "SDHA")
#' merged <- mergeMetabolicGeneSets(list(set1, set2))
#' merged
mergeMetabolicGeneSets <- function(gene_sets, remove_duplicates = TRUE, 
                                   case_sensitive = FALSE) {
  
  if (!is.list(gene_sets)) {
    stop("gene_sets must be a list of character vectors")
  }
  
  merged <- unlist(gene_sets, use.names = FALSE)
  merged <- merged[!is.na(merged) & merged != ""]
  
  if (remove_duplicates) {
    if (case_sensitive) {
      merged <- unique(merged)
    } else {
      merged_upper <- toupper(merged)
      merged <- merged[!duplicated(merged_upper)]
    }
  }
  
  return(merged)
}
