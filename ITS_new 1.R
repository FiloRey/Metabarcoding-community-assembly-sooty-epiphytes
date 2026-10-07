# ===============================
# LIBRARIES -----------------------------------------------------------------
# ===============================
library(broom)
library(tidyverse)      
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(forcats)
library(tidyr)
library(stringr)
library(grid)
library(taxize)
library(gridGraphics)
library(readxl)
library(openxlsx)
library(gt)
library(qiime2R)
library(phyloseq)       
library(ggpubr)         
library(microbiome)     
library(vegan)          
library(cluster)        
library(readr)          
library(knitr)         
library(ggpmisc)
library(kableExtra) 
library(webshot2)       
library(magrittr)       
library(DESeq2)         
library(DirichletMultinomial)
library(RColorBrewer)   
library(reshape2)       
library(gplots)         
library(gridExtra) 
library(ggfortify)      
library(ggforce)        
library(ggrepel)        
library(rstatix)        
library(FSA)            
library(dunn.test)      
library(indicspecies)   
library(pheatmap)
library(BiocParallel)
library(mia)
library(miaViz)
library(scater)
library(sechm)
library(ComplexHeatmap)
library(metagenomeSeq)
library(ggExtra)
library(bluster)
library(dendextend)
library(NMF)
library(ape)
library(ggtree)
library(factoextra)
library(NbClust)
library(ALDEx2)
library(ANCOMBC)
library(seqtime)
library(devtools)  
library(maaslin3)
library(qgraph)
library(shadowtext)
library(SPRING)
library(NetCoMi)
library(igraph)
library(circlize)
library(corpcor)
library(SpiecEasi)
library(DT)
library(picante)
library(lme4)          
library(lmerTest)      
library(emmeans)       
library(car) 
library(biomformat)
library(randomForest)
library(caret)
library(mikropml)
library(pROC)
library(gbm)
library(lefser)
library(betapart)
library(eulerr)
library(limma)
library(DOC)
library(networkD3)
library(writexl)
library(ggdist)
library(ARTool)
library(tibble)
library(permute)        # restricted permutations for blocked PERMANOVA
library(iCAMP)          # phylogenetic-bin null modelling of community assembly
set.seed(423542)


# =============================================================================
# 1) explicit Block x Treatment and PlotID diagnostics;
# 2) low-read + singleton ASV QC with exported reviewer tables;
# 3) treatment-plot-aware random effects in alpha-diversity models;
# 4) split-plot, block-restricted marginal Aitchison PERMANOVA;
# 5) effect-specific PERMDISP permutation restrictions;
# 6) iCAMP bin-size diagnostics + optional ps.bin phylogenetic-signal test;
# 7) automatic mapping of selection-dominated iCAMP bins back to taxa;
# 8) conditional 16S implementation of the same reviewer-revised pipeline.
# =============================================================================

# ===============================
# COMMON PUBLICATION THEME -------------------------------------------------
# ===============================
# Base font size and family
base_text_size <- 20
base_font_family <- "Arial"  # Nature-friendly sans-serif

# Nature-style theme
theme_nature <- theme_bw(base_size = base_text_size, base_family = base_font_family) +
  theme(
    # Plot title
    plot.title = element_text(
      size = 24,
      face = "bold",
      hjust = 0.5,
      margin = ggplot2::margin(t = 0, r = 0, b = 10, l = 0)
    ),
    # Axis titles
    axis.title = element_text(
      size = 20,
      face = "bold",
      color = "black"
    ),
    # Axis text
    axis.text = element_text(
      size = 16,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.8,
      color = "black"
    ),
    # Panel
    panel.background = element_blank(),
    panel.border = element_rect(
      size = 1,
      color = "black",
      fill = NA
    ),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Legend
    legend.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 18
    ),
    # Facet strips
    strip.background = element_rect(
      fill = NA,
      color = NA
    ),
    strip.text = element_text(
      size = 18,
      face = "bold"
    ),
    # Plot margins
    plot.margin = ggplot2::margin(t = 10, r = 10, b = 10, l = 10),
    # Remove plot background for clean look
    plot.background = element_blank()
  )

# Apply globally
theme_set(theme_nature)

# ==============================
# STEP 1: IMPORT AND CLEAN DATA 
# ==============================
#physeq
physeq_asv <- qza_to_phyloseq(
  "table_ITS.qza",
  "rooted-tree_ITS_new.qza",
  "taxonomy_ITS_new.qza",
  "metadata_ITS.txt"
)

# Basic info
ntaxa(physeq_asv)
nsamples(physeq_asv)
sample_names(physeq_asv)
rank_names(physeq_asv)
sample_variables(physeq_asv)
otu_table(physeq_asv)
tax_table(physeq_asv)
phy_tree(physeq_asv)

unique_phyla <- unique(as.vector(tax_table(physeq_asv)[, "Phylum"]))
unique_order <- unique(as.vector(tax_table(physeq_asv)[, "Order"]))
unique_genus <- unique(as.vector(tax_table(physeq_asv)[, "Genus"]))

# FILTER and FORMAT --------------------------
unwanted_phyla <- c(
  NA,
  "Anthophyta", "Chordata", "Eukaryota_phy_Incertae_sedis", "Ochrophyta",
  "Chlorophyta", "Heterolobosa_phy_Incertae_sedis", "Ciliophora",  "Apicomplexa",
  "Discosea_phy_Incertae_sedis", "Protista_phy_Incertae_sedis", "Bryophyta",
  "Viridiplantae_phy_Incertae_sedis", "Nematoda", "Arthropoda", "Cercozoa",
  "Fungi_phy_Incertae_sedis" 
)

physeq_asv_filtered <- subset_taxa(
  physeq_asv,
  !(Phylum %in% unwanted_phyla | is.na(Phylum))
)

# Quick check
ntaxa(physeq_asv_filtered)
nsamples(physeq_asv_filtered)
sample_names(physeq_asv_filtered)
rank_names(physeq_asv_filtered)
sample_variables(physeq_asv_filtered)
otu_table(physeq_asv_filtered)
tax_table(physeq_asv_filtered)
phy_tree(physeq_asv_filtered)

# =============================================================================
# STEP 1B: EXPERIMENTAL-DESIGN + SEQUENCING-DEPTH QC -- REVIEWER-REVISED
# =============================================================================
# Reviewer comments addressed here:
# - explicitly verify how block and treatment are encoded;
# - identify the true treatment-plot experimental unit;
# - quantify low-read samples rather than silently discarding them;
# - quantify global singleton ASVs rather than silently discarding them.
#
# Nothing is filtered in this QC section. The exported tables document the
# dataset actually entering the inferential analyses.
# =============================================================================

get_sample_by_taxa_matrix <- function(ps) {
  x <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) x <- t(x)
  storage.mode(x) <- "numeric"
  x
}

get_taxa_by_sample_matrix <- function(ps) {
  x <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps)) x <- t(x)
  storage.mode(x) <- "numeric"
  x
}

diagnose_experimental_units <- function(ps, prefix = "ITS") {
  
  meta_qc <- data.frame(
    phyloseq::sample_data(ps)
  ) %>%
    tibble::rownames_to_column("SampleID")
  
  required_cols <- c(
    "field",
    "treatment",
    "block",
    "sampling"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(meta_qc)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  meta_qc <- meta_qc %>%
    dplyr::mutate(
      field = factor(field),
      treatment = factor(treatment),
      block = factor(block),
      sampling = factor(
        sampling,
        levels = c("Harvest", "Storage")
      ),
      
      # A/B/C are replicate labels reused within each Treatment.
      # The independent physical experimental unit is therefore:
      # Field x Treatment x Block.
      ExperimentalUnitID = interaction(
        field,
        treatment,
        block,
        drop = TRUE
      )
    )
  
  # One row per independent physical experimental unit.
  unit_structure <- meta_qc %>%
    dplyr::distinct(
      field,
      treatment,
      block,
      ExperimentalUnitID
    ) %>%
    dplyr::arrange(
      field,
      treatment,
      block
    )
  
  # Replication of Treatment at the experimental-unit level.
  treatment_replication <- unit_structure %>%
    dplyr::count(
      field,
      treatment,
      name = "N_independent_units"
    ) %>%
    dplyr::arrange(
      field,
      treatment
    )
  
  # Number of apple samples in each Sampling stage within each unit.
  sampling_structure <- meta_qc %>%
    dplyr::count(
      field,
      treatment,
      block,
      ExperimentalUnitID,
      sampling,
      name = "N_samples"
    ) %>%
    dplyr::arrange(
      field,
      treatment,
      block,
      sampling
    )
  
  # Wide table for checking Harvest/Storage balance.
  unit_sampling_balance <- sampling_structure %>%
    dplyr::select(
      field,
      treatment,
      block,
      ExperimentalUnitID,
      sampling,
      N_samples
    ) %>%
    tidyr::pivot_wider(
      names_from = sampling,
      values_from = N_samples,
      values_fill = 0
    )
  
  cat(
    "\n================ EXPERIMENTAL UNITS:",
    prefix,
    "================\n"
  )
  
  print(unit_structure)
  print(treatment_replication)
  
  cat(
    "\n================ SAMPLING BALANCE:",
    prefix,
    "================\n"
  )
  
  print(unit_sampling_balance)
  
  # Each Field x Treatment should have replicated independent units.
  if (any(treatment_replication$N_independent_units < 2L)) {
    warning(
      "At least one Field x Treatment has fewer than two ",
      "independent experimental units."
    )
  }
  
  # Each experimental unit should contain both sampling stages.
  if (
    all(
      c("Harvest", "Storage") %in%
      colnames(unit_sampling_balance)
    )
  ) {
    if (
      any(unit_sampling_balance$Harvest == 0) ||
      any(unit_sampling_balance$Storage == 0)
    ) {
      warning(
        "At least one experimental unit is missing ",
        "Harvest or Storage samples."
      )
    }
  } else {
    warning(
      "Harvest and/or Storage columns were not found in the ",
      "experimental-unit sampling-balance table."
    )
  }
  
  utils::write.csv(
    unit_structure,
    paste0(prefix, "_experimental_unit_structure.csv"),
    row.names = FALSE
  )
  
  utils::write.csv(
    treatment_replication,
    paste0(prefix, "_treatment_replication.csv"),
    row.names = FALSE
  )
  
  utils::write.csv(
    sampling_structure,
    paste0(prefix, "_experimental_unit_sampling_structure.csv"),
    row.names = FALSE
  )
  
  utils::write.csv(
    unit_sampling_balance,
    paste0(prefix, "_experimental_unit_sampling_balance.csv"),
    row.names = FALSE
  )
  
  list(
    metadata = meta_qc,
    unit_structure = unit_structure,
    treatment_replication = treatment_replication,
    sampling_structure = sampling_structure,
    unit_sampling_balance = unit_sampling_balance
  )
}

run_sequence_depth_qc <- function(
    ps,
    prefix = "ITS",
    low_read_thresholds = c(1000, 5000, 10000, 15000)
) {
  counts_qc <- get_sample_by_taxa_matrix(ps)
  lib_size_qc <- rowSums(counts_qc)
  taxa_totals_qc <- colSums(counts_qc)
  
  if (any(lib_size_qc <= 0)) {
    warning("Zero-read samples are present after taxonomic filtering.")
  }
  
  global_singletons <- taxa_totals_qc == 1
  n_singletons <- sum(global_singletons)
  
  singleton_reads_per_sample <- if (n_singletons > 0) {
    rowSums(counts_qc[, global_singletons, drop = FALSE])
  } else {
    rep(0, nrow(counts_qc))
  }
  
  meta_qc <- data.frame(phyloseq::sample_data(ps)) %>%
    tibble::rownames_to_column("SampleID")
  
  sample_qc <- tibble::tibble(
    SampleID = rownames(counts_qc),
    LibrarySize = lib_size_qc,
    SingletonReads = singleton_reads_per_sample,
    SingletonFraction = ifelse(
      lib_size_qc > 0,
      singleton_reads_per_sample / lib_size_qc,
      NA_real_
    )
  ) %>%
    dplyr::left_join(meta_qc, by = "SampleID")
  
  for (cutoff in low_read_thresholds) {
    sample_qc[[paste0("Below_", cutoff, "_reads")]] <-
      sample_qc$LibrarySize < cutoff
  }
  
  library_summary <- tibble::tibble(
    Metric = c(
      "Minimum", "1%", "5%", "25%", "Median",
      "Mean", "75%", "95%", "99%", "Maximum"
    ),
    Reads = c(
      min(lib_size_qc),
      unname(quantile(lib_size_qc, 0.01)),
      unname(quantile(lib_size_qc, 0.05)),
      unname(quantile(lib_size_qc, 0.25)),
      median(lib_size_qc),
      mean(lib_size_qc),
      unname(quantile(lib_size_qc, 0.75)),
      unname(quantile(lib_size_qc, 0.95)),
      unname(quantile(lib_size_qc, 0.99)),
      max(lib_size_qc)
    )
  )
  
  low_read_summary <- tibble::tibble(
    Threshold = low_read_thresholds,
    N_below = vapply(
      low_read_thresholds,
      function(z) sum(lib_size_qc < z),
      numeric(1)
    ),
    Percent_below = vapply(
      low_read_thresholds,
      function(z) 100 * mean(lib_size_qc < z),
      numeric(1)
    )
  )
  
  singleton_summary <- tibble::tibble(
    Total_ASVs = length(taxa_totals_qc),
    Global_singleton_ASVs = n_singletons,
    Singleton_ASV_percent = 100 * n_singletons / length(taxa_totals_qc),
    Reads_in_singleton_ASVs = sum(taxa_totals_qc[global_singletons]),
    Percent_reads_in_singletons =
      100 * sum(taxa_totals_qc[global_singletons]) / sum(taxa_totals_qc)
  )
  
  cat("\n================ LIBRARY-SIZE QC:", prefix, "================\n")
  print(library_summary)
  print(low_read_summary)
  print(singleton_summary)
  
  write.csv(
    sample_qc,
    paste0(prefix, "_sample_library_size_QC.csv"),
    row.names = FALSE
  )
  write.csv(
    library_summary,
    paste0(prefix, "_library_size_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    low_read_summary,
    paste0(prefix, "_low_read_summary.csv"),
    row.names = FALSE
  )
  write.csv(
    singleton_summary,
    paste0(prefix, "_singleton_summary.csv"),
    row.names = FALSE
  )
  
  list(
    sample_QC = sample_qc,
    library_summary = library_summary,
    low_read_summary = low_read_summary,
    singleton_summary = singleton_summary
  )
}

ITS_design <- diagnose_experimental_units(physeq_asv_filtered, prefix = "ITS")
ITS_QC <- run_sequence_depth_qc(physeq_asv_filtered, prefix = "ITS")

TSE <- convertFromPhyloseq(physeq_asv_filtered)

# Ready to use
tse__fresh <- TSE

# Quick check
tse__fresh

# ===============================
# STEP 2: SUMMARY & QUALITY CHECK ------------------------------------------
# ===============================
summary(tse__fresh, assay.type = "counts")
tse_ <- transformAssay(tse__fresh, MARGIN = "cols", method = "relabundance")

df <- summarizeDominance(
  tse_, group = "sampling", rank = "Genus", assay.type = "relabundance", top = 20
)
uniq <- getUnique(tse_, rank = "Phylum")
prev <- getPrevalent(tse_, rank = "Genus", assay.type = "relabundance",
                     prevalence = 0.75, detection = 0, sort = TRUE)
rare <- getRare(tse__fresh, rank = "Genus", prevalence = 0.75, detection = 0)

# Quality checks
tse_quality <- addPerCellQC(tse__fresh)
p1 <- plotColData(tse_quality, x = "field", y = "total", colour_by = "treatment")
p2 <- plotHistogram(tse_quality, col.var = "total")
p1 + p2

# ===============================
# STEP 3: ABUNDANCE & PREVALENCE -------------------------------------------
# ===============================
tse_ <- transformAssay(tse__fresh, MARGIN = "cols", method = "relabundance")
altExp(tse_) <- agglomerateByRank(tse_, "Genus")

plotAbundanceDensity(
  altExp(tse_), layout = "point",
  assay.type = "relabundance",
  shape_by = 'treatment',
  colour.by = "sampling",
  n = 40, point_size = 1.9, point.shape = 19,
  point.alpha = 0.3
) + scale_x_log10(label = scales::percent)

plotAbundanceDensity(
  altExp(tse_), layout = "density",
  assay.type = "relabundance",
  colour.by = "sampling",
  add_legend = FALSE,
  n = 25, point.alpha = 0.1
) + scale_x_log10()

top <- getTop(altExp(tse_), top = 10L, method = "mean")
plotExpression(
  altExp(tse_), features = top,
  x = "sampling", assay.type = "relabundance",
  point_alpha = 0.01
) + scale_y_log10()

altExp(tse_) <- addPrevalence(altExp(tse_), detection = 0.1/100, as.relative = TRUE)
plotHistogram(altExp(tse_), row.var = "prevalence")
p1 <- plotPrevalentAbundance(altExp(tse_), as.relative = TRUE)
p2 <- plotRowPrevalence(altExp(tse_), as.relative = TRUE) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
p1 + p2
plotPrevalence(altExp(tse_), as.relative = TRUE)

altExp(tse_) <- addPrevalentAbundance(altExp(tse_), prevalence = 50/100, detection = 0.1/100)
plotHistogram(altExp(tse_), col.var = "prevalent_abundance")

# ==============================


# ============================================================================
# STEP 4: ALPHA DIVERSITY -- INTEGRATED REVIEWER-REVISED ANALYSIS
# ============================================================================
#
# Purpose
# -------
# This script integrates TWO complementary approaches to unequal sequencing
# depth without fixed-depth rarefaction:
#
# A) PRIMARY ANALYSIS
#    Coverage-standardized taxonomic and phylogenetic alpha diversity using
#    iNEXT.3D, followed by field-specific mixed-effects models.
#
# B) SENSITIVITY ANALYSIS
#    Alpha-diversity metrics calculated directly from the complete,
#    non-rarefied count table, followed by field-specific mixed-effects models
#    with standardized log10 library size as a covariate.
#
# IMPORTANT: SINGLETONS ARE RETAINED.
# ----------------------------------
# No ASV is removed merely because it is a global or sample singleton.
# The script reports singleton abundance for QC, but NEVER filters singletons.
# This is important because rare ASVs contribute information about sampling
# completeness, especially for richness / q = 0.
#
# Experimental design
# -------------------
# Models are fitted separately within each orchard/field because treatment
# identities differ between orchards.
#
# Fixed effects within each field:
#     Sampling * Treatment
#
# Random effect:
#     (1 | ExperimentalUnitID)
#
# where:
#     ExperimentalUnitID = Field x Treatment x Block
#
# Experimental design:
# - A/B/C are replicate labels reused within each Treatment;
# - each Field x Treatment x Block combination is one independent physical
#   experimental unit;
# - each treatment therefore has three independent replicate units per field;
# - Harvest and Storage apples are sampled within the same experimental unit.
#
# Treatment is a BETWEEN-unit fixed effect and Sampling is a WITHIN-unit fixed
# effect. The raw Block label alone is not used as a random effect because,
# for example, Control-A and Geoxe-A are different physical units.
#
# Primary coverage-standardized metrics do NOT include library size again,
# because sequencing effort has already been standardized by sample coverage.
#
# Observed-table sensitivity metrics and phylogenetic-structure metrics DO
# include LibrarySize_z.
#
# Outputs include:
# - singleton and library-size QC;
# - iNEXT coverage QC;
# - primary and sensitivity alpha-diversity tables;
# - mixed-model Type-III tests;
# - residual/convergence/singularity diagnostics;
# - Harvest vs Storage EMM contrasts for all metrics;
# - Harvest vs Storage contrasts within each treatment;
# - treatment pairwise contrasts within each Sampling stage;
# - a robustness/concordance table comparing iNEXT vs library-size-adjusted
#   counterparts;
# - main alpha-diversity figure using original metric scales and EMM deltas.
# ============================================================================


# ============================================================================
# 0. REQUIRED PACKAGES
# ============================================================================

required_packages <- c(
  "phyloseq",
  "iNEXT.3D",
  "ape",
  "picante",
  "vegan",
  "dplyr",
  "tidyr",
  "tibble",
  "purrr",
  "ggplot2",
  "ggh4x",
  "lme4",
  "lmerTest",
  "emmeans"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Install the following packages before running STEP 4:\n",
    paste(missing_packages, collapse = ", ")
  )
}

# Use the manuscript theme if already defined in the master script.
theme_alpha <- if (exists("theme_nature")) {
  theme_nature
} else {
  ggplot2::theme_bw(base_size = 14)
}

# Reproducibility for any plotting jitter.
set.seed(423542)


# ============================================================================
# 1. FULL NON-RAREFIED ASV TABLE -- SINGLETONS RETAINED
# ============================================================================

alpha_counts <- as(
  phyloseq::otu_table(physeq_asv_filtered),
  "matrix"
)

if (phyloseq::taxa_are_rows(physeq_asv_filtered)) {
  alpha_counts <- t(alpha_counts)
}

storage.mode(alpha_counts) <- "numeric"

# alpha_counts = samples x ASVs
library_size <- rowSums(alpha_counts)
taxa_totals  <- colSums(alpha_counts)

if (any(library_size <= 0)) {
  stop("Samples with zero library size are present after taxonomic filtering.")
}

# Explicit audit: no taxa have been removed here.
if (ncol(alpha_counts) != phyloseq::ntaxa(physeq_asv_filtered)) {
  stop("Unexpected ASV loss while constructing alpha_counts.")
}

# ---------------------------------------------------------------------------
# Singleton QC ONLY -- DO NOT FILTER
# ---------------------------------------------------------------------------

global_singleton <- taxa_totals == 1
n_global_singletons <- sum(global_singleton)

# Number of within-sample count==1 observations (sample singletons).
n_sample_singleton_occurrences <- sum(alpha_counts == 1)

singleton_reads_per_sample <- if (n_global_singletons > 0L) {
  rowSums(alpha_counts[, global_singleton, drop = FALSE])
} else {
  rep(0, nrow(alpha_counts))
}

singleton_summary <- tibble::tibble(
  Total_ASVs = ncol(alpha_counts),
  Global_singleton_ASVs = n_global_singletons,
  Global_singleton_percent = 100 * n_global_singletons / ncol(alpha_counts),
  Reads_in_global_singletons = sum(taxa_totals[global_singleton]),
  Percent_reads_in_global_singletons =
    100 * sum(taxa_totals[global_singleton]) / sum(taxa_totals),
  Sample_singleton_occurrences = n_sample_singleton_occurrences,
  Singleton_filter_applied = FALSE
)

print(singleton_summary)

utils::write.csv(
  singleton_summary,
  "alpha_singleton_QC_RETAINED.csv",
  row.names = FALSE
)

message(
  "Singleton policy: RETAINED. Global singleton ASVs = ",
  n_global_singletons,
  ". No singleton filtering is applied in STEP 4."
)

# iNEXT.3D expects ASVs x samples.
alpha_abun <- t(alpha_counts)


# ============================================================================
# 2. METADATA + EXPERIMENTAL UNITS
# ============================================================================

alpha_meta <- data.frame(
  phyloseq::sample_data(physeq_asv_filtered)
) %>%
  tibble::rownames_to_column("sam_name")


required_metadata <- c(
  "field",
  "treatment",
  "block",
  "sampling"
)

missing_metadata <- setdiff(
  required_metadata,
  colnames(alpha_meta)
)

if (length(missing_metadata) > 0L) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata, collapse = ", ")
  )
}


alpha_meta <- alpha_meta %>%
  dplyr::mutate(
    
    Sampling = factor(
      sampling,
      levels = c("Harvest", "Storage")
    ),
    
    Treatment = factor(treatment),
    
    Field = factor(field),
    
    Block = factor(block),
    
    # A/B/C are replicate labels reused within Treatment.
    # The physical independent experimental unit is therefore
    # Field x Treatment x Block.
    ExperimentalUnitID = interaction(
      Field,
      Treatment,
      Block,
      drop = TRUE
    )
  )


# ============================================================================
# CHECK EXPERIMENTAL DESIGN
# ============================================================================

experimental_unit_structure <- alpha_meta %>%
  dplyr::distinct(
    Field,
    Treatment,
    Block,
    ExperimentalUnitID
  ) %>%
  dplyr::arrange(
    Field,
    Treatment,
    Block
  )

print(experimental_unit_structure)


# Number of independent replicate units per treatment
treatment_replication <- experimental_unit_structure %>%
  dplyr::count(
    Field,
    Treatment,
    name = "N_independent_units"
  )

print(treatment_replication)


if (any(treatment_replication$N_independent_units < 2L)) {
  
  stop(
    "At least one treatment has fewer than two ",
    "independent experimental units."
  )
}


# Check that each experimental unit contains both sampling stages
unit_sampling <- alpha_meta %>%
  dplyr::distinct(
    Field,
    ExperimentalUnitID,
    Sampling
  ) %>%
  dplyr::count(
    Field,
    ExperimentalUnitID,
    name = "N_sampling_stages"
  )

print(unit_sampling)


if (any(unit_sampling$N_sampling_stages != 2L)) {
  
  warning(
    "At least one experimental unit does not contain ",
    "both Harvest and Storage."
  )
}
# ============================================================================
# 3. LIBRARY-SIZE + SINGLETON QC BY SAMPLE
# ============================================================================

library_qc <- tibble::tibble(
  sam_name = rownames(alpha_counts),
  LibrarySize = as.numeric(library_size),
  GlobalSingletonReads = as.numeric(singleton_reads_per_sample),
  GlobalSingletonFraction = singleton_reads_per_sample / library_size,
  SampleSingletonASVs = rowSums(alpha_counts == 1)
) %>%
  dplyr::left_join(
    alpha_meta,
    by = "sam_name"
  )

utils::write.csv(
  library_qc,
  "alpha_library_and_singleton_QC.csv",
  row.names = FALSE
)

p_library_size <- ggplot2::ggplot(
  library_qc,
  ggplot2::aes(
    x = LibrarySize,
    fill = Sampling
  )
) +
  ggplot2::geom_histogram(
    bins = 30,
    alpha = 0.65,
    position = "identity"
  ) +
  ggplot2::scale_x_log10() +
  ggplot2::labs(
    x = "Library size (reads, log10 scale)",
    y = "Number of samples",
    fill = "Sampling"
  ) +
  theme_alpha

print(p_library_size)
saveRDS(p_library_size, "alpha_library_size_QC.RDS")


# ============================================================================
# 4. PREPARE PHYLOGENETIC DATA
# ============================================================================

alpha_tree <- phyloseq::phy_tree(physeq_asv_filtered)

common_taxa_alpha <- intersect(
  colnames(alpha_counts),
  alpha_tree$tip.label
)

if (length(common_taxa_alpha) < 2L) {
  stop("Too few ASVs overlap between count table and phylogenetic tree.")
}

alpha_counts_phy <- alpha_counts[, common_taxa_alpha, drop = FALSE]
alpha_abun_phy   <- t(alpha_counts_phy)
alpha_tree_phy   <- ape::keep.tip(alpha_tree, common_taxa_alpha)

phy_reads <- rowSums(alpha_counts_phy)
phy_fraction <- phy_reads / library_size

message(
  "Median fraction of reads represented in phylogenetic tree = ",
  round(stats::median(phy_fraction, na.rm = TRUE), 4)
)

phy_tree_qc <- tibble::tibble(
  sam_name = rownames(alpha_counts_phy),
  ReadsTotal = as.numeric(library_size[rownames(alpha_counts_phy)]),
  ReadsInTree = as.numeric(phy_reads),
  FractionReadsInTree = as.numeric(phy_fraction)
)

utils::write.csv(
  phy_tree_qc,
  "alpha_phylogenetic_tree_coverage_QC.csv",
  row.names = FALSE
)


# ============================================================================
# PART A -- PRIMARY COVERAGE-STANDARDIZED ANALYSIS
# ============================================================================

# ============================================================================
# 5A. SAMPLE-COVERAGE INFORMATION
# ============================================================================

coverage_info_TD <- iNEXT.3D::DataInfo3D(
  data = alpha_abun,
  diversity = "TD",
  datatype = "abundance"
)

coverage_info_PD <- iNEXT.3D::DataInfo3D(
  data = alpha_abun_phy,
  diversity = "PD",
  datatype = "abundance",
  PDtree = alpha_tree_phy
)

find_coverage_column <- function(dat, doubled = FALSE) {
  pattern <- if (doubled) "^SC\\(2 *n\\)$" else "^SC\\(n\\)$"
  
  hit <- grep(
    pattern,
    colnames(dat),
    value = TRUE
  )
  
  if (length(hit) != 1L) {
    stop(
      "Could not identify sample-coverage column. Columns are: ",
      paste(colnames(dat), collapse = ", ")
    )
  }
  
  hit
}

SC_n_TD_col  <- find_coverage_column(coverage_info_TD, doubled = FALSE)
SC_2n_TD_col <- find_coverage_column(coverage_info_TD, doubled = TRUE)
SC_n_PD_col  <- find_coverage_column(coverage_info_PD, doubled = FALSE)
SC_2n_PD_col <- find_coverage_column(coverage_info_PD, doubled = TRUE)

# "2n" allows up to approximately 2x extrapolation.
# "observed" avoids extrapolation entirely.
INEXT_COVERAGE_RULE <- "2n"

if (INEXT_COVERAGE_RULE == "2n") {
  target_TD <- min(coverage_info_TD[[SC_2n_TD_col]], na.rm = TRUE)
  target_PD <- min(coverage_info_PD[[SC_2n_PD_col]], na.rm = TRUE)
} else if (INEXT_COVERAGE_RULE == "observed") {
  target_TD <- min(coverage_info_TD[[SC_n_TD_col]], na.rm = TRUE)
  target_PD <- min(coverage_info_PD[[SC_n_PD_col]], na.rm = TRUE)
} else {
  stop("INEXT_COVERAGE_RULE must be '2n' or 'observed'.")
}

# One common target across taxonomic and phylogenetic analyses.
target_coverage <- min(target_TD, target_PD, 0.999999)
target_coverage <- floor(target_coverage * 100000) / 100000

message("Common iNEXT target coverage = ", target_coverage)

coverage_qc <- coverage_info_TD %>%
  dplyr::transmute(
    sam_name = as.character(Assemblage),
    LibrarySize = n,
    ObservedRichness = S.obs,
    CoverageObserved = .data[[SC_n_TD_col]],
    Coverage2x = .data[[SC_2n_TD_col]]
  ) %>%
  dplyr::left_join(alpha_meta, by = "sam_name")

utils::write.csv(
  coverage_qc,
  "iNEXT_sample_coverage_QC.csv",
  row.names = FALSE
)

p_coverage <- ggplot2::ggplot(
  coverage_qc,
  ggplot2::aes(
    x = Sampling,
    y = CoverageObserved,
    fill = Sampling
  )
) +
  ggplot2::geom_boxplot(
    width = 0.35,
    outlier.shape = NA,
    alpha = 0.7
  ) +
  ggplot2::geom_jitter(
    width = 0.10,
    size = 1.5,
    alpha = 0.55
  ) +
  ggplot2::geom_hline(
    yintercept = target_coverage,
    linetype = 2
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Estimated sample coverage",
    fill = "Sampling"
  ) +
  theme_alpha +
  ggplot2::theme(legend.position = "none")

print(p_coverage)
saveRDS(p_coverage, "iNEXT_sample_coverage_QC.RDS")


# ============================================================================
# 6A. COVERAGE-STANDARDIZED TAXONOMIC DIVERSITY
# ============================================================================
# q = 0 : richness
# q = 1 : Shannon effective diversity = exp(Shannon entropy)
# q = 2 : inverse-Simpson effective diversity
# q = 2 is NOT evenness.
# ============================================================================

TD_est <- iNEXT.3D::estimate3D(
  data = alpha_abun,
  diversity = "TD",
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  level = target_coverage,
  nboot = 0
)

utils::write.csv(
  TD_est,
  "iNEXT_TD_coverage_standardized_full.csv",
  row.names = FALSE
)

print(
  TD_est %>%
    dplyr::count(Order.q, Method)
)

TD_wide <- TD_est %>%
  dplyr::transmute(
    sam_name = as.character(Assemblage),
    Metric = paste0("TD_q", as.integer(Order.q)),
    Diversity = qTD
  ) %>%
  tidyr::pivot_wider(
    names_from = Metric,
    values_from = Diversity
  )


# ============================================================================
# 7A. COVERAGE-STANDARDIZED PHYLOGENETIC DIVERSITY
# ============================================================================

PD_est <- iNEXT.3D::estimate3D(
  data = alpha_abun_phy,
  diversity = "PD",
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  level = target_coverage,
  nboot = 0,
  PDtree = alpha_tree_phy,
  PDtype = "meanPD"
)

# q=0 total branch-length form, analogous to Faith's PD.
FaithPD_est <- iNEXT.3D::estimate3D(
  data = alpha_abun_phy,
  diversity = "PD",
  q = 0,
  datatype = "abundance",
  base = "coverage",
  level = target_coverage,
  nboot = 0,
  PDtree = alpha_tree_phy,
  PDtype = "PD"
)

utils::write.csv(
  PD_est,
  "iNEXT_PD_meanPD_coverage_standardized_full.csv",
  row.names = FALSE
)

utils::write.csv(
  FaithPD_est,
  "iNEXT_FaithPD_coverage_standardized_full.csv",
  row.names = FALSE
)

PD_wide <- PD_est %>%
  dplyr::transmute(
    sam_name = as.character(Assemblage),
    Metric = paste0("PD_q", as.integer(Order.q)),
    Diversity = qPD
  ) %>%
  tidyr::pivot_wider(
    names_from = Metric,
    values_from = Diversity
  )

FaithPD_cov <- FaithPD_est %>%
  dplyr::transmute(
    sam_name = as.character(Assemblage),
    FaithPD_cov = qPD
  )

# ============================================================================
# 8A. PRIMARY COVERAGE-STANDARDIZED DATASET
# ============================================================================

alpha_primary <- TD_wide %>%
  dplyr::left_join(PD_wide, by = "sam_name") %>%
  dplyr::left_join(FaithPD_cov, by = "sam_name") %>%
  dplyr::left_join(alpha_meta, by = "sam_name") %>%
  dplyr::mutate(
    HillEvenness = dplyr::if_else(
      TD_q0 > 0,
      TD_q1 / TD_q0,
      NA_real_
    )
  )

utils::write.csv(
  alpha_primary,
  "alpha_PRIMARY_iNEXT_coverage_standardized.csv",
  row.names = FALSE
)

primary_metrics <- c(
  "TD_q0",
  "TD_q1",
  "TD_q2",
  "HillEvenness",
  "FaithPD_cov",
  "PD_q1",
  "PD_q2"
)

primary_labels <- c(
  TD_q0 = "Richness (q = 0)",
  TD_q1 = "Shannon effective diversity (q = 1)",
  TD_q2 = "Inverse Simpson effective diversity (q = 2)",
  HillEvenness = "Hill evenness (D1/D0)",
  FaithPD_cov = "Faith's phylogenetic diversity",
  PD_q1 = "Phylogenetic diversity (q = 1)",
  PD_q2 = "Phylogenetic diversity (q = 2)"
)


# ============================================================================
# PART B -- NON-RAREFIED + LIBRARY-SIZE-COVARIATE SENSITIVITY ANALYSIS
# ============================================================================

# ============================================================================
# 5B. OBSERVED TAXONOMIC ALPHA DIVERSITY -- FULL COUNTS
# ============================================================================
# SINGLETONS REMAIN PRESENT.
# ============================================================================

Observed_raw <- vegan::specnumber(alpha_counts)
Shannon_raw <- vegan::diversity(alpha_counts, index = "shannon")
ShannonEffective_raw <- exp(Shannon_raw)
InvSimpson_raw <- vegan::diversity(alpha_counts, index = "invsimpson")
Pielou_raw <- ifelse(
  Observed_raw > 1,
  Shannon_raw / log(Observed_raw),
  NA_real_
)

alpha_tax_raw <- tibble::tibble(
  sam_name = rownames(alpha_counts),
  LibrarySize = as.numeric(library_size),
  Observed_raw = Observed_raw,
  Shannon_raw = Shannon_raw,
  ShannonEffective_raw = ShannonEffective_raw,
  InvSimpson_raw = InvSimpson_raw,
  Pielou_raw = Pielou_raw
)


# ============================================================================
# 6B. OBSERVED PHYLOGENETIC ALPHA DIVERSITY / STRUCTURE
# ============================================================================

faith_pd_raw <- picante::pd(
  alpha_counts_phy,
  alpha_tree_phy,
  include.root = TRUE
)

psv_out <- picante::psv(
  alpha_counts_phy,
  alpha_tree_phy,
  compute.var = FALSE
)

psr_out <- picante::psr(
  alpha_counts_phy,
  alpha_tree_phy,
  compute.var = FALSE
)

pse_out <- picante::pse(
  alpha_counts_phy,
  alpha_tree_phy
)

alpha_phy_raw <- tibble::tibble(
  sam_name = rownames(alpha_counts_phy),
  FaithPD_raw = faith_pd_raw$PD,
  PSV = psv_out$PSV,
  PSR = psr_out$PSR,
  PSE = pse_out$PSE,
  # Explicitly derived convenience index; interpret as relatedness/redundancy,
  # not as an independent standard phylogenetic-diversity estimator.
  PhyloRedundancy = 1 - psv_out$PSV
)


# ============================================================================
# 7B. SENSITIVITY DATASET + LIBRARY-SIZE COVARIATE
# ============================================================================

alpha_sensitivity <- alpha_tax_raw %>%
  dplyr::left_join(alpha_phy_raw, by = "sam_name") %>%
  dplyr::left_join(alpha_meta, by = "sam_name") %>%
  dplyr::mutate(
    LibrarySize_log10 = log10(LibrarySize),
    LibrarySize_z = as.numeric(scale(LibrarySize_log10))
  )

utils::write.csv(
  alpha_sensitivity,
  "alpha_SENSITIVITY_nonrarefied_librarysize_adjusted.csv",
  row.names = FALSE
)

sensitivity_metrics <- c(
  "Observed_raw",
  "Shannon_raw",
  "ShannonEffective_raw",
  "InvSimpson_raw",
  "Pielou_raw",
  "FaithPD_raw",
  "PSV",
  "PSR",
  "PSE",
  "PhyloRedundancy"
)

sensitivity_labels <- c(
  Observed_raw = "Observed richness",
  Shannon_raw = "Shannon entropy",
  ShannonEffective_raw = "Shannon effective diversity",
  InvSimpson_raw = "Inverse Simpson diversity",
  Pielou_raw = "Pielou evenness",
  FaithPD_raw = "Faith's phylogenetic diversity",
  PSV = "Phylogenetic species variability",
  PSR = "Phylogenetic species richness",
  PSE = "Phylogenetic species evenness",
  PhyloRedundancy = "Phylogenetic relatedness (1 - PSV)"
)


# ============================================================================
# 9. FIELD-SPECIFIC MIXED-MODEL FITTING FUNCTIONS
# ============================================================================
# The two orchards have different treatment identities, so models are fitted
# separately within each Field rather than using a global Field x Treatment
# factorial.
#
# The independent physical experimental unit is Field x Treatment x Block.
# The A/B/C labels are replicate labels reused within each Treatment, so the raw
# Block label alone is not a shared blocking factor across treatments.
# Multiple apples and both Sampling stages occur within each ExperimentalUnitID.
# ExperimentalUnitID is therefore included as a random intercept.
# Treatment is a fixed BETWEEN-unit factor and Sampling is a fixed WITHIN-unit
# factor.
# ============================================================================

prepare_field_data <- function(dat, field_i) {
  out <- dat %>%
    dplyr::filter(Field == field_i) %>%
    droplevels()
  
  if (nlevels(out$Sampling) != 2L) {
    stop("Field ", field_i, " does not contain both Sampling levels.")
  }
  
  if (nlevels(out$Treatment) < 2L) {
    stop("Field ", field_i, " contains <2 Treatment levels.")
  }
  
  contrasts(out$Sampling) <- stats::contr.sum(nlevels(out$Sampling))
  contrasts(out$Treatment) <- stats::contr.sum(nlevels(out$Treatment))
  
  out
}

fit_primary_model <- function(field_i, metric_i) {
  dat <- prepare_field_data(alpha_primary, field_i)
  
  f <- stats::as.formula(
    paste0(
      metric_i,
      " ~ Sampling * Treatment + ",
      "(1 | ExperimentalUnitID)"
    )
  )
  
  message("PRIMARY | ", field_i, " | ", metric_i, " | ", deparse(f))
  
  lmerTest::lmer(
    f,
    data = dat,
    REML = TRUE,
    na.action = na.omit,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
}

fit_sensitivity_model <- function(field_i, metric_i) {
  dat <- prepare_field_data(alpha_sensitivity, field_i)
  
  f <- stats::as.formula(
    paste0(
      metric_i,
      " ~ Sampling * Treatment + LibrarySize_z + ",
      "(1 | ExperimentalUnitID)"
    )
  )
  
  message("SENSITIVITY | ", field_i, " | ", metric_i, " | ", deparse(f))
  
  lmerTest::lmer(
    f,
    data = dat,
    REML = TRUE,
    na.action = na.omit,
    control = lme4::lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 200000)
    )
  )
}

fields_alpha <- levels(droplevels(alpha_meta$Field))

primary_registry <- tidyr::expand_grid(
  Field = fields_alpha,
  Metric = primary_metrics
) %>%
  dplyr::mutate(
    Analysis = "Primary_coverage_standardized",
    MetricLabel = unname(primary_labels[Metric]),
    Model = purrr::map2(Field, Metric, fit_primary_model)
  )

sensitivity_registry <- tidyr::expand_grid(
  Field = fields_alpha,
  Metric = sensitivity_metrics
) %>%
  dplyr::mutate(
    Analysis = "Sensitivity_nonrarefied_librarysize",
    MetricLabel = unname(sensitivity_labels[Metric]),
    Model = purrr::map2(Field, Metric, fit_sensitivity_model)
  )

model_registry <- dplyr::bind_rows(
  primary_registry,
  sensitivity_registry
)


# ============================================================================
# 10. MODEL DIAGNOSTICS
# ============================================================================
# Raw response values do not have to be normally distributed.
# For Gaussian LMMs, inspect residuals, QQ plots, heteroscedasticity patterns,
# convergence, extreme residuals and singular random-effect fits.
# ============================================================================

dir.create(
  "alpha_model_diagnostics_integrated",
  showWarnings = FALSE
)

extract_model_diagnostics <- function(
    model_i,
    field_i,
    metric_i,
    analysis_i
) {
  residual_i <- stats::residuals(model_i)
  fitted_i <- stats::fitted(model_i)
  
  std_residual_i <- if (stats::sd(residual_i, na.rm = TRUE) > 0) {
    as.numeric(scale(residual_i))
  } else {
    rep(0, length(residual_i))
  }
  
  if (length(residual_i) >= 3L && length(residual_i) <= 5000L) {
    shapiro_i <- stats::shapiro.test(residual_i)
    Shapiro_W <- unname(shapiro_i$statistic)
    Shapiro_p <- shapiro_i$p.value
  } else {
    Shapiro_W <- NA_real_
    Shapiro_p <- NA_real_
  }
  
  conv_message <- model_i@optinfo$conv$lme4$messages
  if (is.null(conv_message)) {
    conv_message <- NA_character_
  } else {
    conv_message <- paste(conv_message, collapse = " | ")
  }
  
  tibble::tibble(
    Analysis = analysis_i,
    Field = field_i,
    Metric = metric_i,
    N = stats::nobs(model_i),
    AIC = stats::AIC(model_i),
    BIC = stats::BIC(model_i),
    Singular = lme4::isSingular(model_i, tol = 1e-4),
    ConvergenceMessage = conv_message,
    Shapiro_W = Shapiro_W,
    Shapiro_p = Shapiro_p,
    MaxAbsStdResidual = max(abs(std_residual_i), na.rm = TRUE),
    N_StdResidual_gt3 = sum(abs(std_residual_i) > 3, na.rm = TRUE),
    Cor_AbsResidual_Fitted = suppressWarnings(
      stats::cor(abs(residual_i), fitted_i, use = "complete.obs")
    )
  )
}

model_diagnostics <- purrr::pmap_dfr(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$Analysis
  ),
  extract_model_diagnostics
)

model_diagnostics <- model_diagnostics %>%
  dplyr::mutate(
    Flag_Singular = Singular,
    Flag_Convergence = !is.na(ConvergenceMessage),
    Flag_ExtremeResidual = N_StdResidual_gt3 > 0,
    Flag_HeteroscedasticTrend =
      !is.na(Cor_AbsResidual_Fitted) & abs(Cor_AbsResidual_Fitted) > 0.30,
    # Shapiro is reported but deliberately NOT included in the overall flag.
    # A significant Shapiro test alone is not a reason to reject an LMM.
    Flag_Shapiro_only = !is.na(Shapiro_p) & Shapiro_p < 0.05,
    AnyMajorDiagnosticFlag =
      Flag_Singular |
      Flag_Convergence |
      Flag_ExtremeResidual |
      Flag_HeteroscedasticTrend
  )

print(model_diagnostics)

utils::write.csv(
  model_diagnostics,
  "alpha_integrated_model_diagnostics_summary.csv",
  row.names = FALSE
)

save_diagnostic_plot <- function(
    model_i,
    field_i,
    metric_i,
    analysis_i
) {
  safe <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
  
  file_i <- file.path(
    "alpha_model_diagnostics_integrated",
    paste0(
      safe(analysis_i), "__",
      safe(field_i), "__",
      safe(metric_i), ".pdf"
    )
  )
  
  residual_i <- stats::residuals(model_i)
  fitted_i <- stats::fitted(model_i)
  
  grDevices::pdf(file_i, width = 9, height = 4.5)
  graphics::par(mfrow = c(1, 2))
  
  graphics::plot(
    fitted_i,
    residual_i,
    pch = 16,
    cex = 0.7,
    xlab = "Fitted values",
    ylab = "Residuals",
    main = paste(analysis_i, field_i, metric_i, "Residuals vs fitted", sep = "\n")
  )
  graphics::abline(h = 0, lty = 2)
  
  stats::qqnorm(
    residual_i,
    pch = 16,
    cex = 0.7,
    main = paste(analysis_i, field_i, metric_i, "Normal Q-Q", sep = "\n")
  )
  stats::qqline(residual_i)
  
  grDevices::dev.off()
}

purrr::pwalk(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$Analysis
  ),
  save_diagnostic_plot
)


# ============================================================================
# 11. RANDOM-EFFECT VARIANCES
# ============================================================================

random_effect_variances <- purrr::pmap_dfr(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$Analysis
  ),
  function(model_i, field_i, metric_i, analysis_i) {
    as.data.frame(lme4::VarCorr(model_i)) %>%
      dplyr::mutate(
        Analysis = analysis_i,
        Field = field_i,
        Metric = metric_i,
        .before = 1
      )
  }
)

utils::write.csv(
  random_effect_variances,
  "alpha_integrated_random_effect_variances.csv",
  row.names = FALSE
)


# ============================================================================
# 12. TYPE-III OMNIBUS TESTS
# ============================================================================

alpha_type3 <- purrr::pmap_dfr(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$MetricLabel,
    model_registry$Analysis
  ),
  function(model_i, field_i, metric_i, label_i, analysis_i) {
    as.data.frame(
      stats::anova(
        model_i,
        type = 3,
        ddf = "Satterthwaite"
      )
    ) %>%
      tibble::rownames_to_column("Effect") %>%
      dplyr::mutate(
        Analysis = analysis_i,
        Field = field_i,
        Metric = metric_i,
        MetricLabel = label_i,
        .before = 1
      )
  }
)

utils::write.csv(
  alpha_type3,
  "TABLE_alpha_integrated_TypeIII_by_field.csv",
  row.names = FALSE
)

# Explicit library-size effects from sensitivity models.
library_size_effects <- alpha_type3 %>%
  dplyr::filter(
    Analysis == "Sensitivity_nonrarefied_librarysize",
    Effect == "LibrarySize_z"
  )

utils::write.csv(
  library_size_effects,
  "TABLE_alpha_sensitivity_library_size_effects.csv",
  row.names = FALSE
)


# ============================================================================
# 13. POST-HOC HELPER FUNCTIONS
# ============================================================================

extract_stage_overall <- function(
    model_i,
    field_i,
    metric_i,
    label_i,
    analysis_i
) {
  emm_i <- emmeans::emmeans(
    model_i,
    ~ Sampling,
    weights = "equal"
  )
  
  emmeans::contrast(
    emm_i,
    method = list("Storage - Harvest" = c(-1, 1)),
    adjust = "none"
  ) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Analysis = analysis_i,
      Field = field_i,
      Metric = metric_i,
      MetricLabel = label_i,
      .before = 1
    )
}

extract_stage_within_treatment <- function(
    model_i,
    field_i,
    metric_i,
    label_i,
    analysis_i
) {
  emm_i <- emmeans::emmeans(
    model_i,
    ~ Sampling | Treatment
  )
  
  emmeans::contrast(
    emm_i,
    method = list("Storage - Harvest" = c(-1, 1)),
    adjust = "none"
  ) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Analysis = analysis_i,
      Field = field_i,
      Metric = metric_i,
      MetricLabel = label_i,
      .before = 1
    )
}

extract_treatment_pairwise <- function(
    model_i,
    field_i,
    metric_i,
    label_i,
    analysis_i
) {
  emm_i <- emmeans::emmeans(
    model_i,
    ~ Treatment | Sampling
  )
  
  emmeans::contrast(
    emm_i,
    method = "pairwise",
    adjust = "tukey"
  ) %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    dplyr::mutate(
      Analysis = analysis_i,
      Field = field_i,
      Metric = metric_i,
      MetricLabel = label_i,
      .before = 1
    )
}


# ============================================================================
# 14. HARVEST vs STORAGE -- ALL METRICS WITHIN EACH FIELD
# ============================================================================
# estimate = original-scale EMM difference: Storage - Harvest.
# ============================================================================

stage_overall_all <- purrr::pmap_dfr(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$MetricLabel,
    model_registry$Analysis
  ),
  extract_stage_overall
)

# Holm within each Analysis x Field family across reported alpha metrics.
stage_overall_all <- stage_overall_all %>%
  dplyr::group_by(Analysis, Field) %>%
  dplyr::mutate(
    p.value.Holm = stats::p.adjust(p.value, method = "holm")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Significance = dplyr::case_when(
      p.value.Holm < 0.001 ~ "***",
      p.value.Holm < 0.01  ~ "**",
      p.value.Holm < 0.05  ~ "*",
      TRUE                 ~ "ns"
    )
  )

utils::write.csv(
  stage_overall_all,
  "TABLE_alpha_Harvest_vs_Storage_ALL_integrated.csv",
  row.names = FALSE
)

utils::write.csv(
  stage_overall_all %>%
    dplyr::filter(Analysis == "Primary_coverage_standardized"),
  "TABLE_alpha_Harvest_vs_Storage_PRIMARY_iNEXT.csv",
  row.names = FALSE
)

utils::write.csv(
  stage_overall_all %>%
    dplyr::filter(Analysis == "Sensitivity_nonrarefied_librarysize"),
  "TABLE_alpha_Harvest_vs_Storage_SENSITIVITY_librarysize.csv",
  row.names = FALSE
)


# ============================================================================
# 15. HARVEST vs STORAGE WITHIN EACH TREATMENT
# ============================================================================

stage_within_treatment_all <- purrr::pmap_dfr(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$MetricLabel,
    model_registry$Analysis
  ),
  extract_stage_within_treatment
)

# Holm across treatment-specific Harvest-vs-Storage tests within each
# Analysis x Field x Metric family.
stage_within_treatment_all <- stage_within_treatment_all %>%
  dplyr::group_by(Analysis, Field, Metric) %>%
  dplyr::mutate(
    p.value.Holm = stats::p.adjust(p.value, method = "holm")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Significance = dplyr::case_when(
      p.value.Holm < 0.001 ~ "***",
      p.value.Holm < 0.01  ~ "**",
      p.value.Holm < 0.05  ~ "*",
      TRUE                 ~ "ns"
    )
  )

utils::write.csv(
  stage_within_treatment_all,
  "TABLE_alpha_Harvest_vs_Storage_WITHIN_TREATMENT_integrated.csv",
  row.names = FALSE
)


# ============================================================================
# 16. TREATMENT COMPARISONS WITHIN EACH FIELD x SAMPLING STAGE
# ============================================================================

treatment_pairwise_all <- purrr::pmap_dfr(
  list(
    model_registry$Model,
    model_registry$Field,
    model_registry$Metric,
    model_registry$MetricLabel,
    model_registry$Analysis
  ),
  extract_treatment_pairwise
) %>%
  dplyr::mutate(
    Significance = dplyr::case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ "ns"
    )
  )

utils::write.csv(
  treatment_pairwise_all,
  "TABLE_alpha_TREATMENT_PAIRWISE_within_Field_Sampling_integrated.csv",
  row.names = FALSE
)


# ============================================================================
# 17. ROBUSTNESS / CONCORDANCE TABLE:
#     iNEXT COVERAGE STANDARDIZATION vs LIBRARY-SIZE COVARIATE
# ============================================================================
# Conceptually matched pairs:
#   TD_q0       <-> Observed_raw
#   TD_q1       <-> ShannonEffective_raw
#   TD_q2       <-> InvSimpson_raw
#   HillEvenness<-> Pielou_raw        (both evenness, not identical formulae)
#   FaithPD_cov <-> FaithPD_raw
# ============================================================================

metric_crosswalk <- tibble::tribble(
  ~PrimaryMetric, ~SensitivityMetric, ~ComparisonMeaning,
  "TD_q0",        "Observed_raw",         "Richness",
  "TD_q1",        "ShannonEffective_raw", "Shannon effective diversity",
  "TD_q2",        "InvSimpson_raw",       "Inverse Simpson effective diversity",
  "HillEvenness", "Pielou_raw",           "Evenness (related but not identical estimators)",
  "FaithPD_cov",  "FaithPD_raw",          "Faith's phylogenetic diversity"
)

primary_stage <- stage_overall_all %>%
  dplyr::filter(Analysis == "Primary_coverage_standardized") %>%
  dplyr::select(
    Field,
    PrimaryMetric = Metric,
    PrimaryEstimate = estimate,
    PrimarySE = SE,
    PrimaryP = p.value,
    PrimaryPHolm = p.value.Holm,
    PrimarySignificance = Significance
  )

sensitivity_stage <- stage_overall_all %>%
  dplyr::filter(Analysis == "Sensitivity_nonrarefied_librarysize") %>%
  dplyr::select(
    Field,
    SensitivityMetric = Metric,
    SensitivityEstimate = estimate,
    SensitivitySE = SE,
    SensitivityP = p.value,
    SensitivityPHolm = p.value.Holm,
    SensitivitySignificance = Significance
  )

robustness_table <- metric_crosswalk %>%
  dplyr::left_join(primary_stage, by = "PrimaryMetric") %>%
  dplyr::left_join(
    sensitivity_stage,
    by = c("Field", "SensitivityMetric")
  ) %>%
  dplyr::mutate(
    SameDirection =
      sign(PrimaryEstimate) == sign(SensitivityEstimate),
    PrimarySignificant = PrimaryPHolm < 0.05,
    SensitivitySignificant = SensitivityPHolm < 0.05,
    SameSignificanceConclusion =
      PrimarySignificant == SensitivitySignificant
  )

print(robustness_table)

utils::write.csv(
  robustness_table,
  "TABLE_alpha_ROBUSTNESS_iNEXT_vs_librarysize.csv",
  row.names = FALSE
)





# ============================================================================
# 20. SESSION / ANALYSIS SETTINGS SUMMARY
# ============================================================================

analysis_settings <- tibble::tibble(
  Setting = c(
    "Fixed-depth rarefaction",
    "Singleton filtering",
    "Primary depth adjustment",
    "Sensitivity depth adjustment",
    "iNEXT coverage rule",
    "iNEXT target coverage",
    "Primary model",
    "Sensitivity model",
    "Field handling",
    "Stage contrast",
    "Stage multiplicity correction"
  ),
  Value = c(
    "No",
    "No - global and sample singletons retained",
    "Coverage standardization with iNEXT.3D",
    "standardized log10 library size covariate",
    INEXT_COVERAGE_RULE,
    as.character(target_coverage),
    "Metric ~ Sampling * Treatment + (1|ExperimentalUnitID)",
    "Metric ~ Sampling * Treatment + LibrarySize_z + (1|ExperimentalUnitID)",
    "Separate model per Field",
    "Storage - Harvest estimated marginal mean contrast",
    "Holm within Analysis x Field"
  )
)

print(analysis_settings)

utils::write.csv(
  analysis_settings,
  "alpha_integrated_analysis_settings.csv",
  row.names = FALSE
)

capture.output(
  sessionInfo(),
  file = "alpha_integrated_sessionInfo.txt"
)
target_coverage
# ============================================================================
# END STEP 4
# ============================================================================


# STEP 5: BETA DIVERSITY -- REVIEWER-REVISED
# ==============================
# Reviewer comments addressed here:
# - no rarefaction;
# - Aitchison distance replaces Bray-Curtis for inferential community-composition analysis;
# - CLR is scale invariant, so CSS is not applied before Aitchison distance;
# - PERMANOVA uses marginal tests (by = "margin");
# - 9,999 permutations are used;
# - permutations respect block as the independent treatment experimental unit.

# 1) Genus-level counts and CLR transformation ---------------------------
tse__2 <- tse__fresh
tse__2 <- agglomerateByRank(tse__2, rank = "Genus", update.tree = TRUE)

counts_beta <- assay(tse__2, "counts")
keep_samples <- colSums(counts_beta) > 0
counts_beta <- counts_beta[, keep_samples, drop = FALSE]
tse__2 <- tse__2[, keep_samples]

# Remove extremely sparse genera before log-ratio analysis.
# This is a prevalence filter, NOT a library-size normalization.
min_prev <- 0.05
keep_taxa <- rowSums(counts_beta > 0) >= ceiling(min_prev * ncol(counts_beta))
counts_beta <- counts_beta[keep_taxa, , drop = FALSE]
tse__2 <- tse__2[keep_taxa, ]

# CLR requires strictly positive values. A small count-scale pseudocount is
# used only to define log-ratios for zeros; no rarefaction or CSS is applied.
clr_pseudocount <- 0.5
clr_counts <- apply(counts_beta, 2, function(x) {
  lx <- log(x + clr_pseudocount)
  lx - mean(lx)
})
clr_counts <- as.matrix(clr_counts)
rownames(clr_counts) <- rownames(counts_beta)
colnames(clr_counts) <- colnames(counts_beta)

# Euclidean distance in CLR space = Aitchison distance.
aitchison_dist <- stats::dist(t(clr_counts), method = "euclidean")

meta <- as.data.frame(colData(tse__2))
meta <- meta[colnames(clr_counts), , drop = FALSE]
meta$sampling  <- factor(meta$sampling, levels = c("Harvest", "Storage"))
meta$treatment <- factor(meta$treatment)
meta$field     <- factor(meta$field)
meta$block     <- factor(meta$block)

# A/B/C are replicate labels reused within each Treatment.
# The independent physical experimental unit is Field x Treatment x Block.
meta$ExperimentalUnitID <- interaction(
  meta$field,
  meta$treatment,
  meta$block,
  drop = TRUE
)

# 2) Constrained ordination in Aitchison space --------------------------
dbrda_all <- vegan::capscale(
  aitchison_dist ~ sampling + treatment + field,
  data = meta
)

dbrda_scores <- as.data.frame(vegan::scores(dbrda_all, display = "sites"))
dbrda_scores$Sampling  <- meta$sampling
dbrda_scores$Treatment <- meta$treatment
dbrda_scores$Field     <- meta$field

fields <- levels(meta$field)
sampling_levels <- levels(meta$sampling)

# 3) Experimental-unit-aware marginal PERMANOVA + PERMDISP ----------------
PERMUTATIONS_N <- 9999L

# Experimental design within each field
# -------------------------------------
# A/B/C are replicate labels reused within each Treatment:
#
#   Control: A, B, C
#   Geoxe:   A, B, C
#   Ulmasud: A, B, C
#
# Thus "A" under Control and "A" under Geoxe are DIFFERENT physical units.
# The independent unit is:
#
#   ExperimentalUnitID = Field x Treatment x Block
#
# Multiple apples are sampled from each experimental unit at Harvest and
# Storage.
#
# Sampling is therefore a WITHIN-unit factor.
# Treatment is a BETWEEN-unit factor.
#
# Sampling and Treatment require different permutation restrictions. To retain
# marginal tests (by = "margin"), the same model is fitted twice:
#
#   community ~ sampling + treatment
#
# - Sampling P-value: observations permuted within ExperimentalUnitID.
# - Treatment P-value: complete ExperimentalUnitID groups permuted.
#
# Only the valid marginal row for each effect is retained in the combined
# PERMANOVA table.


validate_experimental_unit_design <- function(dat) {
  
  dat <- droplevels(dat)
  
  required_cols <- c(
    "field",
    "treatment",
    "block",
    "sampling",
    "ExperimentalUnitID"
  )
  
  missing_cols <- setdiff(
    required_cols,
    colnames(dat)
  )
  
  if (length(missing_cols) > 0L) {
    stop(
      "Missing metadata columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # One treatment per physical experimental unit.
  unit_treatment <- dat %>%
    dplyr::distinct(
      ExperimentalUnitID,
      treatment
    ) %>%
    dplyr::count(
      ExperimentalUnitID,
      name = "N_treatments"
    )
  
  if (any(unit_treatment$N_treatments != 1L)) {
    stop(
      "At least one ExperimentalUnitID is associated with more than one Treatment."
    )
  }
  
  # Replication at the physical-unit level.
  treatment_replication <- dat %>%
    dplyr::distinct(
      ExperimentalUnitID,
      treatment
    ) %>%
    dplyr::count(
      treatment,
      name = "N_independent_units"
    )
  
  if (any(treatment_replication$N_independent_units < 2L)) {
    stop(
      "At least one Treatment has fewer than two independent experimental units."
    )
  }
  
  invisible(treatment_replication)
}


# -------------------------------------------------------------------------
# Sampling permutation
#
# Harvest/Storage labels are permuted only among observations belonging to
# the same physical ExperimentalUnitID.
# -------------------------------------------------------------------------

make_sampling_permutation <- function(
    dat,
    nperm = PERMUTATIONS_N
) {
  
  dat <- droplevels(dat)
  
  ctrl <- permute::how(
    blocks = dat$ExperimentalUnitID,
    nperm = nperm
  )
  
  possible <- permute::numPerms(
    seq_len(nrow(dat)),
    control = ctrl
  )
  
  message(
    "Available within-unit permutations for Sampling: ",
    possible
  )
  
  ctrl
}


# -------------------------------------------------------------------------
# Treatment permutation
#
# Treatment is assigned between independent experimental units. Therefore
# complete ExperimentalUnitID groups are permuted while their internal sample
# structure is kept intact.
#
# The plot-level permutation machinery in permute is used here only to move
# complete physical units. There is NO additional block stratum because A/B/C
# are replicate labels reused independently within each treatment.
# -------------------------------------------------------------------------

make_treatment_permutation <- function(
    dat,
    nperm = PERMUTATIONS_N
) {
  
  dat <- droplevels(dat)
  
  validate_experimental_unit_design(dat)
  
  unit_n <- table(
    dat$ExperimentalUnitID
  )
  
  if (length(unique(as.numeric(unit_n))) != 1L) {
    
    print(unit_n)
    
    stop(
      "Whole-unit Treatment permutation requires equal numbers of observations ",
      "per ExperimentalUnitID."
    )
  }
  
  # For the combined Harvest + Storage analysis, verify that all units have
  # the same sampling-stage layout.
  if (
    "sampling" %in% colnames(dat) &&
    nlevels(dat$sampling) > 1L
  ) {
    
    unit_sampling <- table(
      dat$ExperimentalUnitID,
      dat$sampling
    )
    
    reference_row <- as.numeric(
      unit_sampling[1, ]
    )
    
    same_layout <- apply(
      unit_sampling,
      1,
      function(x) {
        identical(
          as.numeric(x),
          reference_row
        )
      }
    )
    
    if (!all(same_layout)) {
      
      print(unit_sampling)
      
      stop(
        "Whole-unit Treatment permutation requires the same Harvest/Storage ",
        "sample layout in every ExperimentalUnitID."
      )
    }
  }
  
  ctrl <- permute::how(
    within = permute::Within(
      type = "none"
    ),
    plots = permute::Plots(
      strata = dat$ExperimentalUnitID,
      type = "free"
    ),
    nperm = nperm
  )
  
  possible <- permute::numPerms(
    seq_len(nrow(dat)),
    control = ctrl
  )
  
  message(
    "Available whole-unit permutations for Treatment: ",
    possible
  )
  
  ctrl
}


# -------------------------------------------------------------------------
# Combine valid marginal rows
# -------------------------------------------------------------------------

combine_marginal_permanova <- function(
    sampling_fit,
    treatment_fit
) {
  
  sampling_tab <- as.data.frame(
    sampling_fit
  )
  
  treatment_tab <- as.data.frame(
    treatment_fit
  )
  
  rbind(
    sampling =
      sampling_tab[
        "sampling",
        ,
        drop = FALSE
      ],
    treatment =
      treatment_tab[
        "treatment",
        ,
        drop = FALSE
      ],
    Residual =
      sampling_tab[
        "Residual",
        ,
        drop = FALSE
      ],
    Total =
      sampling_tab[
        "Total",
        ,
        drop = FALSE
      ]
  )
}


# =========================================================================
# 3A) OVERALL ANALYSIS WITHIN EACH FIELD
#
# community ~ sampling + treatment
#
# Sampling:
#   marginal effect, permutations within ExperimentalUnitID
#
# Treatment:
#   marginal effect, whole ExperimentalUnitID permutations
# =========================================================================

field_results <- list()

for (f in fields) {
  
  idx <- meta$field == f
  
  meta_sub <- droplevels(
    meta[
      idx,
      ,
      drop = FALSE
    ]
  )
  
  # Keep observations belonging to a physical unit contiguous and give all
  # units the same internal Harvest/Storage ordering.
  ord <- order(
    meta_sub$ExperimentalUnitID,
    meta_sub$sampling,
    rownames(meta_sub)
  )
  
  meta_sub <- meta_sub[
    ord,
    ,
    drop = FALSE
  ]
  
  clr_sub <- clr_counts[
    ,
    rownames(meta_sub),
    drop = FALSE
  ]
  
  dist_sub <- stats::dist(
    t(clr_sub),
    method = "euclidean"
  )
  
  validate_experimental_unit_design(
    meta_sub
  )
  
  
  # ---------------- Sampling effect ----------------
  
  perm_sampling <- make_sampling_permutation(
    meta_sub
  )
  
  permanova_sampling_full <- vegan::adonis2(
    dist_sub ~ sampling + treatment,
    data = meta_sub,
    permutations = perm_sampling,
    by = "margin"
  )
  
  
  # ---------------- Treatment effect ----------------
  
  perm_treatment <- make_treatment_permutation(
    meta_sub
  )
  
  permanova_treatment_full <- vegan::adonis2(
    dist_sub ~ sampling + treatment,
    data = meta_sub,
    permutations = perm_treatment,
    by = "margin"
  )
  
  
  # One table with the valid P-value for each factor.
  permanova_res <- combine_marginal_permanova(
    sampling_fit = permanova_sampling_full,
    treatment_fit = permanova_treatment_full
  )
  
  
  # ---------------- PERMDISP: Sampling ----------------
  
  permdisp_sampling <- vegan::betadisper(
    dist_sub,
    meta_sub$sampling
  )
  
  permdisp_sampling_test <- vegan::permutest(
    permdisp_sampling,
    permutations = perm_sampling
  )
  
  
  # ---------------- PERMDISP: Treatment ----------------
  
  permdisp_treatment <- vegan::betadisper(
    dist_sub,
    meta_sub$treatment
  )
  
  permdisp_treatment_test <- vegan::permutest(
    permdisp_treatment,
    permutations = perm_treatment
  )
  
  
  field_results[[f]] <- list(
    permanova = permanova_res,
    permanova_sampling_full = permanova_sampling_full,
    permanova_treatment_full = permanova_treatment_full,
    permdisp_sampling = permdisp_sampling_test,
    permdisp_treatment = permdisp_treatment_test,
    metadata = meta_sub
  )
}


# =========================================================================
# 3B) TREATMENT EFFECT WITHIN EACH SAMPLING STAGE
#
# Treatment remains a between-unit effect when Harvest and Storage are
# analysed separately, so complete ExperimentalUnitID groups are permuted.
# =========================================================================

field_sampling_results <- list()

for (f in fields) {
  
  for (smp in sampling_levels) {
    
    key <- paste(
      f,
      smp,
      sep = "_"
    )
    
    idx <-
      meta$field == f &
      meta$sampling == smp
    
    meta_sub <- droplevels(
      meta[
        idx,
        ,
        drop = FALSE
      ]
    )
    
    if (
      nrow(meta_sub) < 4L ||
      nlevels(meta_sub$treatment) < 2L
    ) {
      next
    }
    
    ord <- order(
      meta_sub$ExperimentalUnitID,
      rownames(meta_sub)
    )
    
    meta_sub <- meta_sub[
      ord,
      ,
      drop = FALSE
    ]
    
    clr_sub <- clr_counts[
      ,
      rownames(meta_sub),
      drop = FALSE
    ]
    
    dist_sub <- stats::dist(
      t(clr_sub),
      method = "euclidean"
    )
    
    validate_experimental_unit_design(
      meta_sub
    )
    
    perm_treatment <- make_treatment_permutation(
      meta_sub
    )
    
    permanova_res <- vegan::adonis2(
      dist_sub ~ treatment,
      data = meta_sub,
      permutations = perm_treatment,
      by = "margin"
    )
    
    permdisp_treatment <- vegan::betadisper(
      dist_sub,
      meta_sub$treatment
    )
    
    permdisp_treatment_test <- vegan::permutest(
      permdisp_treatment,
      permutations = perm_treatment
    )
    
    field_sampling_results[[key]] <- list(
      permanova = permanova_res,
      permdisp_treatment = permdisp_treatment_test,
      metadata = meta_sub
    )
  }
}


# =========================================================================
# 3C) PRINT RESULTS
# =========================================================================

# Pfatten/Vadena overall
field_results[["Pfatten/Vadena"]]$permanova
field_results[["Pfatten/Vadena"]]$permdisp_sampling
field_results[["Pfatten/Vadena"]]$permdisp_treatment

# Sinich/Sinigo overall
field_results[["Sinich/Sinigo"]]$permanova
field_results[["Sinich/Sinigo"]]$permdisp_sampling
field_results[["Sinich/Sinigo"]]$permdisp_treatment

# Pfatten/Vadena within Sampling
field_sampling_results[["Pfatten/Vadena_Harvest"]]$permanova
field_sampling_results[["Pfatten/Vadena_Harvest"]]$permdisp_treatment

field_sampling_results[["Pfatten/Vadena_Storage"]]$permanova
field_sampling_results[["Pfatten/Vadena_Storage"]]$permdisp_treatment

# Sinich/Sinigo within Sampling
field_sampling_results[["Sinich/Sinigo_Harvest"]]$permanova
field_sampling_results[["Sinich/Sinigo_Harvest"]]$permdisp_treatment

field_sampling_results[["Sinich/Sinigo_Storage"]]$permanova
field_sampling_results[["Sinich/Sinigo_Storage"]]$permdisp_treatment


saveRDS(
  field_results,
  "PERMANOVA_Aitchison_ExperimentalUnit_by_field.rds"
)

saveRDS(
  field_sampling_results,
  "PERMANOVA_Aitchison_ExperimentalUnit_field_sampling.rds"
)


# 4) Ordination ellipses -------------------------------------------------
vegan_ellipse <- function(df, group_col, axes = c("CAP1", "CAP2"), level = 0.68) {
  groups <- unique(df[[group_col]])
  ellipses <- lapply(groups, function(g) {
    sub <- df[df[[group_col]] == g, axes, drop = FALSE]
    if (nrow(sub) < 3) return(NULL)
    cov_mat <- cov(sub)
    center <- colMeans(sub)
    angles <- seq(0, 2 * pi, length.out = 100)
    ellipse <- t(sapply(angles, function(theta) {
      center + sqrt(qchisq(level, 2)) * t(chol(cov_mat)) %*% c(cos(theta), sin(theta))
    }))
    ellipse <- as.data.frame(ellipse)
    ellipse$Group <- g
    colnames(ellipse)[1:2] <- axes
    ellipse
  })
  bind_rows(ellipses)
}

dbrda_ellipses <- vegan_ellipse(dbrda_scores, group_col = "Sampling", level = 0.68)

dbrda_scores$Treatment <- dplyr::recode(
  as.character(dbrda_scores$Treatment),
  "Control" = "Control",
  "Geoxe" = "Fludioxonil",
  "Ulmasud" = "Acidic clays"
)

dbrda_scores$Treatment <- factor(
  dbrda_scores$Treatment,
  levels = c("Control", "Fludioxonil", "Acidic clays")
)

dbrda_plot <- ggplot(dbrda_scores, aes(x = CAP1, y = CAP2)) +
  geom_point(aes(color = Sampling, shape = Treatment), size = 3) +
  geom_path(
    data = dbrda_ellipses,
    aes(x = CAP1, y = CAP2, color = Group),
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_color_manual(values = c(Harvest = "#E64B35", Storage = "#4DBBD5")) +
  scale_shape_manual(values = c(
    "Control" = 16,
    "Fludioxonil" = 17,
    "Acidic clays" = 3
  )) +
  labs(
    x = "CAP1",
    y = "CAP2",
    color = "Sampling",
    shape = "Treatment",
    subtitle = "Aitchison distance (CLR-transformed genus counts)"
  ) +
  facet_wrap(~Field) +
  theme_nature +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )

dbrda_plot_no_legend <- dbrda_plot + theme(legend.position = "none")
dbrda_plot

# ============================================================
# 2) PCA OF CLR-TRANSFORMED ABUNDANCES
# ============================================================

# Samples x genera
clr_samples <- t(clr_counts)

# PCA in CLR space
# center = TRUE: standard PCA centering across samples
# scale. = FALSE: do NOT rescale genera, as this would alter
# Aitchison geometry
pca_all <- stats::prcomp(
  clr_samples,
  center = TRUE,
  scale. = FALSE
)

# Extract sample scores
pca_scores <- as.data.frame(pca_all$x)

# Add metadata in exactly the same sample order
pca_scores$Sampling <- meta$sampling
pca_scores$Treatment <- meta$treatment
pca_scores$Field <- meta$field
pca_scores$Block <- meta$block
pca_scores$ExperimentalUnitID <- meta$ExperimentalUnitID

# Percentage of variance explained
pca_var <- 100 * pca_all$sdev^2 / sum(pca_all$sdev^2)

PC1_lab <- paste0(
  "PC1 (",
  round(pca_var[1], 1),
  "%)"
)

PC2_lab <- paste0(
  "PC2 (",
  round(pca_var[2], 1),
  "%)"
)


# ============================================================
# Treatment labels
# ============================================================

pca_scores <- pca_scores %>%
  dplyr::mutate(
    Treatment_label = dplyr::case_when(
      
      Field == "Pfatten/Vadena" &
        Treatment == "Control" ~
        "Control",
      
      Field == "Pfatten/Vadena" &
        Treatment == "Geoxe" ~
        "Fludioxonil",
      
      Field == "Pfatten/Vadena" &
        Treatment == "Ulmasud" ~
        "Acidic clays",
      
      Field == "Sinich/Sinigo" &
        Treatment == "Control" ~
        "Control",
      
      Field == "Sinich/Sinigo" &
        Treatment == "Geoxe" ~
        "Captan + Fludioxonil",
      
      Field == "Sinich/Sinigo" &
        Treatment == "Ulmasud" ~
        "Captan + Acidic clays",
      
      TRUE ~ as.character(Treatment)
    ),
    
    Sampling = factor(
      Sampling,
      levels = c("Harvest", "Storage")
    )
  )


# ============================================================
# PCA plot
# ============================================================

pca_plot <- ggplot2::ggplot(
  pca_scores,
  ggplot2::aes(
    x = PC1,
    y = PC2,
    color = Sampling,
    shape = Treatment_label
  )
) +
  
  ggplot2::geom_point(
    size = 3,
    alpha = 0.8
  ) +
  
  # Ellipses for sampling stage
  ggplot2::stat_ellipse(
    ggplot2::aes(
      group = Sampling,
      color = Sampling
    ),
    type = "t",
    level = 0.68,
    linewidth = 1,
    linetype = "dashed",
    show.legend = FALSE
  ) +
  
  ggplot2::facet_wrap(
    ~ Field
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      Harvest = "#E64B35",
      Storage = "#4DBBD5"
    )
  ) +
  
  ggplot2::labs(
    x = PC1_lab,
    y = PC2_lab,
    color = "Sampling",
    shape = "Treatment"
  ) +
  
  theme_nature


pca_plot

legend_only <- cowplot::get_legend(
  pca_plot +
    ggplot2::theme(
      legend.position = "right"
    )
)

legend_only

pca_plot_nolegend <- pca_plot +
  ggplot2::theme(
    legend.position = "none"
  )

pca_plot_nolegend

# ============================================================================
# 18. MAIN FIGURE -- SHANNON EFFECTIVE DIVERSITY
# ============================================================================
# Primary coverage-standardized alpha-diversity analysis.
#
# TD_q1 = Shannon effective diversity (Hill number q = 1).
#
# Delta is NOT calculated from the plotted raw/sample values.
# Delta is the estimated marginal mean contrast from the primary mixed model:
#
#              Delta = Storage - Harvest
#
# Therefore:
#   Delta > 0  -> higher diversity after storage
#   Delta < 0  -> lower diversity after storage
#
# P-values shown are Holm-adjusted post-hoc P-values.
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Plot data
# ----------------------------------------------------------------------------

shannon_plot_data <- alpha_primary %>%
  dplyr::select(
    sam_name,
    Field,
    Sampling,
    TD_q1
  ) %>%
  dplyr::filter(
    !is.na(TD_q1)
  ) %>%
  dplyr::mutate(
    Sampling = factor(
      Sampling,
      levels = c("Harvest", "Storage")
    ),
    Sampling_x = as.numeric(Sampling)
  )


# ----------------------------------------------------------------------------
# 2. Post-hoc contrast used as effect size
# ----------------------------------------------------------------------------

shannon_posthoc <- stage_overall_all %>%
  dplyr::filter(
    Analysis == "Primary_coverage_standardized",
    Metric == "TD_q1",
    contrast == "Storage - Harvest"
  ) %>%
  dplyr::mutate(
    
    p_label = dplyr::case_when(
      p.value.Holm < 0.001 ~ "p < 0.001",
      TRUE ~ paste0(
        "p = ",
        formatC(
          p.value.Holm,
          format = "f",
          digits = 3
        )
      )
    ),
    
    label = paste0(
      "\u0394 = ",
      formatC(
        estimate,
        format = "f",
        digits = 3
      ),
      "\n",
      p_label
    ),
    
    # Put annotation between Harvest and Storage
    Sampling_x = 1.5
  )


# ----------------------------------------------------------------------------
# 3. Figure
# ----------------------------------------------------------------------------

p_alpha_shannon <- ggplot2::ggplot(
  shannon_plot_data,
  ggplot2::aes(
    x = Sampling_x,
    y = TD_q1,
    fill = Sampling,
    color = Sampling
  )
) +
  
  # Violin distribution
  ggplot2::geom_violin(
    ggplot2::aes(
      group = Sampling
    ),
    width = 0.75,
    trim = FALSE,
    alpha = 0.45,
    linewidth = 0.5
  ) +
  
  # Boxplot inside violin
  ggplot2::geom_boxplot(
    ggplot2::aes(
      group = Sampling
    ),
    width = 0.16,
    outlier.shape = NA,
    alpha = 0.75,
    linewidth = 0.5
  ) +
  
  # Individual apples
  ggplot2::geom_point(
    position = ggplot2::position_jitter(
      width = 0.08,
      height = 0
    ),
    size = 1.7,
    alpha = 0.65
  ) +
  
  # Delta from EMM post-hoc contrast
  ggplot2::geom_text(
    data = shannon_posthoc,
    ggplot2::aes(
      x = Sampling_x,
      y = Inf,
      label = label
    ),
    inherit.aes = FALSE,
    vjust = 1.15,
    size = 3.6
  ) +
  
  # One panel per field
  ggplot2::facet_wrap(
    ~ Field,
    nrow = 1,
    scales = "free_y"
  ) +
  
  ggplot2::scale_x_continuous(
    breaks = c(1, 2),
    labels = c("Harvest", "Storage"),
    limits = c(0.6, 2.4)
  ) +
  
  ggplot2::scale_fill_manual(
    values = c(
      Harvest = "#E64B35",
      Storage = "#4DBBD5"
    )
  ) +
  
  ggplot2::scale_color_manual(
    values = c(
      Harvest = "#E64B35",
      Storage = "#4DBBD5"
    )
  ) +
  
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(
      mult = c(0.05, 0.18)
    )
  ) +
  
  ggplot2::labs(
    x = NULL,
    y = "Shannon (q = 1)"
  ) +
  
  theme_alpha +
  
  ggplot2::theme(
    legend.position = "none",
    
    strip.text = ggplot2::element_text(
      face = "bold",
      size = 12
    ),
    
    axis.text.x = ggplot2::element_text(
      size = 10
    ),
    
    axis.title.y = ggplot2::element_text(
      size = 11
    )
  )


print(p_alpha_shannon)

p_alpha_shannon <- p_alpha_shannon + theme_nature + ggplot2::theme(
  legend.position = "none")
p_alpha_shannon


# 5) Combine with Shannon ------------------------------------------------
p_shannon <- p_shannon + theme(axis.title.x = element_blank())

combined_plot <- plot_grid(
  p_alpha_shannon,
  pca_plot_nolegend,
  nrow = 2,
  labels = c(),
  rel_heights = c(1.2, 1.8),
  label_size = 15
)

combined_plot
saveRDS(combined_plot, file = "combined_ITS_plot_reviewer_revised.rds")

legend <- get_legend(pca_plot + theme(legend.position = "bottom"))
saveRDS(legend, file = "combined_ITS_legend_reviewer_revised.rds")

# ==============================
# STEP 6: BETA PARTITIONING 
# ==============================
pseq.rel <- microbiome::transform(physeq_asv_filtered, "compositional")

otu_pa <- as(otu_table(pseq.rel), "matrix")
otu_pa <- 1 * (otu_pa > 0)

if (taxa_are_rows(pseq.rel)) {
  otu_pa <- t(otu_pa)
}

meta <- data.frame(sample_data(pseq.rel))

beta_summary_list <- list()
mantel_list <- list()
beta_long_list <- list()

fields <- unique(meta$field)
treatments <- unique(meta$treatment)

for (f in fields) {
  
  for (t in treatments) {
    
    idx <- meta$field == f & meta$treatment == t
    
    if (sum(idx) < 4) next
    
    otu_sub <- otu_pa[idx, ]
    meta_sub <- meta[idx, ]
    
    stage <- meta_sub$sampling
    
    # ---------------------------
    # Beta partitioning
    # ---------------------------
    
    beta <- beta.pair(otu_sub, index.family = "sorensen")
    
    beta_sor <- as.matrix(beta$beta.sor)
    beta_sim <- as.matrix(beta$beta.sim)
    beta_nes <- as.matrix(beta$beta.sne)
    
    # indices
    H <- which(stage == "Harvest")
    S <- which(stage == "Storage")
    
    # remove diagonal comparisons
    within_H <- beta_sor[H, H][upper.tri(beta_sor[H, H])]
    within_S <- beta_sor[S, S][upper.tri(beta_sor[S, S])]
    
    between <- beta_sor[H, S]
    turnover <- beta_sim[H, S]
    nestedness <- beta_nes[H, S]
    
    key <- paste(f, t, sep = "_")
    
    # ---------------------------
    # Summary values
    # ---------------------------
    
    beta_summary_list[[key]] <- data.frame(
      Field = f,
      Treatment = t,
      betaSOR_within_Harvest = mean(within_H, na.rm = TRUE),
      betaSOR_within_Storage = mean(within_S, na.rm = TRUE),
      betaSOR_between = mean(between, na.rm = TRUE),
      betaSIM_between = mean(turnover, na.rm = TRUE),
      betaNES_between = mean(nestedness, na.rm = TRUE)
    )
    
    # ---------------------------
    # Mantel tests
    # ---------------------------
    
    stage_numeric <- ifelse(stage == "Harvest", 0, 1)
    stage_dist <- dist(stage_numeric)
    
    mantel_sor <- mantel(beta$beta.sor, stage_dist, permutations = 999)
    mantel_sim <- mantel(beta$beta.sim, stage_dist, permutations = 999)
    mantel_nes <- mantel(beta$beta.sne, stage_dist, permutations = 999)
    
    mantel_list[[key]] <- data.frame(
      Field = f,
      Treatment = t,
      Component = c("βSOR", "βSIM", "βNES"),
      R = c(
        mantel_sor$statistic,
        mantel_sim$statistic,
        mantel_nes$statistic
      ),
      P_value = c(
        mantel_sor$signif,
        mantel_sim$signif,
        mantel_nes$signif
      )
    )
    
    # ---------------------------
    # Extract pairwise distances
    # ---------------------------
    
    extract_pairs <- function(mat, stage_vector, stage_name) {
      
      idx_stage <- which(stage_vector == stage_name)
      vals <- mat[idx_stage, idx_stage]
      vals <- vals[upper.tri(vals)]
      
      data.frame(
        Stage = stage_name,
        Value = vals
      )
    }
    
    beta_long <- bind_rows(
      extract_pairs(beta_sor, stage, "Harvest") %>% mutate(Component = "βSOR"),
      extract_pairs(beta_sim, stage, "Harvest") %>% mutate(Component = "βSIM"),
      extract_pairs(beta_nes, stage, "Harvest") %>% mutate(Component = "βNES"),
      extract_pairs(beta_sor, stage, "Storage") %>% mutate(Component = "βSOR"),
      extract_pairs(beta_sim, stage, "Storage") %>% mutate(Component = "βSIM"),
      extract_pairs(beta_nes, stage, "Storage") %>% mutate(Component = "βNES")
    )
    
    beta_long$Field <- f
    beta_long$Treatment <- t
    
    beta_long_list[[key]] <- beta_long
  }
}

beta_summary <- bind_rows(beta_summary_list)
mantel_df <- bind_rows(mantel_list)
beta_long_all <- bind_rows(beta_long_list)

# -------------------------------------------------
# Recode treatment labels for all plots
# -------------------------------------------------
beta_long_all <- beta_long_all %>%
  mutate(
    Treatment_label = case_when(
      Field == "Pfatten/Vadena" & Treatment == "Control" ~ "Control",
      Field == "Pfatten/Vadena" & Treatment == "Geoxe" ~ "Fludioxonil",
      Field == "Pfatten/Vadena" & Treatment == "Ulmasud" ~ "Acidic Clays",
      Field == "Sinich/Sinigo" & Treatment == "Control" ~ "Control",
      Field == "Sinich/Sinigo" & Treatment == "Geoxe" ~ "Captan + Fludioxonil",
      Field == "Sinich/Sinigo" & Treatment == "Ulmasud" ~ "Captan + Acidic Clays",
      TRUE ~ Treatment
    )
  )

mantel_df <- mantel_df %>%
  mutate(
    Treatment_label = case_when(
      Field == "Pfatten/Vadena" & Treatment == "Control" ~ "Control",
      Field == "Pfatten/Vadena" & Treatment == "Geoxe" ~ "Fludioxonil",
      Field == "Pfatten/Vadena" & Treatment == "Ulmasud" ~ "Acidic Clays",
      Field == "Sinich/Sinigo" & Treatment == "Control" ~ "Control",
      Field == "Sinich/Sinigo" & Treatment == "Geoxe" ~ "Captan + Fludioxonil",
      Field == "Sinich/Sinigo" & Treatment == "Ulmasud" ~ "Captan + Acidic Clays",
      TRUE ~ Treatment
    )
  )

beta_summary <- beta_summary %>%
  mutate(
    Treatment_label = case_when(
      Field == "Pfatten/Vadena" & Treatment == "Control" ~ "Control",
      Field == "Pfatten/Vadena" & Treatment == "Geoxe" ~ "Fludioxonil",
      Field == "Pfatten/Vadena" & Treatment == "Ulmasud" ~ "Acidic Clays",
      Field == "Sinich/Sinigo" & Treatment == "Control" ~ "Control",
      Field == "Sinich/Sinigo" & Treatment == "Geoxe" ~ "Captan + Fludioxonil",
      Field == "Sinich/Sinigo" & Treatment == "Ulmasud" ~ "Captan + Acidic Clays",
      TRUE ~ Treatment
    )
  )

# -------------------------------------------------
# Keep facet order
# -------------------------------------------------
treatment_levels <- c(
  "Control",
  "Fludioxonil",
  "Acidic Clays",
  "Captan + Fludioxonil",
  "Captan + Acidic Clays"
)

beta_long_all$Treatment_label <- factor(
  beta_long_all$Treatment_label,
  levels = treatment_levels
)

mantel_df$Treatment_label <- factor(
  mantel_df$Treatment_label,
  levels = treatment_levels
)

beta_summary$Treatment_label <- factor(
  beta_summary$Treatment_label,
  levels = treatment_levels
)

beta_means <- beta_long_all %>%
  group_by(Field, Treatment, Treatment_label, Component, Stage) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

beta_means

mantel_labels <- mantel_df %>%
  mutate(
    label = paste0(
      "R = ", round(R, 2),
      "\np = ", signif(P_value, 2)
    )
  )

label_positions <- beta_long_all %>%
  group_by(Field, Treatment, Treatment_label, Component) %>%
  summarise(
    Stage = "Storage",
    Value = max(Value, na.rm = TRUE) * 0.95,
    .groups = "drop"
  )

mantel_labels <- dplyr::left_join(
  mantel_labels,
  label_positions,
  by = c("Field", "Treatment", "Treatment_label", "Component")
)

p_beta <- ggplot(beta_long_all, aes(Stage, Value)) +
  
  geom_jitter(
    aes(color = Stage),
    width = 0.08,
    alpha = 0.35,
    size = 1.5
  ) +
  
  geom_line(
    data = beta_means,
    aes(group = interaction(Field, Treatment, Component)),
    linewidth = 1.1,
    color = "black"
  ) +
  
  geom_point(
    data = beta_means,
    size = 4,
    shape = 21,
    fill = "white",
    color = "black",
    stroke = 1.2
  ) +
  
  geom_text(
    data = mantel_labels,
    aes(Stage, Value, label = label),
    inherit.aes = FALSE,
    size = 3.5,
    fontface = "bold",
    hjust = 1
  ) +
  
  facet_grid(Field + Treatment_label ~ Component, drop = TRUE) +
  
  scale_color_manual(values = c(
    Harvest = "#E64B35",
    Storage = "#4DBBD5"
  )) +
  
  labs(
    x = NULL,
    y = "Beta diversity"
  ) +
  
  theme_nature

p_beta

# Filter beta_long_all for βSIM component
beta_sim_long <- beta_long_all %>%
  filter(Component == "βSIM")

# Compute mean values per Field × Treatment × Stage
beta_sim_means <- beta_sim_long %>%
  group_by(Field, Treatment, Treatment_label, Stage) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Prepare Mantel labels for βSIM only
mantel_sim_labels <- mantel_df %>%
  filter(Component == "βSIM") %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\np = ", signif(P_value, 2))
  )

# Position labels at 95% of max value per Field × Treatment
label_positions <- beta_sim_long %>%
  group_by(Field, Treatment, Treatment_label) %>%
  summarise(
    Stage = "Storage",
    Value = max(Value, na.rm = TRUE) * 0.95,
    .groups = "drop"
  )

mantel_sim_labels <- dplyr::left_join(
  mantel_sim_labels,
  label_positions,
  by = c("Field", "Treatment", "Treatment_label")
)

beta_pfatten <- beta_sim_long %>%
  dplyr::filter(Field == "Pfatten/Vadena")

means_pfatten <- beta_sim_means %>%
  dplyr::filter(Field == "Pfatten/Vadena")

mantel_pfatten <- mantel_sim_labels %>%
  dplyr::filter(Field == "Pfatten/Vadena")

p_beta_pfatten <- ggplot(beta_pfatten, aes(Stage, Value)) +
  
  geom_jitter(
    aes(color = Stage),
    width = 0.05,
    alpha = 0.5,
    size = 2
  ) +
  
  geom_line(
    data = means_pfatten,
    aes(group = Treatment),
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = means_pfatten,
    shape = 21,
    size = 4,
    fill = "white",
    color = "black",
    stroke = 1.2
  ) +
  
  geom_text(
    data = mantel_pfatten,
    aes(Stage, Value, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 7,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ Treatment_label,
    nrow = 1,
    drop = TRUE
  ) +
  
  scale_color_manual(values = c(
    Harvest = "#E64B35",
    Storage = "#4DBBD5"
  )) +
  
  labs(
    title = "Pfatten/Vadena",
    x = NULL,
    y = expression(beta[SIM] ~ "(turnover)")
  ) +
  
  theme_nature +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  )

p_beta_pfatten

beta_sinich <- beta_sim_long %>%
  dplyr::filter(Field == "Sinich/Sinigo")

means_sinich <- beta_sim_means %>%
  dplyr::filter(Field == "Sinich/Sinigo")

mantel_sinich <- mantel_sim_labels %>%
  dplyr::filter(Field == "Sinich/Sinigo")

p_beta_sinich <- ggplot(beta_sinich, aes(Stage, Value)) +
  
  geom_jitter(
    aes(color = Stage),
    width = 0.05,
    alpha = 0.5,
    size = 2
  ) +
  
  geom_line(
    data = means_sinich,
    aes(group = Treatment),
    color = "black",
    linewidth = 1
  ) +
  
  geom_point(
    data = means_sinich,
    shape = 21,
    size = 4,
    fill = "white",
    color = "black",
    stroke = 1.2
  ) +
  
  geom_text(
    data = mantel_sinich,
    aes(Stage, Value, label = label),
    inherit.aes = FALSE,
    hjust = 1,
    size = 7,
    fontface = "bold"
  ) +
  
  facet_wrap(
    ~ Treatment_label,
    nrow = 1,
    drop = TRUE
  ) +
  
  scale_color_manual(values = c(
    Harvest = "#E64B35",
    Storage = "#4DBBD5"
  )) +
  
  theme_nature +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    title = "Sinich/Sinigo",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_beta_sinich

p_beta_final <- p_beta_pfatten + p_beta_sinich +
  plot_layout(ncol = 2, guides = "collect") +
  theme(legend.position = "none")

p_beta_final

# =========================================================
# STEP 7: LEFSE
# Harvest vs Storage within each Field x Treatment
# CSS normalization performed INSIDE run_lefse()
# =========================================================

library(microbiomeMarker)
library(phyloseq)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(patchwork)
library(openxlsx)


# =========================================================
# SETTINGS
# =========================================================

fields_order <- c(
  "Pfatten/Vadena",
  "Sinich/Sinigo"
)

treatments_order <- c(
  "Control",
  "Geoxe",
  "Ulmasud"
)

treatment_labels <- c(
  "Pfatten/Vadena__Control" = "Control",
  "Pfatten/Vadena__Geoxe"   = "Fludioxonil",
  "Pfatten/Vadena__Ulmasud" = "Acidic Clays",
  
  "Sinich/Sinigo__Control"  = "Control",
  "Sinich/Sinigo__Geoxe"    = "Captan + Fludioxonil",
  "Sinich/Sinigo__Ulmasud"  = "Captan + Acidic Clays"
)


# =========================================================
# LEfSe PARAMETERS
# =========================================================

lefse_kw_cutoff <- 0.05
lefse_wilcoxon_cutoff <- 0.05
lefse_lda_cutoff <- 2

lefse_bootstrap_n <- 100
lefse_bootstrap_fraction <- 2/3
lefse_sample_min <- 5


# =========================================================
# START FROM RAW COUNTS
#
# IMPORTANT:
#
# NO prevalence filtering
# NO minimum abundance filtering
# NO relative abundance transformation
#
# Only:
#   1. subset Field x Treatment
#   2. genus aggregation
#   3. remove zero-total genera
#   4. remove invariant genera
#   5. CSS normalization inside run_lefse()
#   6. Harvest vs Storage LEfSe
# =========================================================

ps_lefse <- physeq_asv_filtered


# =========================================================
# PREPARE METADATA
# =========================================================

meta_lefse <- as.data.frame(
  phyloseq::sample_data(ps_lefse)
)

meta_lefse$sampling <- factor(
  meta_lefse$sampling,
  levels = c(
    "Harvest",
    "Storage"
  )
)

meta_lefse$field <- factor(
  meta_lefse$field,
  levels = fields_order
)

meta_lefse$treatment <- factor(
  meta_lefse$treatment,
  levels = treatments_order
)

phyloseq::sample_data(ps_lefse) <-
  phyloseq::sample_data(meta_lefse)


# =========================================================
# HELPER:
# CONVERT microbiomeMarker marker_table TO NORMAL DATA FRAME
#
# marker_table is an S4 class inheriting from data.frame.
# We explicitly extract its .Data slot so dplyr never
# operates directly on the S4 object.
# =========================================================

marker_table_to_df <- function(mm) {
  
  mt_raw <- microbiomeMarker::marker_table(mm)
  
  mt <- base::as.data.frame(
    methods::slot(
      mt_raw,
      ".Data"
    ),
    stringsAsFactors = FALSE
  )
  
  colnames(mt) <- colnames(
    mt_raw
  )
  
  mt
}


# =========================================================
# HELPER:
# EXTRACT AND FORMAT LEfSe MARKERS
# =========================================================

extract_lefse_markers <- function(mm) {
  
  mt <- marker_table_to_df(
    mm
  )
  
  if (nrow(mt) == 0) {
    
    return(
      tibble::tibble()
    )
  }
  
  
  # -------------------------------------------------------
  # Current microbiomeMarker source uses ef_lda.
  # Keep compatibility with versions returning lda.
  # -------------------------------------------------------
  
  lda_col <- intersect(
    c(
      "ef_lda",
      "lda"
    ),
    colnames(mt)
  )
  
  
  if (length(lda_col) == 0) {
    
    stop(
      "Could not find LEfSe LDA column. Available columns: ",
      paste(
        colnames(mt),
        collapse = ", "
      )
    )
  }
  
  
  lda_col <- lda_col[1]
  
  
  # -------------------------------------------------------
  # p-value column
  # -------------------------------------------------------
  
  p_col <- intersect(
    c(
      "pvalue",
      "p_value",
      "p.value"
    ),
    colnames(mt)
  )
  
  
  if (length(p_col) == 0) {
    
    stop(
      "Could not find LEfSe p-value column. Available columns: ",
      paste(
        colnames(mt),
        collapse = ", "
      )
    )
  }
  
  
  p_col <- p_col[1]
  
  
  # -------------------------------------------------------
  # Extract vectors BEFORE using dplyr
  # -------------------------------------------------------
  
  feature_vec <- as.character(
    mt[["feature"]]
  )
  
  enriched_vec <- as.character(
    mt[["enrich_group"]]
  )
  
  lda_vec <- as.numeric(
    mt[[lda_col]]
  )
  
  p_vec <- as.numeric(
    mt[[p_col]]
  )
  
  
  # -------------------------------------------------------
  # Build completely new tibble
  # -------------------------------------------------------
  
  out <- tibble::tibble(
    
    Taxon = feature_vec,
    
    Enriched = enriched_vec,
    
    LDA = abs(
      lda_vec
    ),
    
    Pvalue = p_vec
  )
  
  
  # -------------------------------------------------------
  # Clean names + signed LDA
  #
  # Negative = Harvest
  # Positive = Storage
  # -------------------------------------------------------
  
  out <- out %>%
    
    dplyr::mutate(
      
      # Keep terminal taxonomy component if necessary
      Taxon = sub(
        "^.*\\|",
        "",
        Taxon
      ),
      
      # Remove genus/rank prefix if present
      Taxon = sub(
        "^[A-Za-z]__",
        "",
        Taxon
      ),
      
      # Previous name cleaning
      Taxon = gsub(
        "_gen_Incertae_sedis",
        "",
        Taxon,
        fixed = TRUE
      ),
      
      SignedLDA = dplyr::case_when(
        
        Enriched == "Harvest" ~
          -LDA,
        
        Enriched == "Storage" ~
          LDA,
        
        TRUE ~
          NA_real_
      )
    ) %>%
    
    dplyr::select(
      Taxon,
      Enriched,
      LDA,
      SignedLDA,
      Pvalue
    )
  
  
  return(out)
}


# =========================================================
# STORAGE OBJECTS
# =========================================================

lefse_long_list <- list()

lefse_models <- list()

lefse_diagnostics_list <- list()


# =========================================================
# RUN SIX LEfSe COMPARISONS
#
# Pfatten/Vadena:
#   Control              Harvest vs Storage
#   Fludioxonil          Harvest vs Storage
#   Acidic Clays         Harvest vs Storage
#
# Sinich/Sinigo:
#   Control                  Harvest vs Storage
#   Captan + Fludioxonil     Harvest vs Storage
#   Captan + Acidic Clays    Harvest vs Storage
# =========================================================

set.seed(423542)


for (f in fields_order) {
  
  for (trt in treatments_order) {
    
    
    # -----------------------------------------------------
    # Unique result key
    # -----------------------------------------------------
    
    key <- paste(
      f,
      trt,
      sep = "__"
    )
    
    
    trt_label <- unname(
      treatment_labels[
        paste(
          f,
          trt,
          sep = "__"
        )
      ]
    )
    
    
    message("")
    message(
      "============================================================"
    )
    
    message(
      "LEfSe: ",
      f,
      " | ",
      trt_label,
      " | Harvest vs Storage"
    )
    
    message(
      "============================================================"
    )
    
    
    # =====================================================
    # 1. SELECT FIELD x TREATMENT SAMPLES
    # =====================================================
    
    meta_now <- as.data.frame(
      phyloseq::sample_data(
        ps_lefse
      )
    )
    
    
    keep_samples <- (
      as.character(meta_now$field) == f &
        as.character(meta_now$treatment) == trt &
        !is.na(meta_now$sampling)
    )
    
    
    keep_samples[
      is.na(keep_samples)
    ] <- FALSE
    
    
    sample_ids <- rownames(
      meta_now
    )[
      keep_samples
    ]
    
    
    if (length(sample_ids) < 2) {
      
      warning(
        "Skipping ",
        key,
        ": too few samples."
      )
      
      next
    }
    
    
    ps_sub <- phyloseq::prune_samples(
      sample_ids,
      ps_lefse
    )
    
    
    # Remove taxa completely absent from this subset
    ps_sub <- phyloseq::prune_taxa(
      phyloseq::taxa_sums(ps_sub) > 0,
      ps_sub
    )
    
    
    # =====================================================
    # 2. SAMPLING FACTOR
    # =====================================================
    
    sd_sub <- as.data.frame(
      phyloseq::sample_data(
        ps_sub
      )
    )
    
    
    sd_sub$sampling <- droplevels(
      factor(
        sd_sub$sampling,
        levels = c(
          "Harvest",
          "Storage"
        )
      )
    )
    
    
    phyloseq::sample_data(ps_sub) <-
      phyloseq::sample_data(
        sd_sub
      )
    
    
    message(
      "Samples:"
    )
    
    
    print(
      table(
        sd_sub$sampling
      )
    )
    
    
    # Must contain exactly Harvest + Storage
    if (
      nlevels(sd_sub$sampling) != 2
    ) {
      
      warning(
        "Skipping ",
        key,
        ": Harvest and Storage are not both represented."
      )
      
      next
    }
    
    
    # =====================================================
    # 3. AGGLOMERATE RAW COUNTS TO GENUS
    # =====================================================
    
    ps_genus <- phyloseq::tax_glom(
      
      ps_sub,
      
      taxrank = "Genus",
      
      NArm = TRUE
    )
    
    
    # Remove zero-total genera
    ps_genus <- phyloseq::prune_taxa(
      
      phyloseq::taxa_sums(
        ps_genus
      ) > 0,
      
      ps_genus
    )
    
    
    if (
      phyloseq::ntaxa(ps_genus) == 0
    ) {
      
      warning(
        "Skipping ",
        key,
        ": no genera after taxonomic aggregation."
      )
      
      next
    }
    
    
    # =====================================================
    # 4. RENAME TAXA TO ACTUAL GENUS NAMES
    #
    # This allows taxa_rank = "none" below.
    # run_lefse therefore analyses this already
    # genus-agglomerated table without aggregating again.
    # =====================================================
    
    tax_genus <- as.data.frame(
      phyloseq::tax_table(
        ps_genus
      )
    )
    
    
    genus_names <- as.character(
      tax_genus$Genus
    )
    
    
    valid_genus <- (
      !is.na(genus_names) &
        trimws(genus_names) != ""
    )
    
    
    ps_genus <- phyloseq::prune_taxa(
      valid_genus,
      ps_genus
    )
    
    
    tax_genus <- as.data.frame(
      phyloseq::tax_table(
        ps_genus
      )
    )
    
    
    genus_names <- as.character(
      tax_genus$Genus
    )
    
    
    # tax_glom should result in unique genera.
    # Stop if that assumption is violated.
    if (
      anyDuplicated(genus_names) > 0
    ) {
      
      stop(
        "Duplicate genus labels remain after tax_glom in ",
        key
      )
    }
    
    
    phyloseq::taxa_names(
      ps_genus
    ) <- genus_names
    
    
    message(
      "Genera before invariant filtering: ",
      phyloseq::ntaxa(
        ps_genus
      )
    )
    
    
    # =====================================================
    # 5. EXTRACT RAW GENUS COUNTS
    # =====================================================
    
    counts_genus <- as(
      phyloseq::otu_table(
        ps_genus
      ),
      "matrix"
    )
    
    
    if (
      !phyloseq::taxa_are_rows(
        ps_genus
      )
    ) {
      
      counts_genus <- t(
        counts_genus
      )
    }
    
    
    # =====================================================
    # 6. REMOVE INVARIANT GENERA
    #
    # ONLY filtering step beyond removal of zero-total taxa.
    #
    # NO prevalence filter
    # NO minimum abundance filter
    # =====================================================
    
    keep_var <- apply(
      
      counts_genus,
      
      1,
      
      function(x) {
        
        length(
          unique(x)
        ) > 1
      }
    )
    
    
    keep_var[
      is.na(keep_var)
    ] <- FALSE
    
    
    n_invariant <- sum(
      !keep_var
    )
    
    
    message(
      "Invariant genera removed: ",
      n_invariant
    )
    
    
    ps_genus <- phyloseq::prune_taxa(
      
      keep_var,
      
      ps_genus
    )
    
    
    if (
      phyloseq::ntaxa(
        ps_genus
      ) == 0
    ) {
      
      warning(
        "Skipping ",
        key,
        ": all genera were invariant."
      )
      
      next
    }
    
    
    message(
      "Genera entering LEfSe: ",
      phyloseq::ntaxa(
        ps_genus
      )
    )
    
    
    # =====================================================
    # 7. RUN LEfSe
    #
    # group = sampling
    #
    # Therefore the ONLY comparison here is:
    #
    #       Harvest vs Storage
    #
    # within the current Field x Treatment.
    #
    # Raw counts enter the function.
    #
    # CSS normalization occurs INSIDE run_lefse().
    #
    # taxa_rank = "none":
    # object is already genus-agglomerated.
    # =====================================================
    
    res_lefse <- tryCatch(
      
      microbiomeMarker::run_lefse(
        
        ps =
          ps_genus,
        
        group =
          "sampling",
        
        subgroup =
          NULL,
        
        taxa_rank =
          "none",
        
        transform =
          "identity",
        
        # -----------------------------------------------
        # CSS NORMALIZATION INSIDE LEfSe
        # -----------------------------------------------
        
        norm =
          "CSS",
        
        norm_para =
          list(),
        
        # -----------------------------------------------
        # LEfSe thresholds
        # -----------------------------------------------
        
        kw_cutoff =
          lefse_kw_cutoff,
        
        wilcoxon_cutoff =
          lefse_wilcoxon_cutoff,
        
        lda_cutoff =
          lefse_lda_cutoff,
        
        bootstrap_n =
          lefse_bootstrap_n,
        
        bootstrap_fraction =
          lefse_bootstrap_fraction,
        
        sample_min =
          lefse_sample_min,
        
        # Two groups only
        multigrp_strat =
          FALSE,
        
        strict =
          "0",
        
        only_same_subgrp =
          FALSE,
        
        curv =
          FALSE
      ),
      
      
      error = function(e) {
        
        warning(
          "LEfSe failed for ",
          key,
          ": ",
          conditionMessage(e)
        )
        
        NULL
      }
    )
    
    
    if (
      is.null(res_lefse)
    ) {
      
      next
    }
    
    
    # =====================================================
    # 8. SAVE FULL microbiomeMarker OBJECT
    # =====================================================
    
    lefse_models[[key]] <- res_lefse
    
    
    # =====================================================
    # 9. EXTRACT SIGNIFICANT MARKERS
    # =====================================================
    
    marker_i <- extract_lefse_markers(
      res_lefse
    )
    
    
    message(
      "Significant LEfSe markers: ",
      nrow(marker_i)
    )
    
    
    if (
      nrow(marker_i) == 0
    ) {
      
      next
    }
    
    
    # =====================================================
    # 10. SAFETY CHECK
    #
    # Enriched group MUST ONLY be Harvest / Storage.
    # =====================================================
    
    unexpected_groups <- setdiff(
      
      unique(
        marker_i$Enriched
      ),
      
      c(
        "Harvest",
        "Storage"
      )
    )
    
    
    if (
      length(unexpected_groups) > 0
    ) {
      
      stop(
        
        "Unexpected LEfSe enriched groups in ",
        key,
        ": ",
        paste(
          unexpected_groups,
          collapse = ", "
        )
      )
    }
    
    
    # =====================================================
    # 11. ADD FIELD / TREATMENT INFORMATION
    # =====================================================
    
    marker_i <- marker_i %>%
      
      dplyr::mutate(
        
        Field =
          f,
        
        Treatment =
          trt,
        
        TreatmentLabel =
          trt_label,
        
        .before =
          1
      )
    
    
    lefse_long_list[[key]] <- marker_i
    
    
    # =====================================================
    # 12. DIAGNOSTIC TABLE
    # =====================================================
    
    lefse_diagnostics_list[[key]] <- tibble::tibble(
      
      Field =
        f,
      
      Treatment =
        trt,
      
      TreatmentLabel =
        trt_label,
      
      Harvest_n =
        sum(
          sd_sub$sampling ==
            "Harvest"
        ),
      
      Storage_n =
        sum(
          sd_sub$sampling ==
            "Storage"
        ),
      
      GeneraBeforeInvariant =
        phyloseq::ntaxa(
          ps_genus
        ) +
        n_invariant,
      
      InvariantRemoved =
        n_invariant,
      
      GeneraEnteringLEfSe =
        phyloseq::ntaxa(
          ps_genus
        ),
      
      SignificantMarkers =
        nrow(
          marker_i
        )
    )
  }
}


# =========================================================
# COMBINE RESULTS FROM SIX LEfSe RUNS
# =========================================================

lefse_long <- dplyr::bind_rows(
  lefse_long_list
)


lefse_diagnostics <- dplyr::bind_rows(
  lefse_diagnostics_list
)


# =========================================================
# CHECK THE ALGORITHM
#
# This MUST return only:
# Harvest
# Storage
# =========================================================

print(
  unique(
    lefse_long$Enriched
  )
)


if (
  !all(
    unique(
      lefse_long$Enriched
    ) %in%
    c(
      "Harvest",
      "Storage"
    )
  )
) {
  
  stop(
    "LEfSe result contains groups other than Harvest/Storage."
  )
}


# =========================================================
# CHECK ALL SIX FIELD x TREATMENT COMPARISONS
# =========================================================

print(
  lefse_diagnostics
)


print(
  lefse_long %>%
    
    dplyr::count(
      Field,
      Treatment,
      Enriched
    )
)


# =========================================================
# COMPLETE RESULTS
# =========================================================

lefse_long %>%
  
  dplyr::arrange(
    
    Field,
    
    Treatment,
    
    dplyr::desc(
      LDA
    )
  ) %>%
  
  print(
    n = Inf
  )


# =========================================================
# EXPORT RESULTS
# =========================================================

openxlsx::write.xlsx(
  
  list(
    
    "LEfSe_markers" =
      lefse_long,
    
    "Diagnostics" =
      lefse_diagnostics
  ),
  
  file =
    "LEfSe_CSS_Harvest_vs_Storage.xlsx",
  
  overwrite =
    TRUE
)


# =========================================================
# FIGURE:
# KEEP STRONGEST TAXA PER FIELD SEPARATELY
#
# This follows your ORIGINAL bubble-plot logic.
# =========================================================

# =========================================================
# LDA FILTER FOR FIGURE
# =========================================================

lda_plot_cutoff <- 2

lefse_long_plot <- lefse_long %>%
  dplyr::filter(
    !is.na(LDA),
    LDA >= lda_plot_cutoff
  )


# =========================================================
# KEEP STRONGEST TAXA PER FIELD
# =========================================================

top_taxa_field <- lefse_long_plot %>%
  dplyr::group_by(
    Field,
    Taxon
  ) %>%
  dplyr::summarise(
    maxLDA = max(
      abs(SignedLDA),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  dplyr::group_by(
    Field
  ) %>%
  dplyr::slice_max(
    order_by = maxLDA,
    n = 15,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()


lefse_top <- lefse_long_plot %>%
  dplyr::semi_join(
    top_taxa_field,
    by = c(
      "Field",
      "Taxon"
    )
  )


# =========================================================
# COMPLETE FIELD x TREATMENT x TAXON COMBINATIONS
#
# Empty cells indicate a genus was not significant
# for that treatment.
# =========================================================

# =========================================================
# COMPLETE FIELD x TREATMENT x TAXON COMBINATIONS
# =========================================================
# Explicitly create every treatment x selected-taxon
# combination within each field.

plot_grid_df <- top_taxa_field %>%
  
  dplyr::select(
    Field,
    Taxon
  ) %>%
  
  dplyr::distinct() %>%
  
  tidyr::crossing(
    Treatment = treatments_order
  ) %>%
  
  dplyr::mutate(
    
    TreatmentLabel =
      unname(
        treatment_labels[
          paste(
            Field,
            Treatment,
            sep = "__"
          )
        ]
      )
  )

# =========================================================
# JOIN LEfSe MARKERS
# =========================================================

lefse_plot_df <- plot_grid_df %>%
  
  dplyr::left_join(
    
    lefse_top,
    
    by = c(
      "Field",
      "Treatment",
      "TreatmentLabel",
      "Taxon"
    )
  ) %>%
  
  dplyr::mutate(
    
    Enriched = factor(
      Enriched,
      levels = c(
        "Harvest",
        "Storage"
      )
    ),
    
    # Only used for ordering.
    # No bubble is drawn if Enriched is NA.
    LDA = ifelse(
      is.na(LDA),
      0,
      LDA
    ),
    
    SignedLDA = ifelse(
      is.na(SignedLDA),
      0,
      SignedLDA
    ),
    
    Taxon =
      stringr::str_replace_all(
        Taxon,
        "_",
        " "
      ),
    
    TreatmentLabel =
      dplyr::case_when(
        
        Field == "Pfatten/Vadena" &
          Treatment == "Control" ~
          "Control",
        
        Field == "Pfatten/Vadena" &
          Treatment == "Geoxe" ~
          "Fludioxonil",
        
        Field == "Pfatten/Vadena" &
          Treatment == "Ulmasud" ~
          "Acidic Clays",
        
        Field == "Sinich/Sinigo" &
          Treatment == "Control" ~
          "Control",
        
        Field == "Sinich/Sinigo" &
          Treatment == "Geoxe" ~
          "Captan +\nFludioxonil",
        
        Field == "Sinich/Sinigo" &
          Treatment == "Ulmasud" ~
          "Captan +\nAcidic Clays",
        
        TRUE ~
          TreatmentLabel
      )
  )


# =========================================================
# SHARED TAXON ORDER ACROSS BOTH FIELDS
#
# Common genera appear on exactly the same horizontal row.
# =========================================================

global_taxon_order <- lefse_plot_df %>%
  
  dplyr::group_by(
    Taxon
  ) %>%
  
  dplyr::summarise(
    
    maxLDA =
      max(
        LDA,
        na.rm = TRUE
      ),
    
    .groups =
      "drop"
  ) %>%
  
  dplyr::arrange(
    maxLDA
  ) %>%
  
  dplyr::pull(
    Taxon
  )


lefse_plot_df <- lefse_plot_df %>%
  
  dplyr::mutate(
    
    Taxon =
      factor(
        Taxon,
        levels =
          global_taxon_order
      ),
    
    TreatmentLabel =
      factor(
        
        TreatmentLabel,
        
        levels = c(
          "Control",
          "Fludioxonil",
          "Acidic Clays",
          "Captan +\nFludioxonil",
          "Captan +\nAcidic Clays"
        )
      ),
    
    Field =
      factor(
        Field,
        levels =
          fields_order
      )
  )


# =========================================================
# COMPACT ALIGNED BUBBLE PLOT
#
# SAME VISUAL STRUCTURE AS YOUR ORIGINAL FIGURE
#
# x    = treatment
# y    = genus
# fill = Harvest / Storage
# size = |LDA|
# =========================================================

p_lefse_bubble_aligned <- ggplot(
  
  lefse_plot_df,
  
  aes(
    x =
      TreatmentLabel,
    
    y =
      Taxon
  )
  
) +
  
  geom_point(
    
    data =
      lefse_plot_df %>%
      
      dplyr::filter(
        !is.na(
          Enriched
        )
      ),
    
    aes(
      
      # Same scaling used in your original figure
      size =
        LDA^2,
      
      fill =
        Enriched
    ),
    
    shape =
      21,
    
    color =
      "black",
    
    stroke =
      0.25,
    
    alpha =
      0.95
  ) +
  
  facet_grid(
    
    . ~ Field,
    
    scales =
      "free_x",
    
    space =
      "free_x"
  ) +
  
  scale_y_discrete(
    drop =
      FALSE
  ) +
  
  scale_fill_manual(
    values = c(
      Harvest = "#E64B35",
      Storage = "#4DBBD5"
    ),
    name = "Enriched in"
  ) +
  ggplot2::scale_size_continuous(
    name = "|LDA|",
    
    range = c(
      1.5,
      15
    ),
    
    # Because the plotted aesthetic is LDA^2:
    # 1^2, 2^2, 3^2, 4^2, 5^2
    breaks = c(
      1,
      4,
      9,
      16,
      25
    ),
    
    labels = c(
      "1",
      "2",
      "3",
      "4",
      "5"
    ),
    
    # Force full 1-5 reference range
    limits = c(
      1,
      25
    )
  ) +
  
  theme_nature +
  
  theme(
    
    strip.text =
      element_blank(),
    
    strip.background =
      element_blank(),
    
    axis.text.x =
      element_text(
        size = 12,
        angle = 0,
        hjust = 0.5
      ),
    
    axis.text.y =
      element_text(
        size = 12
      ),
    
    panel.grid.major =
      element_blank(),
    
    panel.grid.minor =
      element_blank(),
    
    legend.position =
      "bottom",
    
    legend.box =
      "horizontal"
  ) + labs(
    x = NULL,
    y = NULL
  )


p_lefse_bubble_aligned


# =========================================================
# SAVE LEfSe PANEL
# =========================================================

saveRDS(
  
  p_lefse_bubble_aligned,
  
  "LEfSe_CSS_faceted_bubble_plot.RDS"
)


ggsave(
  
  filename =
    "LEfSe_CSS_faceted_bubble_plot.png",
  
  plot =
    p_lefse_bubble_aligned,
  
  width =
    7,
  
  height =
    5,
  
  units =
    "in",
  
  dpi =
    600
)


# =========================================================
# FINAL FIGURE 2
#
# a = beta turnover
# b = LEfSe Harvest vs Storage
# =========================================================

p_beta_final <- p_beta_final &
  
  theme(
    
    strip.text =
      element_text(
        face = "bold",
        size = 16
      ),
    
    strip.background =
      element_blank(),
    
    strip.placement =
      "outside"
  )


p_final_combined <-
  
  patchwork::wrap_elements(
    p_beta_final
  ) +
  
  patchwork::wrap_elements(
    p_lefse_bubble_aligned
  ) +
  
  patchwork::plot_layout(
    
    ncol =
      1,
    
    heights = c(
      0.32,
      0.68
    )
  ) +
  
  patchwork::plot_annotation(
    
    tag_levels =
      list(
        c(
          "a",
          "b"
        )
      )
  ) &
  
  theme(
    
    plot.tag =
      element_text(
        size = 15,
        face = "bold"
      ),
    
    plot.tag.position =
      c(
        0,
        1
      )
  )


p_final_combined


# =========================================================
# SAVE FINAL FIGURE
# =========================================================

ggsave(
  
  filename =
    "final_compact_manuscript_figure_LEfSe_CSS.png",
  
  plot =
    p_final_combined,
  
  width =
    7.0,
  
  height =
    6.2,
  
  units =
    "in",
  
  dpi =
    600
)

# sooty blotch and sooty moulds -------------------------------------------
tse_ <- tse__fresh
tse_ <- agglomerateByRank(tse_, rank = "Genus")
tse_ <- transformAssay(tse_, method = "relabundance")
tse_ <- scale_assay(tse_, format = "percent")

taxa_list <- unique(trimws(read_excel("fungal_genera_taxonomy.xlsx")$Genus))
taxa_list <- gsub("-like$", "", taxa_list)
rowData(tse_)$Genus <- trimws(rowData(tse_)$Genus)
rowData(tse_)$Genus <- gsub("-like$", "", rowData(tse_)$Genus)

genera_df <- data.frame(Genus = taxa_list)
n_cols <- 4

# Split taxa list into multiple columns
genera_matrix <- matrix(taxa_list, 
                        ncol = n_cols, 
                        byrow = TRUE)

# Convert to data frame for gt
genera_df <- as.data.frame(genera_matrix, stringsAsFactors = FALSE)
colnames(genera_df) <- paste0("Genus_", 1:n_cols)

# Create gt table
# Create gt table without column labels
gt_fungal_db <- genera_df %>%
  gt() %>%
  tab_header(
    title = "Genera of Sooty Moulds, Sooty Blotch and Flyspeck"
  ) %>%
  fmt_missing(
    columns = everything(),
    missing_text = "-"
  ) %>%
  opt_row_striping() %>%        # optional, for readability
  tab_options(
    column_labels.hidden = TRUE  # Hides the column headers
  )
gt_fungal_db
gtsave(gt_fungal_db, "fungi.png")

taxa_in_tse <- rowData(tse_)$Genus %in% taxa_list
rel_ab_subset <- assay(tse_[taxa_in_tse, ], assay = "percent", withDimnames = TRUE)
taxon_names <- rowData(tse_)[taxa_in_tse, "Genus"]

df_long <- as.data.frame(rel_ab_subset) %>%
  rownames_to_column(var = "OTU_ID") %>%
  mutate(Taxon = taxon_names) %>%
  pivot_longer(
    cols = -c(OTU_ID, Taxon),
    names_to = "Sample",
    values_to = "RelAbundance"
  ) %>%
  mutate(sampling = colData(tse_)[Sample, "sampling"])

taxa_category_df <- read_excel("fungal_genera_taxonomy.xlsx")
taxa_category_df <- read_excel("fungal_genera_taxonomy.xlsx") %>%
  mutate(
    Genus = trimws(Genus),
    Genus = gsub("-like$", "", Genus)
  )
# Add Category to df_plot


df_long_filtered <- df_long %>%
  group_by(Taxon) %>%
  filter(sum(RelAbundance > 0) > 0 & max(RelAbundance) > 0.00) %>%
  ungroup()

df_plot <- df_long_filtered %>%
  mutate(Genus = Taxon, Group = sampling) %>%
  group_by(Genus, Group) %>%
  summarise(
    mean_perc = mean(RelAbundance, na.rm = TRUE),
    samples_found = sum(RelAbundance > 0),
    .groups = "drop"
  )

df_plot$Group <- factor(df_plot$Group, levels = c("Storage", "Harvest"))

df_plot  <- df_plot %>%
  dplyr::left_join(taxa_category_df, by = "Genus")

n_genera <- length(unique(df_plot$Genus))
n_genera

# Ensure Group is a factor
colnames(df_plot)[colnames(df_plot) == "Group"] <- "Sampling"
df_plot$Sampling <- factor(df_plot$Sampling, levels = c("Harvest", "Storage"))

my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(unique(df_plot$Genus)))
library(tidytext)
df_plot <- df_plot %>%
  mutate(Genus = reorder_within(Genus, mean_perc, Category))  # reorder per facet

# Dumbbell plot with ordered facets
dumbell <- ggplot(df_plot, aes(y = Genus, x = mean_perc, group = Genus)) +
  geom_line(aes(color = "line"), color = "gray70", size = 1) +         
  geom_point(aes(color = Sampling), size = 4) +                           
  geom_text_repel(aes(label = samples_found),
                  size = 6, nudge_x = 0.01, show.legend = FALSE) +     
  facet_wrap(~Category, scales = "free_x") +
  scale_y_reordered() +  # ensures proper order per facet
  scale_color_manual(values = c("Harvest" ="#F8766D",
                                "Storage" =  "#00BFC4", 
                                
                                "line" = "gray70")) +
  coord_flip() +
  theme_minimal(base_size = 16) +
  labs(x = "Mean %", y = NULL, color = "Sampling") +
  theme_nature + 
  theme(
    axis.text.y = element_text(face = "bold", size = 12),
    axis.text.x = element_text(size = 17, angle = 45, hjust = 1),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 16) 
  )
dumbell
saveRDS(dumbell, "dumbell.RDS")


# =========================================================
# STEP 8: COMMUNITY-ASSEMBLY NULL MODELLING WITH iCAMP
# =========================================================
# Reviewer L318/L389:
# - DOC is no longer used as evidence for deterministic/stochastic assembly;
# - iCAMP partitions assembly into heterogeneous selection (HeS),
#   homogeneous selection (HoS), dispersal limitation (DL),
#   homogenizing dispersal (HD), and drift (DR);
# - bin-size diagnostics are evaluated explicitly;
# - bins dominated by selection are mapped back to their ASVs/taxonomy.
#
# iCAMP documentation recommends exploring bin.size.limit values in roughly
# the 12-48 range for real datasets. A formal ps.bin phylogenetic-signal test
# additionally requires a meaningful numeric environmental/niche matrix.
# Categorical Field/Treatment/Sampling labels are NOT converted to arbitrary
# numbers for this purpose.
# =========================================================

prepare_icamp_data <- function(ps) {
  comm <- get_sample_by_taxa_matrix(ps)
  tree <- phyloseq::phy_tree(ps)

  common_taxa <- intersect(colnames(comm), tree$tip.label)
  if (length(common_taxa) < 2L) {
    stop("Too few taxa overlap between the community table and phylogeny.")
  }

  comm <- comm[, common_taxa, drop = FALSE]
  tree <- ape::keep.tip(tree, common_taxa)

  keep_taxa <- colSums(comm) > 0
  comm <- comm[, keep_taxa, drop = FALSE]
  tree <- ape::keep.tip(tree, colnames(comm))

  meta <- data.frame(phyloseq::sample_data(ps))
  meta <- meta[rownames(comm), , drop = FALSE]

  taxonomy <- as.data.frame(phyloseq::tax_table(ps))
  taxonomy <- taxonomy[colnames(comm), , drop = FALSE]

  list(
    comm = comm,
    tree = tree,
    metadata = meta,
    taxonomy = taxonomy
  )
}

run_icamp_bin_diagnostics <- function(
    ps,
    prefix = "ITS",
    bin_sizes = c(12L, 24L, 48L),
    ds = 0.2,
    nworker = 4L,
    env_numeric = NULL,
    abundance_cutoff = 3
) {
  x <- prepare_icamp_data(ps)
  comm <- x$comm
  tree <- x$tree

  outdir <- file.path(getwd(), paste0("iCAMP_", prefix, "_bin_diagnostics"))
  pd_dir <- file.path(outdir, "phylogenetic_distance")

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(pd_dir, showWarnings = FALSE, recursive = TRUE)

  # Compute the large phylogenetic-distance backing files once.
  pd_big <- iCAMP::pdist.big(
    tree = tree,
    wd = pd_dir,
    nworker = nworker
  )

  structural_list <- list()
  signal_list <- list()

  niche_diff <- NULL

  if (!is.null(env_numeric)) {
    if (is.null(rownames(env_numeric))) {
      stop("env_numeric must have sample IDs as row names.")
    }

    if (!all(rownames(comm) %in% rownames(env_numeric))) {
      stop("env_numeric does not contain all sample IDs in the community matrix.")
    }

    env_numeric <- env_numeric[rownames(comm), , drop = FALSE]

    if (!all(vapply(env_numeric, is.numeric, logical(1)))) {
      stop(
        "env_numeric must contain meaningful numeric environmental variables. ",
        "Do not recode Field/Treatment/Sampling factors to arbitrary integers."
      )
    }

    niche_dir <- file.path(outdir, "niche_distance")
    dir.create(niche_dir, showWarnings = FALSE, recursive = TRUE)

    niche_diff <- iCAMP::dniche(
      env = env_numeric,
      comm = comm,
      method = "niche.value",
      nworker = nworker,
      out.dist = FALSE,
      bigmemo = TRUE,
      nd.wd = niche_dir,
      nd.spname.file = paste0(prefix, "_nd.names.csv")
    )
  }

  for (bin_size in bin_sizes) {
    message("\nTesting iCAMP bin.size.limit = ", bin_size)

    phylobin <- iCAMP::taxa.binphy.big(
      tree = tree,
      pd.desc = pd_big$pd.file,
      pd.spname = pd_big$tip.label,
      pd.wd = pd_big$pd.wd,
      ds = ds,
      bin.size.limit = bin_size,
      nworker = nworker
    )

    final_bin <- phylobin$sp.bin[, 3]
    bin_counts <- table(final_bin)

    structural_list[[as.character(bin_size)]] <- tibble::tibble(
      BinSizeLimit = bin_size,
      N_bins = length(bin_counts),
      Minimum_actual_bin_size = min(bin_counts),
      Median_actual_bin_size = median(bin_counts),
      Mean_actual_bin_size = mean(bin_counts),
      Maximum_actual_bin_size = max(bin_counts)
    )

    write.csv(
      as.data.frame(phylobin$state.united),
      file.path(
        outdir,
        paste0(prefix, "_bin_", bin_size, "_state_united.csv")
      ),
      row.names = FALSE
    )

    # Formal within-bin phylogenetic-signal test, only if actual numeric
    # environmental/niche variables are supplied.
    if (!is.null(niche_diff)) {
      sp_bin <- phylobin$sp.bin[, 3, drop = FALSE]
      sp_ra <- colMeans(comm / rowSums(comm))
      spname_use <- colnames(comm)[colSums(comm) >= abundance_cutoff]

      ps_test <- iCAMP::ps.bin(
        sp.bin = sp_bin,
        sp.ra = sp_ra,
        spname.use = spname_use,
        pd.desc = pd_big$pd.file,
        pd.spname = pd_big$tip.label,
        pd.wd = pd_big$pd.wd,
        nd.list = niche_diff$nd,
        nd.spname = niche_diff$names,
        ndbig.wd = niche_diff$nd.wd,
        cor.method = "pearson",
        r.cut = 0.01,
        p.cut = 0.2,
        min.spn = 6
      )

      signal_index <- as.data.frame(ps_test$Index) %>%
        tibble::rownames_to_column("IndexRow") %>%
        mutate(BinSizeLimit = bin_size, .before = 1)

      signal_list[[as.character(bin_size)]] <- signal_index

      if (!is.null(ps_test$detail)) {
        write.csv(
          as.data.frame(ps_test$detail),
          file.path(
            outdir,
            paste0(prefix, "_bin_", bin_size, "_phylogenetic_signal_detail.csv")
          ),
          row.names = TRUE
        )
      }
    }
  }

  structural <- bind_rows(structural_list)
  signal <- bind_rows(signal_list)

  write.csv(
    structural,
    file.path(outdir, paste0(prefix, "_bin_size_diagnostics.csv")),
    row.names = FALSE
  )

  if (nrow(signal) > 0) {
    write.csv(
      signal,
      file.path(outdir, paste0(prefix, "_bin_phylogenetic_signal_summary.csv")),
      row.names = FALSE
    )
  } else {
    message(
      "\nFormal ps.bin was not run because env_numeric = NULL. ",
      "This is intentional: categorical experimental labels are not valid ",
      "numeric niche distances. Structural bin-size diagnostics were exported."
    )
  }

  print(structural)
  if (nrow(signal) > 0) print(signal)

  list(
    structural = structural,
    phylogenetic_signal = signal,
    pd_big = pd_big,
    pd_dir = pd_dir,
    output_dir = outdir
  )
}

run_icamp_reviewer <- function(
    ps,
    prefix = "ITS",
    bin_size = 24L,
    rand = 1000L,
    nworker = 4L,
    ds = 0.2,
    pd_big = NULL,
    pd_dir = NULL
) {
  x <- prepare_icamp_data(ps)
  comm <- x$comm
  tree <- x$tree
  meta <- x$metadata
  taxonomy <- x$taxonomy

  outdir <- file.path(getwd(), paste0("iCAMP_results_", prefix))
  if (is.null(pd_dir)) {
    pd_dir <- file.path(outdir, "phylogenetic_distance")
  }

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  dir.create(pd_dir, showWarnings = FALSE, recursive = TRUE)

  if (is.null(pd_big)) {
    pd_big <- iCAMP::pdist.big(
      tree = tree,
      wd = pd_dir,
      nworker = nworker
    )
  }

  groups <- data.frame(
    Sampling = as.character(meta$sampling),
    Field = as.character(meta$field),
    Treatment = as.character(meta$treatment),
    Field_Treatment = interaction(meta$field, meta$treatment, drop = TRUE),
    Field_Sampling = interaction(meta$field, meta$sampling, drop = TRUE),
    Field_Treatment_Sampling = interaction(
      meta$field,
      meta$treatment,
      meta$sampling,
      drop = TRUE
    ),
    stringsAsFactors = FALSE,
    row.names = rownames(meta)
  )

  set.seed(423542)

  icamp_out <- iCAMP::icamp.big(
    comm = comm,
    tree = tree,
    pd.desc = pd_big$pd.file,
    pd.spname = pd_big$tip.label,
    pd.wd = pd_big$pd.wd,
    output.wd = outdir,
    rand = rand,
    prefix = paste0(prefix, "_apple_sooty_epiphytes"),
    ds = ds,
    phylo.rand.scale = "within.bin",
    taxa.rand.scale = "across.all",
    phylo.metric = "bMPD",
    sig.index = "SES.RC",
    bin.size.limit = bin_size,
    nworker = nworker,
    detail.save = TRUE,
    qp.save = TRUE,
    detail.null = FALSE
  )

  saveRDS(
    icamp_out,
    file.path(outdir, paste0(prefix, "_icamp_full_output.rds"))
  )

  icamp_bins <- iCAMP::icamp.bins(
    icamp.detail = icamp_out,
    treat = groups,
    clas = taxonomy,
    boot = TRUE,
    rand.time = rand,
    between.group = TRUE
  )

  saveRDS(
    icamp_bins,
    file.path(outdir, paste0(prefix, "_icamp_bin_summary.rds"))
  )

  # Export all table-like components for Supplementary material.
  export_icamp_component <- function(x, name) {
    if (is.data.frame(x) || is.matrix(x)) {
      write.csv(
        as.data.frame(x),
        file.path(outdir, paste0(name, ".csv")),
        row.names = TRUE
      )
    } else if (is.list(x)) {
      for (i in seq_along(x)) {
        xi <- x[[i]]
        nm <- names(x)[i]
        if (is.null(nm) || nm == "") nm <- paste0("part", i)

        if (is.data.frame(xi) || is.matrix(xi)) {
          write.csv(
            as.data.frame(xi),
            file.path(
              outdir,
              paste0(name, "_", make.names(nm), ".csv")
            ),
            row.names = TRUE
          )
        }
      }
    }
  }

  for (nm in names(icamp_bins)) {
    export_icamp_component(icamp_bins[[nm]], paste0(prefix, "_iCAMP_", nm))
  }

  list(
    output = icamp_out,
    bins = icamp_bins,
    groups = groups,
    taxonomy = taxonomy,
    output_directory = outdir
  )
}

extract_icamp_selected_taxa <- function(icamp_run, prefix = "ITS") {
  ptk_raw <- icamp_run$bins$Ptk
  class_bin <- as.data.frame(
    icamp_run$bins$Class.Bin,
    check.names = FALSE
  )

  if (is.null(ptk_raw)) {
    stop("icamp_bins$Ptk is absent; cannot identify selection-dominated bins.")
  }

  # Ptk can differ slightly across iCAMP versions/settings. Convert either
  # one table or a list of tables to a named list and process each safely.
  ptk_list <- if (is.data.frame(ptk_raw) || is.matrix(ptk_raw)) {
    list(Ptk = as.data.frame(ptk_raw, check.names = FALSE))
  } else if (is.list(ptk_raw)) {
    tmp <- lapply(
      ptk_raw,
      function(z) {
        if (is.data.frame(z) || is.matrix(z)) {
          as.data.frame(z, check.names = FALSE)
        } else {
          NULL
        }
      }
    )
    tmp[!vapply(tmp, is.null, logical(1))]
  } else {
    stop("Unsupported Ptk object class: ", paste(class(ptk_raw), collapse = ", "))
  }

  find_one <- function(pattern, x) {
    hit <- grep(pattern, x, ignore.case = TRUE, value = TRUE)
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  class_bin_col <- find_one("bin", colnames(class_bin))
  if (is.na(class_bin_col)) {
    stop(
      "Could not identify the bin-ID column in Class.Bin. Columns: ",
      paste(colnames(class_bin), collapse = ", ")
    )
  }

  class_bin$BinID_JOIN <- as.character(class_bin[[class_bin_col]])

  selected_bins_all <- list()
  selected_taxa_all <- list()

  for (nm in names(ptk_list)) {
    ptk <- ptk_list[[nm]]

    HoS_col <- find_one("(^HoS$|homogeneous.*selection)", colnames(ptk))
    HeS_col <- find_one("(^HeS$|heterogeneous.*selection)", colnames(ptk))
    DL_col  <- find_one("(^DL$|dispersal.*limitation)", colnames(ptk))
    HD_col  <- find_one("(^HD$|homogenizing.*dispersal)", colnames(ptk))
    DR_col  <- find_one("(^DR$|drift)", colnames(ptk))
    bin_col <- find_one("bin", colnames(ptk))

    if (is.na(HoS_col) || is.na(HeS_col) || is.na(bin_col)) {
      warning(
        "Skipping Ptk component '", nm,
        "' because HoS/HeS/bin columns could not be detected. Columns: ",
        paste(colnames(ptk), collapse = ", ")
      )
      next
    }

    process_map <- c(
      HoS = HoS_col,
      HeS = HeS_col,
      DL = DL_col,
      HD = HD_col,
      DR = DR_col
    )
    process_map <- process_map[!is.na(process_map)]

    process_matrix <- ptk[, unname(process_map), drop = FALSE]
    process_matrix[] <- lapply(process_matrix, as.numeric)
    colnames(process_matrix) <- names(process_map)

    ptk$HomogeneousSelection <- process_matrix$HoS
    ptk$HeterogeneousSelection <- process_matrix$HeS
    ptk$TotalSelection <- process_matrix$HoS + process_matrix$HeS

    ptk$DominantProcess <- apply(
      process_matrix,
      1,
      function(z) {
        if (all(is.na(z))) return(NA_character_)
        names(z)[which.max(z)]
      }
    )

    ptk$SelectionDominated <- ptk$DominantProcess %in% c("HoS", "HeS")
    ptk$BinID_JOIN <- as.character(ptk[[bin_col]])
    ptk$PtkComponent <- nm

    selected_bins <- ptk %>%
      filter(SelectionDominated) %>%
      arrange(desc(TotalSelection))

    selected_taxa <- selected_bins %>%
      inner_join(
        class_bin,
        by = "BinID_JOIN",
        suffix = c("_Process", "_Taxonomy")
      ) %>%
      arrange(desc(TotalSelection))

    selected_bins_all[[nm]] <- selected_bins
    selected_taxa_all[[nm]] <- selected_taxa
  }

  selected_bins_df <- bind_rows(selected_bins_all)
  selected_taxa_df <- bind_rows(selected_taxa_all)

  outdir <- icamp_run$output_directory

  write.csv(
    selected_bins_df,
    file.path(outdir, paste0(prefix, "_selection_dominated_bins.csv")),
    row.names = FALSE
  )

  write.csv(
    selected_taxa_df,
    file.path(outdir, paste0(prefix, "_taxa_in_selection_dominated_bins.csv")),
    row.names = FALSE
  )

  # Genus-level reviewer table when genus taxonomy is available.
  genus_col <- grep(
    "^Genus$",
    colnames(selected_taxa_df),
    ignore.case = TRUE,
    value = TRUE
  )

  if (length(genus_col) == 1L && nrow(selected_taxa_df) > 0) {
    selected_genera <- selected_taxa_df %>%
      filter(!is.na(.data[[genus_col]]), .data[[genus_col]] != "") %>%
      group_by(
        PtkComponent,
        DominantProcess,
        .data[[genus_col]]
      ) %>%
      summarise(
        N_ASVs = n(),
        Mean_TotalSelection = mean(TotalSelection, na.rm = TRUE),
        Max_TotalSelection = max(TotalSelection, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(Max_TotalSelection))

    write.csv(
      selected_genera,
      file.path(outdir, paste0(prefix, "_selected_genera_summary.csv")),
      row.names = FALSE
    )
  }

  list(
    selected_bins = selected_bins_df,
    selected_taxa = selected_taxa_df
  )
}

# -----------------------------------------------------------------------------
# ITS iCAMP
# -----------------------------------------------------------------------------
# No meaningful continuous environmental/niche variables are currently supplied,
# therefore ps.bin is intentionally skipped. If you have measured continuous
# variables, create a numeric data.frame with Sample IDs as row names and pass it
# as env_numeric below.
ICAMP_ENV_NUMERIC_ITS <- NULL
ICAMP_BIN_CANDIDATES <- c(12L, 24L, 48L)
ICAMP_BIN_SIZE_ITS <- 24L
ICAMP_RAND <- 1000L
ICAMP_WORKERS <- 4L

ITS_bincheck <- run_icamp_bin_diagnostics(
  physeq_asv_filtered,
  prefix = "ITS",
  bin_sizes = ICAMP_BIN_CANDIDATES,
  nworker = ICAMP_WORKERS,
  env_numeric = ICAMP_ENV_NUMERIC_ITS
)

ITS_icamp <- run_icamp_reviewer(
  physeq_asv_filtered,
  prefix = "ITS",
  bin_size = ICAMP_BIN_SIZE_ITS,
  rand = ICAMP_RAND,
  nworker = ICAMP_WORKERS,
  pd_big = ITS_bincheck$pd_big,
  pd_dir = ITS_bincheck$pd_dir
)

ITS_selected_taxa <- extract_icamp_selected_taxa(
  ITS_icamp,
  prefix = "ITS"
)

print(ITS_icamp$bins$Pt)
print(head(ITS_icamp$bins$Ptk))
print(head(ITS_icamp$bins$Class.Bin))
print(head(ITS_selected_taxa$selected_taxa))

# =========================================================
# PLOT ASSEMBLY MECHANISMS AT HARVEST VS STORAGE
# =========================================================
#
# Goal:
# show the relative importance of the five assembly processes
# separately at Harvest and at Storage.
#
# Each panel = one Field x Treatment combination
# Each panel contains two stacked bars:
#   Harvest
#   Storage
# =========================================================


# ---------------------------------------------------------
# Helper: convert Pt to data.frame
# ---------------------------------------------------------
icamp_component_to_df <- function(x, component_name = "iCAMP component") {
  
  if (is.data.frame(x) || is.matrix(x)) {
    return(
      as.data.frame(
        x,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    )
  }
  
  if (is.list(x)) {
    
    keep <- vapply(
      x,
      function(z) is.data.frame(z) || is.matrix(z),
      logical(1)
    )
    
    x <- x[keep]
    
    if (length(x) == 0L) {
      stop(component_name, " contains no table-like elements.")
    }
    
    nm <- names(x)
    if (is.null(nm)) nm <- paste0("part", seq_along(x))
    nm[nm == ""] <- paste0("part", which(nm == ""))
    
    pieces <- lapply(
      seq_along(x),
      function(i) {
        z <- as.data.frame(
          x[[i]],
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        z$.iCAMP_component <- nm[i]
        z
      }
    )
    
    return(dplyr::bind_rows(pieces))
  }
  
  stop(
    "Unsupported class for ",
    component_name,
    ": ",
    paste(class(x), collapse = ", ")
  )
}


# ---------------------------------------------------------
# Helper: identify process columns
# ---------------------------------------------------------
detect_icamp_process_columns <- function(dat) {
  
  find_one <- function(pattern) {
    hit <- grep(
      pattern,
      colnames(dat),
      ignore.case = TRUE,
      value = TRUE
    )
    if (length(hit) == 0L) return(NA_character_)
    hit[1]
  }
  
  out <- c(
    HeS = find_one("(^HeS$|heterogeneous.*selection)"),
    HoS = find_one("(^HoS$|homogeneous.*selection)"),
    DL  = find_one("(^DL$|dispersal.*limitation)"),
    HD  = find_one("(^HD$|homogenizing.*dispersal)"),
    DR  = find_one("(^DR$|drift)")
  )
  
  if (any(is.na(out))) {
    stop(
      paste0(
        "Could not identify all process columns.\n",
        "Detected:\n",
        paste(names(out), "=", out, collapse = "\n"),
        "\n\nAvailable columns:\n",
        paste(colnames(dat), collapse = ", ")
      )
    )
  }
  
  out
}


# ---------------------------------------------------------
# Helper: treatment labels
#
# Works for BOTH:
#   1. stagewise tables containing Sampling
#   2. Harvest-vs-Storage transition tables without Sampling
# ---------------------------------------------------------

add_icamp_treatment_labels <- function(dat) {
  
  dat <- dat %>%
    dplyr::mutate(
      
      Field =
        as.character(.data$Field),
      
      Treatment =
        as.character(.data$Treatment),
      
      Treatment_label =
        dplyr::case_when(
          
          Field == "Pfatten/Vadena" &
            Treatment == "Control" ~
            "Control",
          
          Field == "Pfatten/Vadena" &
            Treatment == "Geoxe" ~
            "Fludioxonil",
          
          Field == "Pfatten/Vadena" &
            Treatment == "Ulmasud" ~
            "Acidic clays",
          
          Field == "Sinich/Sinigo" &
            Treatment == "Control" ~
            "Control",
          
          Field == "Sinich/Sinigo" &
            Treatment == "Geoxe" ~
            "Captan + Fludioxonil",
          
          Field == "Sinich/Sinigo" &
            Treatment == "Ulmasud" ~
            "Captan + Acidic clays",
          
          TRUE ~
            Treatment
        ),
      
      Field =
        factor(
          Field,
          levels = c(
            "Pfatten/Vadena",
            "Sinich/Sinigo"
          )
        )
    )
  
  
  # -------------------------------------------------------
  # Only standardize Sampling when the table actually
  # contains a Sampling column.
  #
  # BPtk Harvest-vs-Storage transition tables do NOT have
  # one single Sampling value, so this must be skipped.
  # -------------------------------------------------------
  
  if ("Sampling" %in% colnames(dat)) {
    
    dat <- dat %>%
      dplyr::mutate(
        
        Sampling =
          factor(
            as.character(.data$Sampling),
            levels = c(
              "Harvest",
              "Storage"
            )
          )
      )
  }
  
  
  dat
}

# ---------------------------------------------------------
# Extract WITHIN-GROUP process weights
# ---------------------------------------------------------
extract_icamp_process_weights_by_stage <- function(
    icamp_run,
    prefix = "ITS"
) {
  
  pt <- icamp_component_to_df(
    icamp_run$bins$Pt,
    component_name = "Pt"
  )
  
  # Keep Field_Treatment_Sampling summaries
  pt_stage <- pt %>%
    dplyr::mutate(
      GroupBasedOn = as.character(GroupBasedOn),
      Group = as.character(Group)
    ) %>%
    dplyr::filter(
      GroupBasedOn == "Field_Treatment_Sampling"
    )
  
  # We want SINGLE groups, not Harvest_vs_Storage transitions
  pt_stage <- pt_stage %>%
    dplyr::filter(
      !grepl("_vs_", Group)
    )
  
  # If Field/Treatment/Sampling columns are not already present, recover them from icamp_run$groups
  if (!all(c("Field", "Treatment", "Sampling") %in% colnames(pt_stage))) {
    
    lookup <- icamp_run$groups %>%
      dplyr::mutate(
        Group = as.character(Field_Treatment_Sampling),
        Field = as.character(Field),
        Treatment = as.character(Treatment),
        Sampling = as.character(Sampling)
      ) %>%
      dplyr::distinct(Group, Field, Treatment, Sampling)
    
    pt_stage <- pt_stage %>%
      dplyr::left_join(
        lookup,
        by = "Group"
      )
  }
  
  if (!all(c("Field", "Treatment", "Sampling") %in% colnames(pt_stage))) {
    stop("Could not recover Field, Treatment and Sampling information from Pt.")
  }
  
  process_map <- detect_icamp_process_columns(pt_stage)
  col_to_process <- stats::setNames(names(process_map), unname(process_map))
  
  pt_long <- pt_stage %>%
    tidyr::pivot_longer(
      cols = dplyr::all_of(unname(process_map)),
      names_to = "ProcessColumn",
      values_to = "Weight"
    ) %>%
    dplyr::mutate(
      Weight = as.numeric(Weight),
      ProcessCode = unname(col_to_process[ProcessColumn]),
      Process = dplyr::recode(
        ProcessCode,
        "HeS" = "Heterogeneous selection",
        "HoS" = "Homogeneous selection",
        "DL"  = "Dispersal limitation",
        "HD"  = "Homogenizing dispersal",
        "DR"  = "Drift"
      ),
      Process = factor(
        Process,
        levels = c(
          "Heterogeneous selection",
          "Homogeneous selection",
          "Dispersal limitation",
          "Homogenizing dispersal",
          "Drift"
        )
      )
    ) %>%
    add_icamp_treatment_labels()
  
  # useful check
  sum_check <- pt_long %>%
    dplyr::group_by(Field, Treatment_label, Sampling) %>%
    dplyr::summarise(
      Total = sum(Weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  print(sum_check)
  
  utils::write.csv(
    pt_long,
    file.path(
      icamp_run$output_directory,
      paste0(prefix, "_iCAMP_process_weights_Harvest_vs_Storage_stagewise.csv")
    ),
    row.names = FALSE
  )
  
  utils::write.csv(
    sum_check,
    file.path(
      icamp_run$output_directory,
      paste0(prefix, "_iCAMP_process_weights_stagewise_sumcheck.csv")
    ),
    row.names = FALSE
  )
  
  pt_long
}


ITS_process_stagewise <- extract_icamp_process_weights_by_stage(
  ITS_icamp,
  prefix = "ITS"
)


# ---------------------------------------------------------
# Colours
# ---------------------------------------------------------
icamp_process_colours <- c(
  "Heterogeneous selection" = "#D55E00",
  "Homogeneous selection" = "#E69F00",
  "Dispersal limitation" = "#0072B2",
  "Homogenizing dispersal" = "#56B4E9",
  "Drift" = "#999999"
)


# ---------------------------------------------------------
# Final plot: Harvest vs Storage
# ---------------------------------------------------------
ITS_process_stagewise <- ITS_process_stagewise %>%
  dplyr::mutate(
    Treatment_label = factor(
      Treatment_label,
      levels = c(
        "Control",
        "Fludioxonil",
        "Acidic clays",
        "Captan + Fludioxonil",
        "Captan + Acidic clays"
      )
    )
  )

# ---------------------------------------------------------
# Remove empty panels by plotting each orchard separately
# ---------------------------------------------------------

library(patchwork)

ITS_process_stagewise_plot <- ITS_process_stagewise %>%
  dplyr::mutate(
    Treatment_label = as.character(Treatment_label),
    Sampling = factor(Sampling, levels = c("Harvest", "Storage"))
  )

# Field-specific treatment order
pfatten_levels <- c(
  "Control",
  "Fludioxonil",
  "Acidic clays"
)

sinich_levels <- c(
  "Control",
  "Captan + Fludioxonil",
  "Captan + Acidic clays"
)

pfatten_dat <- ITS_process_stagewise_plot %>%
  dplyr::filter(Field == "Pfatten/Vadena") %>%
  dplyr::mutate(
    Treatment_label = factor(
      Treatment_label,
      levels = pfatten_levels
    )
  )

sinich_dat <- ITS_process_stagewise_plot %>%
  dplyr::filter(Field == "Sinich/Sinigo") %>%
  dplyr::mutate(
    Treatment_label = factor(
      Treatment_label,
      levels = sinich_levels
    )
  )

# ---------------------------------------------------------
# Pfatten/Vadena
# ---------------------------------------------------------

p_icamp_pfatten <- ggplot2::ggplot(
  pfatten_dat,
  ggplot2::aes(
    x = Sampling,
    y = Weight,
    fill = Process
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  ggplot2::facet_wrap(
    ~ Treatment_label,
    nrow = 1,
    scales = "free_x"
  ) +
  ggplot2::scale_fill_manual(
    values = icamp_process_colours,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Pfatten/Vadena",
    x = NULL,
    y = NULL
  ) +
  theme_nature +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1,
      size = 12
    ),
    strip.text = ggplot2::element_text(
      size = 13,
      face = "bold"
    ),
    legend.position = "none"
  )


# ---------------------------------------------------------
# Sinich/Sinigo
# ---------------------------------------------------------

p_icamp_sinich <- ggplot2::ggplot(
  sinich_dat,
  ggplot2::aes(
    x = Sampling,
    y = Weight,
    fill = Process
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.25
  ) +
  ggplot2::facet_wrap(
    ~ Treatment_label,
    nrow = 1,
    scales = "free_x"
  ) +
  ggplot2::scale_fill_manual(
    values = icamp_process_colours,
    drop = FALSE
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    title = "Sinich/Sinigo",
    x = NULL,
    y = NULL
  ) +
  theme_nature +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1,
      size = 12
    ),
    strip.text = ggplot2::element_text(
      size = 13,
      face = "bold"
    ),
    legend.position = "none"
  )


# ---------------------------------------------------------
# Shared legend
# ---------------------------------------------------------

legend_icamp <- cowplot::get_legend(
  p_icamp_pfatten +
    ggplot2::theme(
      legend.position = "bottom",
      legend.title = ggplot2::element_text(
        size = 12,
        face = "bold"
      ),
      legend.text = ggplot2::element_text(
        size = 11
      )
    )
)


# ---------------------------------------------------------
# Combine orchard panels
# NO overall title
# NO subtitle
# NO repeated y-axis title
# ---------------------------------------------------------

p_icamp_core <- cowplot::plot_grid(
  p_icamp_pfatten,
  p_icamp_sinich,
  ncol = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)


# ---------------------------------------------------------
# Add ONE shared y-axis title
# ---------------------------------------------------------

p_icamp_process_stagewise_clean <- cowplot::ggdraw() +
  
  cowplot::draw_plot(
    p_icamp_core,
    x = 0.07,
    y = 0.12,
    width = 0.93,
    height = 0.88
  ) +
  
  cowplot::draw_plot(
    legend_icamp,
    x = 0.07,
    y = 0.00,
    width = 0.93,
    height = 0.12
  ) +
  
  cowplot::draw_label(
    "Relative importance of assembly processes",
    x = 0.018,
    y = 0.56,
    angle = 90,
    fontface = "bold",
    size = 18
  )


p_icamp_process_stagewise_clean
# =========================================================
# 9B. BIN CONTRIBUTIONS TO EACH PROCESS: BPtk
# =========================================================
extract_icamp_bin_contributions_HS <- function(
    icamp_run,
    prefix = "ITS"
) {
  
  bptk <- icamp_component_to_df(
    icamp_run$bins$BPtk,
    component_name = "BPtk"
  )
  
  
  bptk_hs <- filter_icamp_HS_transition(
    bptk,
    icamp_run
  )
  
  
  required <- c(
    "Process",
    "Field",
    "Treatment"
  )
  
  
  missing <- setdiff(
    required,
    colnames(bptk_hs)
  )
  
  
  if (length(missing) > 0L) {
    
    stop(
      "BPtk Harvest-vs-Storage table is missing: ",
      paste(missing, collapse = ", ")
    )
  }
  
  
  non_bin_columns <- c(
    
    "Method",
    
    "GroupBasedOn",
    
    "Group",
    
    "Process",
    
    "Field",
    
    "Treatment",
    
    ".iCAMP_component"
  )
  
  
  candidate_bins <- setdiff(
    colnames(bptk_hs),
    non_bin_columns
  )
  
  
  # Retain columns that are genuinely numeric bin-contribution columns.
  numericish <- vapply(
    
    bptk_hs[
      candidate_bins
    ],
    
    function(z) {
      
      z_num <-
        suppressWarnings(
          as.numeric(
            as.character(
              z
            )
          )
        )
      
      any(
        !is.na(
          z_num
        )
      ) &&
        all(
          is.na(z) |
            !is.na(z_num)
        )
    },
    
    logical(1)
  )
  
  
  bin_columns <-
    candidate_bins[
      numericish
    ]
  
  
  if (length(bin_columns) == 0L) {
    
    stop(
      paste0(
        "No numeric bin columns were identified in BPtk.\n",
        "Available columns:\n",
        paste(
          colnames(bptk_hs),
          collapse = ", "
        )
      )
    )
  }
  
  
  bptk_long <- bptk_hs %>%
    tidyr::pivot_longer(
      
      cols =
        dplyr::all_of(
          bin_columns
        ),
      
      names_to =
        "Bin",
      
      values_to =
        "Contribution"
    ) %>%
    dplyr::mutate(
      
      Contribution =
        as.numeric(
          Contribution
        ),
      
      ProcessCode =
        standardize_icamp_process(
          Process
        ),
      
      Process =
        dplyr::recode(
          
          ProcessCode,
          
          "HeS" =
            "Heterogeneous selection",
          
          "HoS" =
            "Homogeneous selection",
          
          "DL" =
            "Dispersal limitation",
          
          "HD" =
            "Homogenizing dispersal",
          
          "DR" =
            "Drift"
        ),
      
      Process =
        factor(
          Process,
          levels = c(
            "Heterogeneous selection",
            "Homogeneous selection",
            "Dispersal limitation",
            "Homogenizing dispersal",
            "Drift"
          )
        )
    ) %>%
    add_icamp_treatment_labels() %>%
    dplyr::mutate(
      
      Panel =
        paste0(
          as.character(
            Field
          ),
          "\n",
          Treatment_label
        )
    )
  
  
  # Order bins numerically when possible.
  bin_number <-
    suppressWarnings(
      readr::parse_number(
        as.character(
          bptk_long$Bin
        )
      )
    )
  
  
  bin_order_df <- tibble::tibble(
    
    Bin =
      as.character(
        bptk_long$Bin
      ),
    
    BinNumber =
      bin_number
  ) %>%
    dplyr::distinct() %>%
    dplyr::arrange(
      is.na(BinNumber),
      BinNumber,
      Bin
    )
  
  
  bptk_long$Bin <-
    factor(
      
      bptk_long$Bin,
      
      levels =
        bin_order_df$Bin
    )
  
  
  # -------------------------------------------------------
  # Total assembly contribution of each bin
  # -------------------------------------------------------
  bin_total <- bptk_long %>%
    dplyr::group_by(
      Field,
      Treatment,
      Treatment_label,
      Panel,
      Bin
    ) %>%
    dplyr::summarise(
      
      TotalBinAssemblyContribution =
        sum(
          Contribution,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
    )
  
  
  # -------------------------------------------------------
  # Selection contribution of each bin
  #
  # Absolute contribution:
  #   HeS + HoS contribution to total community assembly.
  #
  # SelectionShare:
  #   fraction of ALL selection attributable to that bin.
  # -------------------------------------------------------
  selection_bins <- bptk_long %>%
    dplyr::filter(
      ProcessCode %in%
        c(
          "HeS",
          "HoS"
        )
    ) %>%
    dplyr::group_by(
      Field,
      Treatment,
      Treatment_label,
      Panel,
      Bin
    ) %>%
    dplyr::summarise(
      
      HeS_contribution =
        sum(
          Contribution[
            ProcessCode ==
              "HeS"
          ],
          na.rm = TRUE
        ),
      
      HoS_contribution =
        sum(
          Contribution[
            ProcessCode ==
              "HoS"
          ],
          na.rm = TRUE
        ),
      
      TotalSelectionContribution =
        sum(
          Contribution,
          na.rm = TRUE
        ),
      
      .groups =
        "drop"
    ) %>%
    dplyr::group_by(
      Field,
      Treatment
    ) %>%
    dplyr::mutate(
      
      TotalSelection =
        sum(
          TotalSelectionContribution,
          na.rm = TRUE
        ),
      
      SelectionShare =
        dplyr::if_else(
          
          TotalSelection > 0,
          
          TotalSelectionContribution /
            TotalSelection,
          
          NA_real_
        )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(
      Field,
      Treatment,
      dplyr::desc(
        TotalSelectionContribution
      )
    )
  
  
  outdir <-
    icamp_run$output_directory
  
  
  utils::write.csv(
    
    bptk_long,
    
    file.path(
      outdir,
      paste0(
        prefix,
        "_iCAMP_BPtk_bin_process_contributions_Harvest_vs_Storage.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  utils::write.csv(
    
    bin_total,
    
    file.path(
      outdir,
      paste0(
        prefix,
        "_iCAMP_total_bin_contribution_Harvest_vs_Storage.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  utils::write.csv(
    
    selection_bins,
    
    file.path(
      outdir,
      paste0(
        prefix,
        "_iCAMP_selection_bin_contribution_Harvest_vs_Storage.csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  list(
    
    long =
      bptk_long,
    
    bin_total =
      bin_total,
    
    selection =
      selection_bins
  )
}


ITS_bin_contribution_HS <-
  extract_icamp_bin_contributions_HS(
    ITS_icamp,
    prefix = "ITS"
  )



# ---------------------------------------------------------
# BIN FIGURE
#
# Each stacked bar = one phylogenetic bin.
#
# Height:
#   contribution of that bin to total community assembly.
#
# Colours:
#   process(es) to which that bin contributed across the
#   Harvest-vs-Storage turnovers.
# ---------------------------------------------------------



print(colnames(class_bin))
print(head(class_bin))

p_icamp_bins <- ggplot2::ggplot(
  
  ITS_bin_contribution_HS$long,
  
  ggplot2::aes(
    x = Bin,
    y = Contribution,
    fill = Process
  )
  
) +
  
  ggplot2::geom_col(
    
    width =
      0.82,
    
    color =
      "black",
    
    linewidth =
      0.15
  ) +
  
  ggplot2::facet_wrap(
    
    ~ Panel,
    
    ncol =
      3
  ) +
  
  ggplot2::scale_fill_manual(
    
    values =
      icamp_process_colours,
    
    drop =
      FALSE
  ) +
  
  ggplot2::scale_y_continuous(
    
    labels =
      scales::percent_format(
        accuracy = 0.1
      ),
    
    expand =
      ggplot2::expansion(
        mult = c(
          0,
          0.05
        )
      )
  ) +
  
  ggplot2::labs(
    
    x =
      "Phylogenetic bin",
    
    y =
      "Bin contribution to community assembly",
    
    fill =
      "Assembly process",
    
    subtitle =
      "Bin-level contribution to Harvest–Storage turnover"
  ) +
  
  theme_nature +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 9
      ),
    
    strip.text =
      ggplot2::element_text(
        size = 13,
        face = "bold"
      ),
    
    legend.position =
      "bottom",
    
    legend.title =
      ggplot2::element_text(
        size = 13,
        face = "bold"
      ),
    
    legend.text =
      ggplot2::element_text(
        size = 12
      ),
    
    plot.subtitle =
      ggplot2::element_text(
        hjust = 0.5,
        size = 14
      )
  )


p_icamp_bins


saveRDS(
  
  p_icamp_bins,
  
  file.path(
    ITS_icamp$output_directory,
    "ITS_iCAMP_bin_contributions_Harvest_vs_Storage_plot.RDS"
  )
)


ggplot2::ggsave(
  
  filename =
    file.path(
      ITS_icamp$output_directory,
      "ITS_iCAMP_bin_contributions_Harvest_vs_Storage.png"
    ),
  
  plot =
    p_icamp_bins,
  
  width =
    13,
  
  height =
    8,
  
  units =
    "in",
  
  dpi =
    600,
  
  bg =
    "white"
)



# =========================================================
# 9C. TABLE TO USE MANUALLY FOR SELECTION TEXT
# =========================================================
#
# This is the table I recommend using when writing the
# Results paragraph about bins under selection.
#
# TotalSelectionContribution:
#   absolute contribution of that bin to community assembly
#   through HeS + HoS.
#
# SelectionShare:
#   proportion of total selection attributable to that bin.
#
# Example interpretation:
#   SelectionShare = 0.42
#   -> that bin accounts for 42% of the selection signal
#      in that Field x Treatment Harvest-vs-Storage turnover.
# =========================================================

ITS_selection_bin_table <-
  ITS_bin_contribution_HS$selection %>%
  dplyr::mutate(
    
    TotalSelectionContribution_percent =
      100 *
      TotalSelectionContribution,
    
    SelectionShare_percent =
      100 *
      SelectionShare
  )


print(
  ITS_selection_bin_table
)


utils::write.csv(
  
  ITS_selection_bin_table,
  
  file.path(
    ITS_icamp$output_directory,
    "ITS_iCAMP_SELECTION_BINS_FOR_MANUSCRIPT.csv"
  ),
  
  row.names = FALSE
)



# =========================================================
# 9D. OPTIONAL SELECTION-ONLY BIN PLOT
# =========================================================
#
# Not required for the main manuscript figure.
# Useful for inspection / Supplementary material.
# =========================================================

p_icamp_selection_bins <- ggplot2::ggplot(
  
  ITS_selection_bin_table,
  
  ggplot2::aes(
    x = Bin,
    y = TotalSelectionContribution,
    fill = Treatment_label
  )
  
) +
  
  ggplot2::geom_col(
    width = 0.8
  ) +
  
  ggplot2::facet_wrap(
    ~ Field,
    nrow = 2,
    scales = "free_x"
  ) +
  
  ggplot2::scale_y_continuous(
    labels =
      scales::percent_format(
        accuracy = 0.1
      )
  ) +
  
  ggplot2::labs(
    
    x =
      "Phylogenetic bin",
    
    y =
      "Contribution to selection\n(HeS + HoS)",
    
    fill =
      "Treatment",
    
    subtitle =
      "Selection component of Harvest–Storage turnover"
  ) +
  
  theme_nature +
  
  ggplot2::theme(
    
    axis.text.x =
      ggplot2::element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 9
      ),
    
    legend.position =
      "bottom"
  )


p_icamp_selection_bins


saveRDS(
  
  p_icamp_selection_bins,
  
  file.path(
    ITS_icamp$output_directory,
    "ITS_iCAMP_selection_bin_contributions_plot.RDS"
  )
)

citation("iCAMP")

# =========================================================
# STATISTICAL TEST:
# Harvest versus Storage for assembly processes
# =========================================================

# Community-level iCAMP result already calculated
icamp_process <- ITS_icamp$output$bNRIiRCa


# Group samples by orchard x treatment x sampling
stage_groups <- data.frame(
  Group = paste(
    ITS_icamp$groups$Field,
    ITS_icamp$groups$Treatment,
    ITS_icamp$groups$Sampling,
    sep = "__"
  ),
  row.names = rownames(ITS_icamp$groups)
)


# Bootstrap comparison
set.seed(423542)

ITS_process_boot <- iCAMP::icamp.boot(
  icamp.result = icamp_process,
  treat = stage_groups,
  rand.time = 1000,
  compare = TRUE,
  between.group = FALSE,
  ST.estimation = FALSE
)

ITS_process_boot$summary
ITS_process_boot$compare

# =========================================================
# MANUSCRIPT TABLE:
# Harvest vs Storage comparison of iCAMP processes
# =========================================================

library(dplyr)
library(tidyr)
library(tibble)


# ---------------------------------------------------------
# 1. Lookup table for the groups used in icamp.boot()
# ---------------------------------------------------------

group_lookup <- ITS_icamp$groups %>%
  as.data.frame() %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::transmute(
    Group = paste(
      Field,
      Treatment,
      Sampling,
      sep = "__"
    ),
    Field = as.character(Field),
    Treatment = as.character(Treatment),
    Sampling = as.character(Sampling)
  ) %>%
  dplyr::distinct() %>%
  dplyr::mutate(
    Treatment_label = dplyr::case_when(
      
      Field == "Pfatten/Vadena" &
        Treatment == "Control" ~
        "Control",
      
      Field == "Pfatten/Vadena" &
        Treatment == "Geoxe" ~
        "Fludioxonil",
      
      Field == "Pfatten/Vadena" &
        Treatment == "Ulmasud" ~
        "Acidic clays",
      
      Field == "Sinich/Sinigo" &
        Treatment == "Control" ~
        "Control",
      
      Field == "Sinich/Sinigo" &
        Treatment == "Geoxe" ~
        "Captan + Fludioxonil",
      
      Field == "Sinich/Sinigo" &
        Treatment == "Ulmasud" ~
        "Captan + Acidic clays",
      
      TRUE ~ Treatment
    )
  )


# ---------------------------------------------------------
# 2. Observed process importance:
#    Harvest and Storage
# ---------------------------------------------------------

# ---------------------------------------------------------
# 2. Observed process importance:
#    Harvest and Storage
# ---------------------------------------------------------

boot_summary <- as.data.frame(
  ITS_process_boot$summary,
  check.names = FALSE
)

# iCAMP returns duplicated "Median" column names.
# Make all names unique so dplyr can work with the table.
names(boot_summary) <- make.unique(
  names(boot_summary),
  sep = "_"
)

# We only need these three columns anyway
boot_summary <- boot_summary %>%
  dplyr::select(
    Group,
    Process,
    Observed
  ) %>%
  dplyr::mutate(
    Observed = as.numeric(Observed),
    ProcessCode = standardize_icamp_process(Process)
  ) %>%
  dplyr::left_join(
    group_lookup,
    by = "Group"
  )


observed_table <- boot_summary %>%
  dplyr::filter(
    Sampling %in% c("Harvest", "Storage")
  ) %>%
  dplyr::select(
    Field,
    Treatment,
    Treatment_label,
    Sampling,
    ProcessCode,
    Observed
  ) %>%
  tidyr::pivot_wider(
    names_from = Sampling,
    values_from = Observed
  )

# ---------------------------------------------------------
# 3. Convert icamp.boot() pairwise comparisons to long form
# ---------------------------------------------------------

boot_compare <- ITS_process_boot$compare %>%
  
  tidyr::pivot_longer(
    cols = -c(Group1, Group2),
    names_to = c(
      "Process_raw",
      ".value"
    ),
    names_pattern =
      "^(.*)_(Relative\\.Diff|Cohen\\.d|Effect\\.Size|P\\.value)$"
  ) %>%
  
  dplyr::mutate(
    Relative.Diff =
      as.numeric(Relative.Diff),
    
    Cohen.d =
      as.numeric(Cohen.d),
    
    P.value =
      as.numeric(P.value),
    
    ProcessCode =
      standardize_icamp_process(Process_raw)
  )


# ---------------------------------------------------------
# 4. Attach metadata to Group1 and Group2
# ---------------------------------------------------------

lookup1 <- group_lookup %>%
  dplyr::select(
    Group1 = Group,
    Field1 = Field,
    Treatment1 = Treatment,
    Sampling1 = Sampling
  )


lookup2 <- group_lookup %>%
  dplyr::select(
    Group2 = Group,
    Field2 = Field,
    Treatment2 = Treatment,
    Sampling2 = Sampling
  )


boot_compare <- boot_compare %>%
  dplyr::left_join(
    lookup1,
    by = "Group1"
  ) %>%
  dplyr::left_join(
    lookup2,
    by = "Group2"
  )


# ---------------------------------------------------------
# 5. Keep ONLY Harvest vs Storage within the SAME
#    orchard and treatment
# ---------------------------------------------------------

HS_compare <- boot_compare %>%
  dplyr::filter(
    
    Field1 == Field2,
    
    Treatment1 == Treatment2,
    
    Sampling1 != Sampling2,
    
    Sampling1 %in% c(
      "Harvest",
      "Storage"
    ),
    
    Sampling2 %in% c(
      "Harvest",
      "Storage"
    )
    
  ) %>%
  
  # Orient effect sizes so POSITIVE always means
  # higher in Storage than Harvest
  dplyr::mutate(
    
    Cohen_d_Storage_minus_Harvest =
      dplyr::if_else(
        Sampling1 == "Storage",
        Cohen.d,
        -Cohen.d
      ),
    
    RelativeDiff_Storage_minus_Harvest =
      dplyr::if_else(
        Sampling1 == "Storage",
        Relative.Diff,
        -Relative.Diff
      )
  ) %>%
  
  dplyr::select(
    Field = Field1,
    Treatment = Treatment1,
    ProcessCode,
    Cohen_d_Storage_minus_Harvest,
    Effect.Size,
    P_value = P.value
  )


# ---------------------------------------------------------
# 6. Join observed values and statistics
# ---------------------------------------------------------

ITS_Harvest_Storage_statistics <- observed_table %>%
  
  dplyr::left_join(
    HS_compare,
    by = c(
      "Field",
      "Treatment",
      "ProcessCode"
    )
  ) %>%
  
  dplyr::mutate(
    
    Process =
      dplyr::recode(
        ProcessCode,
        "HeS" =
          "Heterogeneous selection",
        "HoS" =
          "Homogeneous selection",
        "DL" =
          "Dispersal limitation",
        "HD" =
          "Homogenizing dispersal",
        "DR" =
          "Drift"
      ),
    
    Harvest_percent =
      100 * Harvest,
    
    Storage_percent =
      100 * Storage,
    
    Change_percentage_points =
      100 * (Storage - Harvest),
    
    # Correct across all Harvest-vs-Storage
    # process comparisons
    P_FDR =
      p.adjust(
        P_value,
        method = "BH"
      ),
    
    Significance =
      dplyr::case_when(
        P_FDR < 0.001 ~ "***",
        P_FDR < 0.01  ~ "**",
        P_FDR < 0.05  ~ "*",
        TRUE ~ "ns"
      ),
    
    Direction =
      dplyr::case_when(
        
        P_FDR < 0.05 &
          Change_percentage_points > 0 ~
          "Higher after storage",
        
        P_FDR < 0.05 &
          Change_percentage_points < 0 ~
          "Lower after storage",
        
        TRUE ~
          "No significant change"
      )
  ) %>%
  
  dplyr::select(
    Field,
    Treatment = Treatment_label,
    Process,
    Harvest_percent,
    Storage_percent,
    Change_percentage_points,
    Cohen_d = Cohen_d_Storage_minus_Harvest,
    Effect_size = Effect.Size,
    P_value,
    P_FDR,
    Significance,
    Direction
  ) %>%
  
  dplyr::mutate(
    dplyr::across(
      c(
        Harvest_percent,
        Storage_percent,
        Change_percentage_points,
        Cohen_d
      ),
      ~ round(.x, 2)
    ),
    
    P_value =
      signif(
        P_value,
        3
      ),
    
    P_FDR =
      signif(
        P_FDR,
        3
      )
  )


ITS_Harvest_Storage_statistics
writexl::write_xlsx(
  ITS_Harvest_Storage_statistics,
  file.path(
    ITS_icamp$output_directory,
    "ITS_iCAMP_Harvest_vs_Storage_statistics.xlsx"
  )
)
