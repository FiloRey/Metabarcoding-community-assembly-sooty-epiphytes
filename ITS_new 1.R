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
set.seed(423542)

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
# STEP 4: ALPHA DIVERSITY ANALYSIS 
# ==============================
# 1) Rarefaction curves -------------------------------------------------
asv <- t(abundances(physeq_asv_filtered))
rarefy_depth <- min(rowSums(asv))
cat("Rarefying to:", rarefy_depth, "reads per sample\n")
rarefy_depth <- 15000
p_raref <- rarecurve(asv, step = 100, sample = rarefy_depth, tidy = TRUE)
rare_ITS <- ggplot(p_raref, aes(x = Sample, y = Species, group = Site, color = Site)) +
  geom_line(linewidth = 0.8, alpha = 0.7) +
  geom_vline(xintercept = rarefy_depth, linetype = "dashed", color = "black") +
  labs(
    x = "Sequencing Depth (reads)", 
    y = "Observed Species (ASVs/OTUs)"
  )

rare_ITS <- rare_ITS + theme(legend.position = "none")
rare_ITS
saveRDS(rare_ITS, "rare_ITS.RDS")

# 2) Alpha diversity calculation ----------------------------------------
ps0.rar <- rarefy_even_depth(
  physeq_asv_filtered, 
  sample.size = rarefy_depth, 
  rngseed = 123, 
  verbose = FALSE
)

hmp.div <- microbiome::alpha(ps0.rar, index = "all")
hmp.meta <- microbiome::meta(ps0.rar)
hmp.meta$sam_name <- rownames(hmp.meta)
hmp.div$sam_name <- rownames(hmp.div)

# Merge metrics with metadata
div.df <- merge(hmp.div, hmp.meta, by = "sam_name")

# Select and rename metrics
metrics <- c("Chao1","Shannon","Pielou","Dominance","Rarity")
div.df2 <- div.df[, c("sam_name","sampling","treatment","field","class", "block",
                      "chao1","diversity_shannon","evenness_pielou","dominance_relative","rarity_rare_abundance")]
colnames(div.df2) <- c("sam_name","Sampling","Treatment","Field","Class", "Block",
                       "Chao1","Shannon","Pielou","Dominance","Rarity")

# Ensure correct types
div.df2 <- div.df2 %>%
  mutate(
    Sampling  = factor(Sampling),
    Treatment = factor(Treatment),
    Block = factor(Block),
    Field     = factor(Field),
    Class     = factor(Class)
  )

# 3) Raincloud plot: All metrics, Harvest vs Storage ----------------------

metrics <- c("Chao1","Shannon","Pielou","Dominance","Rarity")

#Mixed Models
div_sub <- div.df2 %>%
  dplyr::select(sam_name, Treatment, Sampling, Field, Block, Chao1, Pielou, Rarity, Shannon, Dominance)

m_Chao1 <- lm(
  Chao1 ~ Sampling + Field + Treatment,
  data = div_sub
)

summary(m_Chao1)
anova(m_Chao1)

# shannon
m_shannon <- lm(
  Shannon ~ Sampling + Field + Treatment,
  data = div_sub
)

summary(m_shannon)
anova(m_shannon)

m_Pielou <- lm(
  Pielou ~ Sampling + Field + Treatment,
  data = div_sub
)

summary(m_Pielou)
anova(m_Pielou)


# Dominance
m_dominance <- lm(
  Dominance ~ Sampling + Field + Treatment,
  data = div_sub
)

summary(m_dominance)
anova(m_dominance)

# Dominance
m_Rarity <- lm(
  Rarity ~ Sampling + Field + Treatment,
  data = div_sub
)

summary(m_Rarity)
anova(m_Rarity)

par(mfrow = c(2,2))

plot(m_shannon, which = 1)
qqnorm(residuals(m_shannon)); qqline(residuals(m_shannon))

plot(m_dominance, which = 1)
qqnorm(residuals(m_dominance)); qqline(residuals(m_dominance))

par(mfrow = c(1,1))


# Function to collapse each metric into a single string with newlines
make_label <- function(txt_vec) {
  paste(txt_vec, collapse = "\n")
}


shannon_results <- div_sub %>%
  group_by(Field) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(Shannon ~ Sampling + Treatment, data = .x)),
    tidy  = map(model, broom::tidy)
  ) %>%
  unnest(tidy)

shannon_results

shannon_sampling <- shannon_results %>%
  filter(term == "SamplingStorage") %>%
  dplyr::select(Field, estimate, p.value)

shannon_sampling
labels_shannon <- shannon_sampling %>%
  mutate(
    label = paste0(
      "Δ = ", round(estimate,3),
      "\np = ", signif(p.value,3)
    )
  )
p_shannon <- ggplot(div_sub, aes(x = Sampling, y = Shannon, fill = Sampling)) +
  stat_halfeye(
    adjust = 0.6,
    width = 0.6,
    justification = -0.25,
    alpha = 0.6,
    slab_color = NA
  ) +
  geom_boxplot(width = 0.2, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, size = 1.6, alpha = 0.6) +
  facet_wrap(~Field) +
  geom_text(
    data = labels_shannon,
    aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE,
    vjust = 1.5
  ) +
  scale_fill_manual(values = c(Harvest = "#E64B35", Storage = "#4DBBD5")) +
  labs(x = "Sampling", y = "Shannon") +
  theme_nature +
  theme(legend.position = "none")

p_shannon

# -------------------------------
# Subset data per field
# -------------------------------
shannon_pfatten <- div.df2 %>% 
  filter(Field == "Pfatten/Vadena") %>%
  dplyr::select(Sampling, Treatment, Shannon) %>%
  dplyr::rename(Value = Shannon)

shannon_sinich <- div.df2 %>% 
  filter(Field == "Sinich/Sinigo") %>%
  dplyr::select(Sampling, Treatment, Shannon) %>%
  dplyr::rename(Value = Shannon)

# -------------------------------
# Pfatten/Vadena: ART ANOVA
# -------------------------------
# Harvest
shannon_pfatten_harvest <- shannon_pfatten %>% filter(Sampling == "Harvest")
art_pf_harvest <- art(Value ~ Treatment, data = shannon_pfatten_harvest)
aov_pf_harvest <- anova(art_pf_harvest)
aov_pf_harvest  # view results

# Posthoc
ph_pf_harvest <- shannon_pfatten_harvest %>%
  pairwise_wilcox_test(Value ~ Treatment)
ph_pf_harvest

# Storage
shannon_pfatten_storage <- shannon_pfatten %>% filter(Sampling == "Storage")
art_pf_storage <- art(Value ~ Treatment, data = shannon_pfatten_storage)
aov_pf_storage <- anova(art_pf_storage)
aov_pf_storage

ph_pf_storage <- shannon_pfatten_storage %>%
  pairwise_wilcox_test(Value ~ Treatment)
ph_pf_storage

# -------------------------------
# Sinich/Sinigo: ART ANOVA
# -------------------------------
# Harvest
shannon_sinich_harvest <- shannon_sinich %>% filter(Sampling == "Harvest")
art_si_harvest <- art(Value ~ Treatment, data = shannon_sinich_harvest)
aov_si_harvest <- anova(art_si_harvest)
aov_si_harvest

ph_si_harvest <- shannon_sinich_harvest %>%
  pairwise_wilcox_test(Value ~ Treatment, p.adjust.method = "bonferroni")
ph_si_harvest

shannon_sinich_harvest %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean_shannon = mean(Value, na.rm = TRUE),
    median_shannon = median(Value, na.rm = TRUE),
    sd_shannon = sd(Value, na.rm = TRUE)
  )
# Storage
shannon_sinich_storage <- shannon_sinich %>% filter(Sampling == "Storage")
art_si_storage <- art(Value ~ Treatment, data = shannon_sinich_storage)
aov_si_storage <- anova(art_si_storage)
aov_si_storage

ph_si_storage <- shannon_sinich_storage %>%
  pairwise_wilcox_test(Value ~ Treatment, p.adjust.method = "bonferroni")
ph_si_storage

# ==============================
# STEP 5: BETA DIVERSITY 
# ==============================

# 1) CSS normalization ----------------------------------------------------
tse__2 <- tse__fresh
# Agglomerate to Genus 
tse__2 <- agglomerateByRank(tse__2, rank = "Genus", update.tree = TRUE)
# Filter low-abundance samples
counts <- assay(tse__2, "counts") 
keep_samples <- colSums(counts > 0) > 1 
counts <- counts[, keep_samples] 
tse__2 <- tse__2[, keep_samples]
# Create MRexperiment and CSS normalization 
mr_obj <- newMRexperiment(counts)
p <- cumNormStatFast(mr_obj) 
mr_obj <- cumNorm(mr_obj, p = p)
# Extract normalized log2 counts
norm_counts <- MRcounts(mr_obj, norm = TRUE, log = F) 
norm_counts <- norm_counts[rownames(tse__2), colnames(tse__2)] 
assay(tse__2, "css_norm") <- norm_counts

meta <- as.data.frame(colData(tse__2))
meta <- meta[colnames(norm_counts), ]
meta$sampling  <- factor(meta$sampling)
meta$treatment <- factor(meta$treatment)
meta$field     <- factor(meta$field)


# 2) dbRDA -------------------------------------------------------------------
dbrda_all <- capscale(t(norm_counts) ~ sampling + treatment + field,
                      data = meta, distance = "bray")

dbrda_scores <- as.data.frame(scores(dbrda_all, display = "sites"))
dbrda_scores$Sampling  <- meta$sampling
dbrda_scores$Treatment <- meta$treatment
dbrda_scores$Field     <- meta$field


fields <- levels(meta$field)
sampling_levels <- levels(meta$sampling)


# 3) PERMANOVA ------------------------------------------------------------
field_results <- list()

for(f in fields){
  
  idx <- meta$field == f
  counts_sub <- norm_counts[, idx]
  meta_sub <- meta[idx, ]
  
  dist_sub <- vegdist(t(counts_sub), method = "bray")
  
  # PERMANOVA: sampling + treatment
  permanova_res <- adonis2(
    dist_sub ~ sampling + treatment,
    data = meta_sub,
    permutations = 999,
    by = "margin"
  )
  
  # PERMDISP for sampling
  permdisp_sampling <- betadisper(dist_sub, meta_sub$sampling)
  perm_sampling <- permutest(permdisp_sampling)
  
  # PERMDISP for treatment
  permdisp_treatment <- betadisper(dist_sub, meta_sub$treatment)
  perm_treatment <- permutest(permdisp_treatment)
  
  field_results[[f]] <- list(
    permanova = permanova_res,
    permdisp_sampling = perm_sampling,
    permdisp_treatment = perm_treatment
  )
}


field_sampling_results <- list()

for(f in fields){
  for(s in sampling_levels){
    key <- paste(f, s, sep="_")
    
    idx <- meta$field == f & meta$sampling == s
    counts_sub <- norm_counts[, idx]
    meta_sub <- meta[idx, ]
    
    dist_sub <- vegdist(t(counts_sub), method = "bray")
    
    # PERMANOVA: treatment only
    permanova_res <- adonis2(
      dist_sub ~ treatment,
      data = meta_sub,
      permutations = 999,
      by = "margin"
    )
    
    # PERMDISP: treatment only
    permdisp_treatment <- betadisper(dist_sub, meta_sub$treatment)
    perm_treatment <- permutest(permdisp_treatment)
    
    field_sampling_results[[key]] <- list(
      permanova = permanova_res,
      permdisp_treatment = perm_treatment
    )
  }
}


# Example: PERMANOVA + PERMDISP for Pfatten/Vadena (sampling + treatment)
field_results[["Pfatten/Vadena"]]$permanova
field_results[["Pfatten/Vadena"]]$permdisp_sampling
field_results[["Pfatten/Vadena"]]$permdisp_treatment

field_results[["Sinich/Sinigo"]]$permanova
field_results[["Sinich/Sinigo"]]$permdisp_sampling
field_results[["Sinich/Sinigo"]]$permdisp_treatment

# Example: PERMANOVA + PERMDISP for Pfatten/Vadena Harvest (treatment only)
field_sampling_results[["Pfatten/Vadena_Harvest"]]$permanova
field_sampling_results[["Pfatten/Vadena_Harvest"]]$permdisp_treatment
field_sampling_results[["Sinich/Sinigo_Harvest"]]$permanova
field_sampling_results[["Sinich/Sinigo_Harvest"]]$permdisp_treatment
# Example: PERMANOVA + PERMDISP for Sinich/Sinigo Storage (treatment only)
field_sampling_results[["Pfatten/Vadena_Storage"]]$permanova
field_sampling_results[["Pfatten/Vadena_Storage"]]$permdisp_treatment
field_sampling_results[["Sinich/Sinigo_Storage"]]$permanova
field_sampling_results[["Sinich/Sinigo_Storage"]]$permdisp_treatment



# 4) Ellipses -------------------------------------------------------------
vegan_ellipse <- function(df, group_col, axes = c("CAP1","CAP2"), level = 0.68){
  groups <- levels(df[[group_col]])
  ellipses <- lapply(groups, function(g){
    sub <- df[df[[group_col]]==g, axes]
    cov_mat <- cov(sub)
    center <- colMeans(sub)
    angles <- seq(0, 2*pi, length.out=100)
    ellipse <- t(sapply(angles, function(theta){
      center + sqrt(qchisq(level, 2)) * t(chol(cov_mat)) %*% c(cos(theta), sin(theta))
    }))
    ellipse <- as.data.frame(ellipse)
    ellipse$Group <- g
    colnames(ellipse)[1:2] <- axes
    ellipse
  })
  do.call(rbind, ellipses)
}

dbrda_ellipses <- vegan_ellipse(dbrda_scores, group_col = "Sampling", level = 0.68)

# PLOT --------------------------------------------------------------------
dbrda_plot <- ggplot(dbrda_scores, aes(x = CAP1, y = CAP2)) +
  geom_point(aes(color = Sampling, shape = Treatment), size = 3) +
  geom_path(
    data = dbrda_ellipses,
    aes(x = CAP1, y = CAP2, color = Group),
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_color_manual(values = c(Harvest = "#E64B35", Storage = "#4DBBD5")) +
  scale_shape_manual(values = 1:length(levels(meta$treatment))) +
  labs(x = "CAP1", y = "CAP2", color = "Sampling", shape = "Treatment") +
  facet_wrap(~ Field) +
  theme_minimal(base_size = 14) +
  theme_nature +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank()
  )



dbrda_plot_no_legend <- dbrda_plot + theme(legend.position = "none")
dbrda_plot


# 5) combine with Shannon -------------------------------------------------
p_shannon <- p_shannon +
  theme(axis.title.x = element_blank())
p_shannon

combined_plot <- plot_grid(
  p_shannon,  # previously prepared alpha-diversity plot
  dbrda_plot_no_legend,
  nrow = 2,
  labels = c("a","c"),
  rel_heights = c(1.2,1.8)
)

combined_plot

saveRDS(combined_plot, file = "combined_ITS_plot.rds")

legend <- get_legend(
 dbrda_plot + theme(legend.position = "bottom")
)
legend
# Save legend as RDS
saveRDS(legend, file = "combined_ITS_legend.rds")

# ==============================
# STEP 6: BETA PARTITIONING 
# ==============================
pseq.rel <- microbiome::transform(physeq_asv_filtered, "compositional")

otu_pa <- as(otu_table(pseq.rel), "matrix")
otu_pa <- 1 * (otu_pa > 0)

if(taxa_are_rows(pseq.rel)) {
  otu_pa <- t(otu_pa)
}

meta <- data.frame(sample_data(pseq.rel))

beta_summary_list <- list()
mantel_list <- list()
beta_long_list <- list()

fields <- unique(meta$field)
treatments <- unique(meta$treatment)

for(f in fields){
  
  for(t in treatments){
    
    idx <- meta$field == f & meta$treatment == t
    
    if(sum(idx) < 4) next
    
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
    within_H <- beta_sor[H, H][upper.tri(beta_sor[H,H])]
    within_S <- beta_sor[S, S][upper.tri(beta_sor[S,S])]
    
    between <- beta_sor[H, S]
    turnover <- beta_sim[H, S]
    nestedness <- beta_nes[H, S]
    
    key <- paste(f,t,sep="_")
    
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
    
    stage_numeric <- ifelse(stage=="Harvest",0,1)
    stage_dist <- dist(stage_numeric)
    
    mantel_sor <- mantel(beta$beta.sor, stage_dist, permutations=999)
    mantel_sim <- mantel(beta$beta.sim, stage_dist, permutations=999)
    mantel_nes <- mantel(beta$beta.sne, stage_dist, permutations=999)
    
    mantel_list[[key]] <- data.frame(
      Field = f,
      Treatment = t,
      Component = c("βSOR","βSIM","βNES"),
      R = c(mantel_sor$statistic,
            mantel_sim$statistic,
            mantel_nes$statistic),
      P_value = c(mantel_sor$signif,
                  mantel_sim$signif,
                  mantel_nes$signif)
    )
    
    # ---------------------------
    # Extract pairwise distances
    # ---------------------------
    
    extract_pairs <- function(mat, stage_vector, stage_name){
      
      idx_stage <- which(stage_vector == stage_name)
      vals <- mat[idx_stage, idx_stage]
      vals <- vals[upper.tri(vals)]
      
      data.frame(
        Stage = stage_name,
        Value = vals
      )
    }
    
    beta_long <- bind_rows(
      extract_pairs(beta_sor, stage, "Harvest") %>% mutate(Component="βSOR"),
      extract_pairs(beta_sim, stage, "Harvest") %>% mutate(Component="βSIM"),
      extract_pairs(beta_nes, stage, "Harvest") %>% mutate(Component="βNES"),
      extract_pairs(beta_sor, stage, "Storage") %>% mutate(Component="βSOR"),
      extract_pairs(beta_sim, stage, "Storage") %>% mutate(Component="βSIM"),
      extract_pairs(beta_nes, stage, "Storage") %>% mutate(Component="βNES")
    )
    
    beta_long$Field <- f
    beta_long$Treatment <- t
    
    beta_long_list[[key]] <- beta_long
    
  }
}

beta_summary <- bind_rows(beta_summary_list)
mantel_df <- bind_rows(mantel_list)
beta_long_all <- bind_rows(beta_long_list)

beta_means <- beta_long_all %>%
  group_by(Field, Treatment, Component, Stage) %>%
  summarise(Value = mean(Value, na.rm=TRUE), .groups="drop")
beta_means
mantel_labels <- mantel_df %>%
  mutate(
    label = paste0(
      "R = ", round(R,2),
      "\np = ", signif(P_value,2)
    )
  )

label_positions <- beta_long_all %>%
  group_by(Field,Treatment,Component) %>%
  summarise(
    Stage="Storage",
    Value=max(Value,na.rm=TRUE)*0.95,
    .groups="drop"
  )

mantel_labels <- dplyr::left_join(mantel_labels,label_positions,
                                  by=c("Field","Treatment","Component"))

p_beta <- ggplot(beta_long_all, aes(Stage, Value)) +
  
  geom_jitter(aes(color=Stage),
              width=0.08,
              alpha=0.35,
              size=1.5) +
  
  geom_line(data=beta_means,
            aes(group=interaction(Field,Treatment,Component)),
            linewidth=1.1,
            color="black") +
  
  geom_point(data=beta_means,
             size=4,
             shape=21,
             fill="white",
             color="black",
             stroke=1.2) +
  
  geom_text(data=mantel_labels,
            aes(Stage, Value, label=label),
            inherit.aes=FALSE,
            size=3.5,
            fontface="bold",
            hjust=1) +
  
  facet_grid(Field + Treatment ~ Component) +
  
  scale_color_manual(values=c(
    Harvest="#E64B35",
    Storage="#4DBBD5"
  )) +
  
  labs(
    x=NULL,
    y="Beta diversity"
  ) +
  
  theme_nature
p_beta

# Filter beta_long_all for βSIM component
beta_sim_long <- beta_long_all %>%
  filter(Component == "βSIM")

# Compute mean values per Field × Treatment × Stage
beta_sim_means <- beta_sim_long %>%
  group_by(Field, Treatment, Stage) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop")

# Prepare Mantel labels for βSIM only
mantel_sim_labels <- mantel_df %>%
  filter(Component == "βSIM") %>%
  mutate(
    label = paste0("R = ", round(R, 2), "\np = ", signif(P_value, 2))
  )

# Position labels at 95% of max value per Field × Treatment
label_positions <- beta_sim_long %>%
  group_by(Field, Treatment) %>%
  summarise(
    Stage = "Storage",
    Value = max(Value, na.rm = TRUE) * 0.95,
    .groups = "drop"
  )

mantel_sim_labels <- dplyr::left_join(mantel_sim_labels, label_positions, 
                                      by = c("Field", "Treatment"))


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
    ~ Treatment,
    nrow = 1,
    labeller = labeller(
      Treatment = c(
        Control = "Control",
        Geoxe = "Fludioxonil",
        Ulmasud = "Plant fortifier"
      )
    )
  ) +
  
  scale_color_manual(values = c(
    Harvest = "#E64B35",
    Storage = "#4DBBD5"
  )) +
  
  labs(
    title = "Pfatten/Vadena",
    x = NULL,
    y = expression(beta[SIM]~"(turnover)")
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
    ~ Treatment,
    nrow = 1,
    labeller = labeller(
      Treatment = c(
        Control = "Control",
        Geoxe = "Captan + Fludioxonil",
        Ulmasud = "Captan + Plant fortifier"
      )
    )
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
  ) + labs(
    title = "Sinich/Sinigo",
    x = NULL,
    y = NULL) +
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
# =========================================================

fields_order <- c("Pfatten/Vadena", "Sinich/Sinigo")
treatments_order <- c("Control", "Geoxe", "Ulmasud")

treatment_labels <- c(
  "Pfatten/Vadena__Control" = "Control",
  "Pfatten/Vadena__Geoxe"   = "Fludioxonil",
  "Pfatten/Vadena__Ulmasud" = "Plant fortifier",
  "Sinich/Sinigo__Control"  = "Control",
  "Sinich/Sinigo__Geoxe"    = "Captan + Fludioxonil",
  "Sinich/Sinigo__Ulmasud"  = "Captan + Plant fortifier"
)

lefse_long_list <- list()

for(f in fields_order){
  for(t in treatments_order){
    
    key <- paste(f, t, sep = "_")
    
    tse_treat <- tse__fresh[, tse__fresh$field == f & tse__fresh$treatment == t]
    if(ncol(tse_treat) < 3) next
    
    tse_genus <- agglomerateByRank(tse_treat, rank = "Genus", update.tree = TRUE)
    
    keep_var <- apply(assay(tse_genus), 1, var) > 0
    tse_genus <- tse_genus[keep_var, ]
    if(nrow(tse_genus) == 0) next
    
    prevalence <- rowSums(assay(tse_genus) > 0) / ncol(tse_genus)
    tse_genus <- tse_genus[prevalence >= 0.7, ]
    if(nrow(tse_genus) == 0) next
    
    tse_ra <- relativeAb(tse_genus)
    max_abund <- apply(assay(tse_ra), 1, max)
    tse_ra <- tse_ra[max_abund >= 0.001, ]
    if(nrow(tse_ra) == 0) next
    
    tn <- get_terminal_nodes(rownames(tse_ra))
    tse_final <- tse_ra[tn, ]
    if(nrow(tse_final) == 0) next
    
    sampling_vals <- colData(tse_final)$sampling
    valid <- !is.na(sampling_vals)
    tse_final <- tse_final[, valid]
    sampling_vals <- sampling_vals[valid]
    
    if(length(unique(sampling_vals)) < 2) next
    
    res_lefse <- tryCatch(
      lefser(tse_final, classCol = "sampling"),
      error = function(e) NULL
    )
    
    if(is.null(res_lefse) || nrow(res_lefse) == 0) next
    
    res_lefse$features <- gsub("_gen_Incertae_sedis", "", res_lefse$features, fixed = TRUE)
    
    res_lefse <- res_lefse %>%
      mutate(
        Field = f,
        Treatment = t,
        TreatmentLabel = treatment_labels[paste(f, t, sep = "__")],
        Enriched = ifelse(scores < 0, "Harvest", "Storage"),
        LDA = abs(scores),
        SignedLDA = ifelse(scores < 0, -abs(scores), abs(scores)),
        Taxon = features
      ) %>%
      dplyr::select(Field, Treatment, TreatmentLabel, Taxon, Enriched, LDA, SignedLDA)
    
    lefse_long_list[[key]] <- res_lefse
  }
}

lefse_long <- bind_rows(lefse_long_list)

# Keep strongest taxa per field separately

top_taxa_field <- lefse_long %>%
  group_by(Field, Taxon) %>%
  summarise(maxLDA = max(abs(SignedLDA), na.rm = TRUE), .groups = "drop") %>%
  group_by(Field) %>%
  slice_max(order_by = maxLDA, n = 15, with_ties = FALSE) %>%
  ungroup()

lefse_top <- lefse_long %>%
  semi_join(top_taxa_field, by = c("Field", "Taxon"))

# Complete field × treatment × taxon combinations
#    using the top taxa chosen separately for each field

plot_grid <- expand.grid(
  Field = fields_order,
  Treatment = treatments_order,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    top_taxa_field %>% dplyr::select(Field, Taxon),
    by = "Field"
  ) %>%
  mutate(
    TreatmentLabel = treatment_labels[paste(Field, Treatment, sep = "__")]
  )

lefse_plot_df <- plot_grid %>%
  dplyr::left_join(
    lefse_top,
    by = c("Field", "Treatment", "TreatmentLabel", "Taxon")
  ) %>%
  mutate(
    Enriched = factor(Enriched, levels = c("Harvest", "Storage")),
    LDA = ifelse(is.na(LDA), 0, LDA),
    SignedLDA = ifelse(is.na(SignedLDA), 0, SignedLDA),
    Taxon = str_replace_all(Taxon, "_", " "),
    TreatmentLabel = case_when(
      Field == "Pfatten/Vadena" & Treatment == "Control" ~ "Control",
      Field == "Pfatten/Vadena" & Treatment == "Geoxe" ~ "Fludioxonil",
      Field == "Pfatten/Vadena" & Treatment == "Ulmasud" ~ "Plant fortifier",
      Field == "Sinich/Sinigo" & Treatment == "Control" ~ "Control",
      Field == "Sinich/Sinigo" & Treatment == "Geoxe" ~ "Captan +\nFludioxonil",
      Field == "Sinich/Sinigo" & Treatment == "Ulmasud" ~ "Captan +\nPlant fortifier",
      TRUE ~ TreatmentLabel
    )
  )

# Shared taxon order across both fields
#    This keeps common taxa on the same line

global_taxon_order <- lefse_plot_df %>%
  group_by(Taxon) %>%
  summarise(maxLDA = max(LDA, na.rm = TRUE), .groups = "drop") %>%
  arrange(maxLDA) %>%
  pull(Taxon)

lefse_plot_df <- lefse_plot_df %>%
  mutate(
    Taxon = factor(Taxon, levels = global_taxon_order),
    TreatmentLabel = factor(
      TreatmentLabel,
      levels = c(
        "Control", "Fludioxonil", "Plant fortifier",
        "Captan +\nFludioxonil", "Captan +\nPlant fortifier"
      )
    ),
    Field = factor(Field, levels = c("Pfatten/Vadena", "Sinich/Sinigo"))
  )

# Compact aligned bubble plot

p_lefse_bubble_aligned <- ggplot(
  lefse_plot_df,
  aes(x = TreatmentLabel, y = Taxon)
) +
  geom_point(
    data = lefse_plot_df %>% filter(!is.na(Enriched)),
    aes(size = pmax(LDA, 2), fill = Enriched),
    shape = 21,
    color = "black",
    stroke = 0.25,
    alpha = 0.95
  ) +
  facet_grid(. ~ Field, scales = "free_x", space = "free_x") +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(
    values = c(Harvest = "#E64B35", Storage = "#4DBBD5"),
    breaks = c("Harvest", "Storage"),
    na.value = "transparent"
  ) +
  scale_size_continuous(
    range = c(4, 12),
    breaks = c(2, 3, 4, 5),
    name = "|LDA|"
  ) +
  labs(
    x = NULL,
    y = NULL,
    fill = "Enriched in",
    size = "|LDA|"
  ) +
  theme_nature +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p_lefse_bubble_aligned
top_taxa_field

# Final combined figure
p_beta_final <- p_beta_final &
  theme(
    strip.text = element_text(face = "bold", size = 16),
    strip.background = element_blank(),
    strip.placement = "outside"
  )


p_final_combined <-
  wrap_elements(p_beta_final) +
  wrap_elements(p_lefse_bubble_aligned) +
  plot_layout(ncol = 1, heights = c(0.32, 0.68)) +
  plot_annotation(tag_levels = list(c("a", "b"))) &
  theme(
    plot.tag = element_text(size = 12, face = "bold"),
    plot.tag.position = c(0,1)
  )

p_final_combined

#  Save
ggsave(
  "final_compact_manuscript_figure.png",
  p_final_combined,
  width = 7.0,
  height = 6.2,
  units = "in",
  dpi = 300
)


lefse_all_list <- list()

tse_all <- tse__fresh

# Agglomerate to Genus
tse_genus <- agglomerateByRank(tse_all, rank = "Genus", update.tree = TRUE)

# Remove zero-variance taxa
keep_var <- apply(assay(tse_genus), 1, var) > 0
tse_genus <- tse_genus[keep_var, ]
if(nrow(tse_genus) == 0) stop("No taxa left after variance filtering")

# Prevalence filter
prevalence <- rowSums(assay(tse_genus) > 0) / ncol(tse_genus)
tse_genus <- tse_genus[prevalence >= 0.7, ]
if(nrow(tse_genus) == 0) stop("No taxa left after prevalence filtering")

# Relative abundance
tse_ra <- relativeAb(tse_genus)

# Minimum abundance filter
max_abund <- apply(assay(tse_ra), 1, max)
tse_ra <- tse_ra[max_abund >= 0.001, ]
if(nrow(tse_ra) == 0) stop("No taxa left after abundance filtering")

# Keep terminal nodes only
tn <- get_terminal_nodes(rownames(tse_ra))
tse_final <- tse_ra[tn, ]
if(nrow(tse_final) == 0) stop("No taxa left after terminal node filtering")

# Remove samples with missing sampling annotation
sampling_vals <- colData(tse_final)$sampling
valid <- !is.na(sampling_vals)
tse_final <- tse_final[, valid]
sampling_vals <- sampling_vals[valid]

if(length(unique(sampling_vals)) < 2) {
  stop("Need at least two sampling groups for LEfSe")
}

# Run LEfSe on all data
res_lefse <- tryCatch(
  lefser(tse_final, classCol = "sampling"),
  error = function(e) NULL
)

if(is.null(res_lefse) || nrow(res_lefse) == 0) {
  stop("LEfSe returned no significant features")
}

# Clean and format output
res_lefse$features <- gsub("_gen_Incertae_sedis", "", res_lefse$features, fixed = TRUE)

lefse_all <- res_lefse %>%
  mutate(
    Enriched = ifelse(scores < 0, "Harvest", "Storage"),
    LDA = abs(scores),
    SignedLDA = ifelse(scores < 0, -abs(scores), abs(scores)),
    Taxon = features
  ) %>%
  dplyr::select(Taxon, Enriched, LDA, SignedLDA)

lefse_all

lefse_all_plot <- lefse_all %>%
  filter(LDA > 0) %>%
  mutate(
    Enriched = factor(Enriched, levels = c("Harvest", "Storage")),
    Taxon = fct_reorder(Taxon, SignedLDA)
  )

p_lefse_all <- ggplot(lefse_all_plot, aes(x = Taxon, y = SignedLDA, fill = Enriched)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_fill_manual(values = c(
    Harvest = "#E64B35",
    Storage = "#4DBBD5"
  )) +
  labs(
    x = NULL,
    y = "Signed LEfSe score (LDA)",
    fill = NULL
  ) +
  theme_nature +
  theme(
    legend.position = "top"
  )

p_lefse_all


saveRDS(lefse_all, "lefse_result.rds")


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


# DOC ---------------------------------------------------------------------
# Transform to compositional 
physeq <- physeq_asv_filtered
physeq <- tax_glom(physeq, taxrank = "Genus")
physeq <- prune_taxa(taxa_sums(physeq) > 100, physeq)
physeq_rel <- microbiome::transform(physeq, "compositional")

# --- Subset samples ---

Harvest <- subset_samples(physeq_rel, sampling == "Harvest")
Storage <- subset_samples(physeq_rel, sampling == "Storage")

Harvest_none <- subset_samples(Harvest, treatment == "Control")
Harvest_ulmasud <- subset_samples(Harvest, treatment == "Ulmasud")
Harvest_geoxe <- subset_samples(Harvest, treatment == "Geoxe")

Storage_none <- subset_samples(Storage, treatment == "Control")
Storage_ulmasud <- subset_samples(Storage, treatment == "Ulmasud")
Storage_geoxe <- subset_samples(Storage, treatment == "Geoxe")


doc <- DOC(otu_table(physeq_rel), R =1000, iterations = 10, cores = 6)
doc_Harvest <- DOC(otu_table(Harvest), R =1000, iterations = 10, cores = 6)  
doc_Storage <- DOC(otu_table(Storage), R =1000, iterations = 10, cores = 6)  

doc_Harvest_none <- DOC(otu_table(Harvest_none), R =1000, iterations = 10, cores = 6)  
doc_Harvest_ulmasud <- DOC(otu_table(Harvest_ulmasud), R =1000, iterations = 10, cores = 6)  
doc_Harvest_geoxe <- DOC(otu_table(Harvest_geoxe), R =1000, iterations = 10, cores = 6)  

doc_Storage_none <- DOC(otu_table(Storage_none), R =1000, iterations = 10, cores = 6)  
doc_Storage_ulmasud <- DOC(otu_table(Storage_ulmasud), R =1000, iterations = 10, cores = 6)  
doc_Storage_geoxe <- DOC(otu_table(Storage_geoxe), R =1000, iterations = 10, cores = 6) 

plot(doc)
plot(doc_Harvest)
plot(doc_Storage)

results.null_Harvest <- DOC.null(otu_table(Harvest))
results.null <- DOC.null(otu_table(Storage))

# Inspect structure of doc objects returned by the package's DOC()
str(doc_Harvest)
str(doc_Storage)
# Check DO and CI names
head(doc_Harvest$DO)
head(doc_Harvest$CI)
names(doc_Harvest)

analyze_DOC_combined_fromDOC <- function(doc1, doc2,
                                         title1 = "Harvest", title2 = "Storage") {
  # helper to find CI column names robustly
  find_CI_cols <- function(CI) {
    cn <- colnames(CI)
    overlap_col <- if ("Overlap" %in% cn) "Overlap" else cn[1]
    # try to find lower (2.5), median (50), upper (97.5)
    lower_col <- cn[grep("2\\.5|2_5|2p5|X2\\.5", cn)[1]]
    median_col <- cn[grep("(^|[^0-9])50([^0-9]|$)|X50|50\\.", cn)[1]]
    upper_col <- cn[grep("97\\.5|97_5|97p5|X97\\.5", cn)[1]]
    # fallbacks
    if (is.na(lower_col) && length(cn) >= 2) lower_col <- cn[2]
    if (is.na(median_col) && length(cn) >= 3) median_col <- cn[3]
    if (is.na(upper_col) && length(cn) >= 4) upper_col <- cn[4]
    list(overlap = overlap_col, lower = lower_col, median = median_col, upper = upper_col)
  }
  
  # safe extractor for metrics from a DOC object
  extract_metrics_from_DOC <- function(doc_obj) {
    # checks
    if (is.null(doc_obj$DO)) stop("doc_obj has no DO component")
    DO <- as.data.frame(doc_obj$DO)
    CI <- if (!is.null(doc_obj$CI)) as.data.frame(doc_obj$CI) else NULL
    LOW <- if (!is.null(doc_obj$LOWESS)) as.data.frame(doc_obj$LOWESS) else NULL
    LME  <- if (!is.null(doc_obj$LME)) as.data.frame(doc_obj$LME) else NULL
    NEG  <- if (!is.null(doc_obj$NEG)) as.data.frame(doc_obj$NEG) else NULL
    FNS  <- if (!is.null(doc_obj$FNS)) as.data.frame(doc_obj$FNS) else NULL
    
    # Oc: prefer median of NEG$Neg.Slope (bootstrap of Oc). fallback: compute from LOWY/LOWESS.
    Oc <- NA_real_
    if (!is.null(NEG)) {
      # find any numeric column likely to be Neg.Slope
      neg_col <- grep("Neg|neg|Oc|Oc\\.?|Neg.Slope", colnames(NEG), value = TRUE)[1]
      if (!is.na(neg_col)) Oc <- median(NEG[[neg_col]], na.rm = TRUE)
      else Oc <- median(unlist(NEG[, sapply(NEG, is.numeric), drop = FALSE]), na.rm = TRUE)
      
    } else if (!is.null(LOW)) {
      lowx <- LOW[[1]]
      
      # Detect LOWY or LOWESS column dynamically
      lowy <- if (any(grepl("LOWY", colnames(LOW), ignore.case = TRUE))) {
        LOW[[grep("LOWY", colnames(LOW), ignore.case = TRUE)]]
      } else if (any(grepl("LOWESS", colnames(LOW), ignore.case = TRUE))) {
        LOW[[grep("LOWESS", colnames(LOW), ignore.case = TRUE)]]
      } else {
        LOW[[2]]
      }
      
      # Smooth and compute slopes
      ma5 <- stats::filter(lowy, rep(1/5, 5), sides = 2)
      ma5[is.na(ma5)] <- lowy[is.na(ma5)]
      svec <- diff(as.numeric(ma5)) / diff(as.numeric(lowx))
      
      # Find where slope transitions from negative to non-negative (flattening)
      neg_to_flat <- which(diff(sign(svec)) > 0)
      if (length(neg_to_flat)) {
        Oc <- as.numeric(lowx[neg_to_flat[1] + 1])
      } else {
        Oc <- max(lowx, na.rm = TRUE)  # fallback if always negative
      }
    }
    
    # fNS: prefer median of FNS$Fns
    fNS_median <- if (!is.null(FNS)) {
      fn_col <- grep("Fn|fn|FNS|Fns", colnames(FNS), value = TRUE)[1]
      if (!is.na(fn_col)) median(FNS[[fn_col]], na.rm = TRUE) else median(unlist(FNS[, sapply(FNS, is.numeric), drop = FALSE]), na.rm = TRUE)
    } else NA_real_
    
    # slopes: prefer LME$Slope
    slope_median <- NA_real_; p_val <- NA_real_; slopes <- NULL
    if (!is.null(LME)) {
      slope_col <- grep("Slope|slope|beta", colnames(LME), value = TRUE)[1]
      if (!is.na(slope_col)) {
        slopes <- as.numeric(LME[[slope_col]])
        slope_median <- median(slopes, na.rm = TRUE)
        p_val <- mean(slopes >= 0, na.rm = TRUE)  # one-tailed per paper
      } else {
        numeric_cols <- unlist(LME[, sapply(LME, is.numeric), drop = FALSE])
        slopes <- as.numeric(numeric_cols)
        slope_median <- median(slopes, na.rm = TRUE)
        p_val <- mean(slopes >= 0, na.rm = TRUE)
      }
    }
    
    # craft CI frame: ensure column names workable
    CI_frame <- NULL
    if (!is.null(CI)) {
      cis <- find_CI_cols(CI)
      CI_frame <- data.frame(Overlap = CI[[cis$overlap]],
                             ymin = if (!is.na(cis$lower)) CI[[cis$lower]] else NA,
                             y = if (!is.na(cis$median)) CI[[cis$median]] else NA,
                             ymax = if (!is.na(cis$upper)) CI[[cis$upper]] else NA)
    } else if (!is.null(LOW)) {
      # create CI-like frame from LOWESS line if CI missing
      CI_frame <- data.frame(Overlap = LOW[[1]],
                             ymin = NA_real_, y = if ("LOWESS" %in% colnames(LOW)) LOW$LOWESS else LOW[[2]],
                             ymax = NA_real_)
    } else {
      CI_frame <- data.frame(Overlap = DO$Overlap, ymin = NA_real_, y = NA_real_, ymax = NA_real_)
    }
    
    list(DO = DO, CI = CI_frame, Oc = Oc, fNS_median = fNS_median,
         slope_median = slope_median, p_value = p_val, slopes = slopes)
  }
  
  m1 <- extract_metrics_from_DOC(doc1)
  m2 <- extract_metrics_from_DOC(doc2)
  
  # tag datasets
  m1$DO$Dataset <- title1
  m2$DO$Dataset <- title2
  m1$CI$Dataset <- title1
  m2$CI$Dataset <- title2
  
  combined_DO <- bind_rows(m1$DO, m2$DO)
  combined_CI <- bind_rows(m1$CI, m2$CI)
  
  # densities scaled for x-axis
  y_max <- max(combined_DO$rJSD, na.rm = TRUE)
  density_scale <- if (is.finite(y_max) && y_max > 0) 0.15 * y_max else 0.15
  densities <- combined_DO %>%
    group_by(Dataset) %>%
    do({
      od <- .$Overlap
      if (length(od) < 2 || all(is.na(od))) return(data.frame(x = numeric(0), y = numeric(0)))
      dens <- stats::density(od, from = min(combined_DO$Overlap, na.rm = TRUE),
                             to = max(combined_DO$Overlap, na.rm = TRUE))
      data.frame(x = dens$x, y = dens$y / max(dens$y) * density_scale)
    }) %>% ungroup()
  
  # colors (user can override)
  col <- c(Harvest = "#F8766D", Storage = "#00BFC4")
  if (!title1 %in% names(col)) col[title1] <- "#F8766D"
  if (!title2 %in% names(col)) col[title2] <- "#00BFC4"
  
  # plotting using known CI column names (y, ymin, ymax)
  p <- ggplot() +
    geom_ribbon(data = combined_CI, aes(x = Overlap, ymin = ymin, ymax = ymax, fill = Dataset), alpha = 0.25) +
    geom_line(data = combined_CI, aes(x = Overlap, y = y, color = Dataset), size = 1) +
    geom_point(data = combined_DO, aes(x = Overlap, y = rJSD, color = Dataset), size = 0.8, alpha = 0.6) +
    geom_vline(xintercept = m1$Oc, color = col[title1], linetype = "dashed") +
    geom_vline(xintercept = m2$Oc, color = col[title2], linetype = "dashed") +
    geom_area(data = densities, aes(x = x, y = y, fill = Dataset), alpha = 0.3, position = "identity") +
    scale_color_manual(values = col, name = "Sampling") +  
    scale_fill_manual(values = col, name = "Sampling") +    
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    xlab("Overlap") + ylab("Dissimilarity (rJSD)") +
    coord_cartesian(ylim = c(0, 1))
  
  # annotation text (robust formatting / NA safe)
  fmt <- function(x, d = 2) if (is.na(x)) "NA" else formatC(x, digits = d, format = "f")
  summary_text <- paste0(
    sprintf("%s: Oc=%s | fNS=%s |  p=%s",
            title1, fmt(m1$Oc), fmt(m1$fNS_median), fmt(m1$p_value)),
    "\n",
    sprintf("%s: Oc=%s | fNS=%s |  p=%s",
            title2, fmt(m2$Oc), fmt(m2$fNS_median), fmt(m2$p_value))
  )
  
  p <- p + annotate("text",
                    x = min(combined_DO$Overlap, na.rm = TRUE) + 0.02 * diff(range(combined_DO$Overlap, na.rm = TRUE)),
                    y = 0.98, label = summary_text, hjust = 0, vjust = 1,
                    size = 4, fontface = "bold", lineheight = 1.05)
  
  print(p)
  
  return(list(metrics1 = m1, metrics2 = m2, plot = p))
}

# assuming doc_Harvest and doc_Storage exist
results_combined <- analyze_DOC_combined_fromDOC(doc_Harvest, doc_Storage,
                                                 title1 = "Harvest", title2 = "Storage")
results <- results_combined$plot + theme_nature
results
saveRDS(results, "doc_its.rds")

plot_DOC_single <- function(doc_obj, title = "DOC Analysis") {
  # --- Extract components safely ---
  DO <- if (!is.null(doc_obj$DO)) as.data.frame(doc_obj$DO) else stop("DOC object has no DO component")
  CI <- if (!is.null(doc_obj$CI)) {
    CI_df <- as.data.frame(doc_obj$CI)
    # Ensure standard column names for plotting
    colnames(CI_df) <- c("Overlap", "ymin", "y", "ymax")[1:ncol(CI_df)]
    CI_df
  } else {
    # fallback if CI missing
    data.frame(Overlap = DO$Overlap, ymin = NA, y = NA, ymax = NA)
  }
  
  # --- Extract Oc and fNS ---
  Oc <- if (!is.null(doc_obj$NEG)) {
    neg_col <- grep("Neg|neg|Oc|Neg.Slope", colnames(doc_obj$NEG), value = TRUE)[1]
    if (!is.na(neg_col)) median(doc_obj$NEG[[neg_col]], na.rm = TRUE) else NA
  } else if (!is.null(doc_obj$LOWESS)) {
    lowy <- doc_obj$LOWESS[,1]
    Oc <- median(lowy, na.rm = TRUE)
  } else NA_real_
  
  fNS <- if (!is.null(doc_obj$FNS)) {
    fn_col <- grep("Fn|FNS|Fns", colnames(doc_obj$FNS), value = TRUE)[1]
    if (!is.na(fn_col)) median(doc_obj$FNS[[fn_col]], na.rm = TRUE) else NA
  } else NA_real_
  
  # --- Density for visual context ---
  density_scale <- 0.15 * max(DO$rJSD, na.rm = TRUE)
  dens <- if (length(DO$Overlap) > 1) {
    dens_tmp <- density(DO$Overlap, from = min(DO$Overlap, na.rm = TRUE),
                        to = max(DO$Overlap, na.rm = TRUE))
    data.frame(x = dens_tmp$x, y = dens_tmp$y / max(dens_tmp$y) * density_scale)
  } else data.frame(x = numeric(0), y = numeric(0))
  
  # --- Plot ---
  p <- ggplot() +
    geom_ribbon(data = CI, aes(x = Overlap, ymin = ymin, ymax = ymax), alpha = 0.25, fill = "steelblue") +
    geom_line(data = CI, aes(x = Overlap, y = y), color = "steelblue", size = 1) +
    geom_point(data = DO, aes(x = Overlap, y = rJSD), size = 1, alpha = 0.6, color = "darkred") +
    geom_vline(xintercept = Oc, linetype = "dashed", color = "darkgreen") +
    geom_area(data = dens, aes(x = x, y = y), alpha = 0.3, fill = "grey") +
    theme_bw() +
    xlab("Overlap") +
    ylab("Dissimilarity (rJSD)") +
    ggtitle(title) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))
  
  # --- Annotate Oc and fNS ---
  summary_text <- paste0("Oc = ", round(Oc, 3), " | fNS = ", round(fNS, 3))
  p <- p + annotate("text",
                    x = min(DO$Overlap, na.rm = TRUE) + 0.02 * diff(range(DO$Overlap, na.rm = TRUE)),
                    y = max(DO$rJSD, na.rm = TRUE) * 0.95,
                    label = summary_text,
                    hjust = 0, vjust = 1,
                    size = 4, fontface = "bold")
  
  print(p)
  return(p)
}

plot_DOC_single(doc, title = "Fungi")

# ---- (1) package-level summaries ----
summary_pkg <- list(
  Oc_median_pkg  = median(doc_Harvest$NEG$Neg.Slope, na.rm = TRUE),
  fNS_median_pkg = median(doc_Harvest$FNS$Fns, na.rm = TRUE),
  slope_median_pkg = median(doc_Harvest$LME$Slope, na.rm = TRUE),
  p_value_pkg = mean(doc_Harvest$LME$Slope >= 0, na.rm = TRUE)
)

summary_pkg

# ---- (2) metrics from your combined analysis ----
summary_new <- with(results_combined$metrics1,
                    c(Oc = Oc, fNS = fNS_median,
                      slope_median = slope_median, p = p_value))
summary_new

# ---- (3) compare numerically ----
rbind(Package = summary_pkg, New_Function = summary_new)

# ---- (1) package-level summaries ----
summary_pkg_Storage <- list(
  Oc_median_pkg  = median(doc_Storage$NEG$Neg.Slope, na.rm = TRUE),
  fNS_median_pkg = median(doc_Storage$FNS$Fns, na.rm = TRUE),
  slope_median_pkg = median(doc_Storage$LME$Slope, na.rm = TRUE),
  p_value_pkg = mean(doc_Storage$LME$Slope >= 0, na.rm = TRUE)
)

summary_pkg_Storage

# ---- (2) metrics from your combined analysis ----
summary_new_Storage <- with(results_combined$metrics2,
                            c(Oc = Oc, fNS = fNS_median,
                              slope_median = slope_median, p = p_value))
summary_new

# ---- (3) compare numerically ----
rbind(Package = summary_pkg_Storage, New_Function = summary_new_Storage)

# Bootstrap check ---------------------------------------------------------

par(mfrow = c(1,3))
hist(doc_Harvest$NEG$Neg.Slope, main = "Bootstrap Oc (Harvest)",
     xlab = "Oc", col = "skyblue")
hist(doc_Harvest$FNS$Fns, main = "Bootstrap fNS (Harvest)",
     xlab = "fNS", col = "lightgreen")
hist(doc_Harvest$LME$Slope, main = "Bootstrap slopes (Harvest)",
     xlab = "Slope", col = "lightcoral")

par(mfrow = c(1,3))
hist(doc_Storage$NEG$Neg.Slope, main = "Bootstrap Oc (Storage)",
     xlab = "Oc", col = "skyblue")
hist(doc_Storage$FNS$Fns, main = "Bootstrap fNS (Storage)",
     xlab = "fNS", col = "lightgreen")
hist(doc_Storage$LME$Slope, main = "Bootstrap slopes (Storage)",
     xlab = "Slope", col = "lightcoral")

head(results_combined$metrics1$CI)
quantile(doc_Harvest$BOOT$rJSD.Boot.1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
head(results_combined$metrics2$CI)
quantile(doc_Storage$BOOT$rJSD.Boot.1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)

compare_DOC_distributions <- function(doc1, doc2,
                                      title1 = "Harvest", title2 = "Storage",
                                      method = "wilcox.test",
                                      p.adjust.method = "BH") {
  # helper to safely extract numeric bootstrap columns
  safe_extract <- function(df, pattern) {
    if (is.null(df)) return(NA_real_)
    col <- grep(pattern, colnames(df), value = TRUE)[1]
    if (is.na(col)) {
      nums <- df[, sapply(df, is.numeric), drop = FALSE]
      return(unlist(nums))
    } else {
      return(as.numeric(df[[col]]))
    }
  }
  
  # extract numeric vectors
  Oc1 <- safe_extract(doc1$NEG, "Neg|Oc|Neg.Slope")
  Oc2 <- safe_extract(doc2$NEG, "Neg|Oc|Neg.Slope")
  
  fNS1 <- safe_extract(doc1$FNS, "Fn|FNS")
  fNS2 <- safe_extract(doc2$FNS, "Fn|FNS")
  
  slope1 <- safe_extract(doc1$LME, "Slope|slope|beta")
  slope2 <- safe_extract(doc2$LME, "Slope|slope|beta")
  
  # define comparison function
  cmp <- function(x, y, method) {
    if (all(is.na(x)) || all(is.na(y))) return(NA_real_)
    x <- x[is.finite(x)]
    y <- y[is.finite(y)]
    if (length(x) < 5 || length(y) < 5) return(NA_real_)
    if (method == "wilcox.test") {
      res <- try(stats::wilcox.test(x, y)$p.value, silent = TRUE)
      if (inherits(res, "try-error")) return(NA_real_) else return(res)
    } else if (method == "t.test") {
      res <- try(stats::t.test(x, y)$p.value, silent = TRUE)
      if (inherits(res, "try-error")) return(NA_real_) else return(res)
    } else {
      stop("Unsupported test method")
    }
  }
  
  # compute raw p-values
  p_Oc <- cmp(Oc1, Oc2, method)
  p_fNS <- cmp(fNS1, fNS2, method)
  p_slope <- cmp(slope1, slope2, method)
  
  # multiple-testing correction
  pvals <- c(p_Oc, p_fNS, p_slope)
  padj <- stats::p.adjust(pvals, method = p.adjust.method)
  
  # summary table
  results <- data.frame(
    Metric = c("Oc", "fNS", "Slope"),
    p_raw = signif(pvals, 3),
    p_adj = signif(padj, 3),
    Median_1 = c(median(Oc1, na.rm = TRUE),
                 median(fNS1, na.rm = TRUE),
                 median(slope1, na.rm = TRUE)),
    Median_2 = c(median(Oc2, na.rm = TRUE),
                 median(fNS2, na.rm = TRUE),
                 median(slope2, na.rm = TRUE)),
    Dataset_1 = title1,
    Dataset_2 = title2
  )
  
  return(results)
}

comparison_results <- compare_DOC_distributions(doc_Harvest, doc_Storage,
                                                title1 = "Harvest",
                                                title2 = "Storage",
                                                method = "wilcox.test",
                                                p.adjust.method = "BH")

print(comparison_results)

results_combined <- comparison_results$plot + theme_nature
results_combined
saveRDS(results_combined, "its_doc.rds")

analyze_DOC_multi_fromDOC <- function(doc_list, 
                                      titles = c("Geoxe", "Control", "Ulmasud"),
                                      phase = "Harvest") {
  
  # ---- helper: robust CI column finder ----
  find_CI_cols <- function(CI) {
    cn <- colnames(CI)
    overlap_col <- if ("Overlap" %in% cn) "Overlap" else cn[1]
    lower_col <- cn[grep("2\\.5|2_5|2p5|X2\\.5", cn)[1]]
    median_col <- cn[grep("(^|[^0-9])50([^0-9]|$)|X50|50\\.", cn)[1]]
    upper_col <- cn[grep("97\\.5|97_5|97p5|X97\\.5", cn)[1]]
    if (is.na(lower_col) && length(cn) >= 2) lower_col <- cn[2]
    if (is.na(median_col) && length(cn) >= 3) median_col <- cn[3]
    if (is.na(upper_col) && length(cn) >= 4) upper_col <- cn[4]
    list(overlap = overlap_col, lower = lower_col, median = median_col, upper = upper_col)
  }
  
  # ---- helper: metric extraction (exactly as your original) ----
  extract_metrics_from_DOC <- function(doc_obj) {
    DO <- as.data.frame(doc_obj$DO)
    CI <- if (!is.null(doc_obj$CI)) as.data.frame(doc_obj$CI) else NULL
    LOW <- if (!is.null(doc_obj$LOWESS)) as.data.frame(doc_obj$LOWESS) else NULL
    LME  <- if (!is.null(doc_obj$LME)) as.data.frame(doc_obj$LME) else NULL
    NEG  <- if (!is.null(doc_obj$NEG)) as.data.frame(doc_obj$NEG) else NULL
    FNS  <- if (!is.null(doc_obj$FNS)) as.data.frame(doc_obj$FNS) else NULL
    
    # Oc: prefer NEG$Neg.Slope, fallback to LOWESS/LOWY
    Oc <- NA_real_
    if (!is.null(NEG)) {
      neg_col <- grep("Neg|neg|Oc|Oc\\.?|Neg.Slope", colnames(NEG), value = TRUE)[1]
      if (!is.na(neg_col)) Oc <- median(NEG[[neg_col]], na.rm = TRUE)
      else Oc <- median(unlist(NEG[, sapply(NEG, is.numeric), drop = FALSE]), na.rm = TRUE)
    } else if (!is.null(LOW)) {
      lowx <- LOW[[1]]
      lowy <- if (any(grepl("LOWY", colnames(LOW), ignore.case = TRUE))) {
        LOW[[grep("LOWY", colnames(LOW), ignore.case = TRUE)]]
      } else if (any(grepl("LOWESS", colnames(LOW), ignore.case = TRUE))) {
        LOW[[grep("LOWESS", colnames(LOW), ignore.case = TRUE)]]
      } else LOW[[2]]
      ma5 <- stats::filter(lowy, rep(1/5, 5), sides = 2)
      ma5[is.na(ma5)] <- lowy[is.na(ma5)]
      svec <- diff(as.numeric(ma5)) / diff(as.numeric(lowx))
      neg_to_flat <- which(diff(sign(svec)) > 0)
      if (length(neg_to_flat)) Oc <- as.numeric(lowx[neg_to_flat[1] + 1])
      else Oc <- max(lowx, na.rm = TRUE)
    }
    
    # fNS and slope metrics
    fNS_median <- if (!is.null(FNS)) {
      fn_col <- grep("Fn|fn|FNS|Fns", colnames(FNS), value = TRUE)[1]
      if (!is.na(fn_col)) median(FNS[[fn_col]], na.rm = TRUE)
      else median(unlist(FNS[, sapply(FNS, is.numeric), drop = FALSE]), na.rm = TRUE)
    } else NA_real_
    
    slope_median <- NA_real_; p_val <- NA_real_; slopes <- NULL
    if (!is.null(LME)) {
      slope_col <- grep("Slope|slope|beta", colnames(LME), value = TRUE)[1]
      if (!is.na(slope_col)) {
        slopes <- as.numeric(LME[[slope_col]])
        slope_median <- median(slopes, na.rm = TRUE)
        p_val <- mean(slopes >= 0, na.rm = TRUE)
      } else {
        numeric_cols <- unlist(LME[, sapply(LME, is.numeric), drop = FALSE])
        slopes <- as.numeric(numeric_cols)
        slope_median <- median(slopes, na.rm = TRUE)
        p_val <- mean(slopes >= 0, na.rm = TRUE)
      }
    }
    
    # Craft CI
    CI_frame <- NULL
    if (!is.null(CI)) {
      cis <- find_CI_cols(CI)
      CI_frame <- data.frame(Overlap = CI[[cis$overlap]],
                             ymin = if (!is.na(cis$lower)) CI[[cis$lower]] else NA,
                             y = if (!is.na(cis$median)) CI[[cis$median]] else NA,
                             ymax = if (!is.na(cis$upper)) CI[[cis$upper]] else NA)
    } else if (!is.null(LOW)) {
      CI_frame <- data.frame(Overlap = LOW[[1]],
                             ymin = NA_real_,
                             y = if ("LOWESS" %in% colnames(LOW)) LOW$LOWESS else LOW[[2]],
                             ymax = NA_real_)
    } else {
      CI_frame <- data.frame(Overlap = DO$Overlap, ymin = NA_real_, y = NA_real_, ymax = NA_real_)
    }
    
    list(DO = DO, CI = CI_frame, Oc = Oc, fNS_median = fNS_median,
         slope_median = slope_median, p_value = p_val, slopes = slopes)
  }
  
  # ---- Extract for all treatments ----
  metrics <- lapply(doc_list, extract_metrics_from_DOC)
  names(metrics) <- titles
  
  for (i in seq_along(titles)) {
    metrics[[i]]$DO$Dataset <- titles[i]
    metrics[[i]]$CI$Dataset <- titles[i]
  }
  
  combined_DO <- dplyr::bind_rows(lapply(metrics, `[[`, "DO"))
  combined_CI <- dplyr::bind_rows(lapply(metrics, `[[`, "CI"))
  
  # ---- Density (same logic as yours) ----
  y_max <- max(combined_DO$rJSD, na.rm = TRUE)
  density_scale <- if (is.finite(y_max) && y_max > 0) 0.15 * y_max else 0.15
  densities <- combined_DO %>%
    group_by(Dataset) %>%
    do({
      od <- .$Overlap
      if (length(od) < 2 || all(is.na(od))) return(data.frame(x = numeric(0), y = numeric(0)))
      dens <- stats::density(od,
                             from = min(combined_DO$Overlap, na.rm = TRUE),
                             to = max(combined_DO$Overlap, na.rm = TRUE))
      data.frame(x = dens$x, y = dens$y / max(dens$y) * density_scale)
    }) %>% ungroup()
  
  # ---- Color palette ----
  col <- c(Control = "#80b1d3", Geoxe = "#fb8072", Ulmasud = "#b3de69")
  
  # ---- Plot ----
  p <- ggplot() +
    geom_ribbon(data = combined_CI,
                aes(x = Overlap, ymin = ymin, ymax = ymax, fill = Dataset),
                alpha = 0.25) +
    geom_line(data = combined_CI,
              aes(x = Overlap, y = y, color = Dataset), size = 1) +
    geom_point(data = combined_DO,
               aes(x = Overlap, y = rJSD, color = Dataset), size = 0.8, alpha = 0.6) +
    geom_area(data = densities, aes(x = x, y = y, fill = Dataset),
              alpha = 0.3, position = "identity") +
    scale_color_manual(values = col, name = "Treatment") +
    scale_fill_manual(values = col, name = "Treatment") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top") +
    xlab("Overlap") + ylab("Dissimilarity (rJSD)") +
    coord_cartesian(ylim = c(0, 1)) +
    ggtitle(paste("DOC curves -", phase))
  
  # ---- Annotate summary ----
  fmt <- function(x, d = 2) if (is.na(x)) "NA" else formatC(x, digits = d, format = "f")
  summary_lines <- sapply(seq_along(titles), function(i) {
    m <- metrics[[i]]
    sprintf("%s: Oc=%s | fNS=%s | p=%s",
            titles[i], fmt(m$Oc), fmt(m$fNS_median), fmt(m$p_value))
  })
  summary_text <- paste(summary_lines, collapse = "\n")
  
  p <- p + annotate("text",
                    x = min(combined_DO$Overlap, na.rm = TRUE) + 
                      0.02 * diff(range(combined_DO$Overlap, na.rm = TRUE)),
                    y = 0.98, label = summary_text,
                    hjust = 0, vjust = 1, size = 4,
                    fontface = "bold", lineheight = 1.05)
  
  print(p)
  return(list(metrics = metrics, plot = p))
}


# Harvest phase
results_Harvest <- analyze_DOC_multi_fromDOC(
  doc_list = list(doc_Harvest_none, doc_Harvest_ulmasud, doc_Harvest_geoxe),
  titles = c("Control", "Geoxe", "Ulmasud"),
  phase = "Harvest"
)

# Storage phase
results_Storage <- analyze_DOC_multi_fromDOC(
  doc_list = list(doc_Storage_none, doc_Storage_ulmasud, doc_Storage_geoxe),
  titles = c("Control", "Geoxe", "Ulmasud"),
  phase = "Storage"
)
results_harvest <- results_Harvest$plot + theme_nature
results_storage <- results_Storage$plot + theme_nature
its_combined <- wrap_plots(results_harvest, results_storage, ncol = 2, guides = "collect") &
  theme(legend.position = "top")  # single legend

# Remove any internal titles to avoid conflicts
its_combined <- its_combined 

its_combined

saveRDS(its_combined, "DOC_Harvest_Storage_combined.rds")
saveRDS(results_harvest, "DOC_Harvest_treatments.rds")
saveRDS(results_storage, "DOC_Storage_treatments.rds")






