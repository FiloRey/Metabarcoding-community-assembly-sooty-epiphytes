# libraries ----------------------------------------------------------------
library(gridGraphics)
library(cowplot)
library(readxl)
library(openxlsx)
library(gt)
library(qiime2R)
library(tidyverse)      # For data manipulation and visualization
library(phyloseq)       # For handling phylogenetic data
library(ggpubr)         # For arranging ggplot objects
library(microbiome)     # For microbiome analysis
library(vegan)          # For community ecology analyses
library(cluster)        # For clustering methods
library(dplyr)         # Already included in tidyverse, but kept for clarity
library(readr)         # Already included in tidyverse, but kept for clarity
library(knitr)         
library(ggpmisc)
library(kableExtra) 
library(webshot2)       # For rendering markdown tables
library(magrittr)       # For the pipe operator
library(DESeq2)         # For differential abundance analysis
library(DirichletMultinomial) # For analyzing driving ASVs in community structures
library(RColorBrewer)   # For color palettes
library(reshape2) # For handling matrices and other structures
library(gplots)        # For standard plots
library(gridExtra) 
library(ggfortify) # For arranging plots
library(ggforce)       # For adding polygons to ggplots
library(ggrepel)       # For repelling text labels in ggplots
library(rstatix)       # For statistical analysis
library(FSA)           # For fisheries stock assessment
library(dunn.test)     # For Dunn's test
library(indicspecies)  # For indicator species analysis
library(openxlsx)
library(ggplot2)       # Already included in tidyverse, but kept for clarity
library(pheatmap)
library(BiocParallel)
library(mia)
library(miaViz)
library(scater)
library(grid)
library(sechm)
library(patchwork)
library(ComplexHeatmap)
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
library(seqtime)  
library(maaslin3)
library(qgraph)
library(shadowtext)
library(SPRING)
library(NetCoMi)
library(igraph)
library(circlize)
library(RColorBrewer)
library(corpcor)
library(SpiecEasi)
library(DT)
library(picante)
library(lme4)         # Mixed model
library(lmerTest)     # P-values for lmer
library(emmeans)      # Estimated marginal means and pairwise tests
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

# ===============================
# COMMON PUBLICATION THEME -------------------------------------------------
# ===============================
# Base font size and family
base_text_size <- 25
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

# INTERKINGDOM NETWORKS NETWORKS ---------------------------------------------------------
# NETWORK COMPARISON ---------------------------------------------
physeq_healthy_merged <- NetCoMi::renameTaxa(physeq_healthy_merged)
physeq_infected_merged <- NetCoMi::renameTaxa(physeq_infected_merged)
sum(duplicated(tax_table(physeq_healthy_merged)[, "Genus"]))
sum(duplicated(tax_table(physeq_infected_merged)[, "Genus"]))
dup_healthy <- tax_table(physeq_healthy_merged)[duplicated(tax_table(physeq_healthy_merged)[, "Genus"]), "Genus"]
dup_infected <- tax_table(physeq_infected_merged)[duplicated(tax_table(physeq_infected_merged)[, "Genus"]), "Genus"]

dup_healthy
dup_infected

tax_table(physeq_healthy_merged)[, "Genus"] <- make.unique(tax_table(physeq_healthy_merged)[, "Genus"])
tax_table(physeq_infected_merged)[, "Genus"] <- make.unique(tax_table(physeq_infected_merged)[, "Genus"])

view(otu_table(physeq_healthy_merged))

spring_net <- netConstruct(
  data = physeq_healthy_merged,
  data2 = physeq_infected_merged,
  
  #dataType = "relabundance",
  taxRank = "Genus",
  
  # Preprocessing:
  filtSamp = "none",
  
  zeroMethod = "none",
  normMethod = "none",
  
  filtTax = "highestVar",
  filtTaxPar = list(highestVar  = 80),
  measure = "spring",
  measurePar = list(
    nlambda = 20,
    rep.num = 20,
    thresh = 0.05,
    Rmethod = "approx"),
  sparsMethod = "none",
  dissFunc = "signed",
  verbose = 3,
  seed = 13075)

spring_netprop <- netAnalyze(
  spring_net,
  centrLCC = FALSE,
  avDissIgnoreInf = TRUE,
  sPathNorm = FALSE,
  clustMethod = "cluster_fast_greedy",
  hubPar = c("degree", "betwenness"),
  hubQuant = 0.9,
  lnormFit = TRUE,
  normDeg = FALSE,
  normBetw = FALSE,
  normClose = FALSE,
  normEigen = T)

summary(spring_netprop, groupNames = c("Healthy", "Infected"))
sum_text <- capture.output(
  summary(
    spring_netprop,
    groupNames = c("Healthy", "Infected")
  )
)

sum_df <- data.frame(Output = sum_text)

library(writexl)
write_xlsx(sum_df, "spring_netprop_summary.xlsx")

# Healthy edges
healthy_edges <- spring_net$edgelist1
colnames(healthy_edges) <- c("source", "target", "weight")

write.csv(healthy_edges, "healthy_edges.csv", row.names = FALSE)

# Infected edges
infected_edges <- spring_net$edgelist2
colnames(infected_edges) <- c("source", "target", "weight")

write.csv(infected_edges, "infected_edges.csv", row.names = FALSE)

nodes_healthy <- unique(c(
  healthy_edges$source,
  healthy_edges$target
))

healthy_nodes <- data.frame(
  name = nodes_healthy
)

cent_h <- spring_netprop$centralities

healthy_nodes$degree  <- cent_h$degree1[healthy_nodes$name]
healthy_nodes$between <- cent_h$between1[healthy_nodes$name]
healthy_nodes$close   <- cent_h$close1[healthy_nodes$name]
healthy_nodes$eigenv  <- cent_h$eigenv1[healthy_nodes$name]

healthy_nodes$hub <- healthy_nodes$name %in% spring_netprop$hubs[[1]]

healthy_nodes$name <- gsub('"', '', healthy_nodes$name)
healthy_edges$source <- gsub('"', '', healthy_edges$source)
healthy_edges$target <- gsub('"', '', healthy_edges$target)

write.csv(healthy_nodes, "healthy_nodes.csv",
          row.names = FALSE, quote = FALSE)

write.csv(healthy_edges, "healthy_edges.csv",
          row.names = FALSE, quote = FALSE)



nodes_infected <- unique(c(
  infected_edges$source,
  infected_edges$target
))

infected_nodes <- data.frame(
  name = nodes_infected
)

cent_h <- spring_netprop$centralities

infected_nodes$degree  <- cent_h$degree2[infected_nodes$name]
infected_nodes$between <- cent_h$between2[infected_nodes$name]
infected_nodes$close   <- cent_h$close2[infected_nodes$name]
infected_nodes$eigenv  <- cent_h$eigenv2[infected_nodes$name]

infected_nodes$hub <- infected_nodes$name %in% spring_netprop$hubs[[1]]

infected_nodes$name <- gsub('"', '', infected_nodes$name)
infected_edges$source <- gsub('"', '', infected_edges$source)
infected_edges$target <- gsub('"', '', infected_edges$target)

write.csv(infected_nodes, "infected_nodes.csv",
          row.names = FALSE, quote = FALSE)

write.csv(infected_edges, "infected_edges.csv",
          row.names = FALSE, quote = FALSE)


# NETWORK COMPARISON ---------------------------------------------
spring_net_fun <- netConstruct(
  data = physeq_fun_healthy,
  data2 = physeq_fun_infected,
  
  #dataType = "relabundance",
  taxRank = "Genus",
  
  # Preprocessing:
  filtSamp = "none",
  
  zeroMethod = "none",
  normMethod = "none",
  
  filtTax = "highestVar",
  filtTaxPar = list(highestVar  = 80),
  measure = "spring",
  measurePar = list(
    nlambda = 20,
    rep.num = 20,
    thresh = 0.05,
    Rmethod = "approx"),
  sparsMethod = "none",
  dissFunc = "signed",
  verbose = 3,
  seed = 13075)


 spring_netprop_fun <- netAnalyze(
  spring_net_fun,
  centrLCC = FALSE,
  avDissIgnoreInf = TRUE,
  sPathNorm = FALSE,
  clustMethod = "cluster_fast_greedy",
  hubPar = c("degree", "betwenness"),
  hubQuant = 0.9,
  lnormFit = TRUE,
  normDeg = FALSE,
  normBetw = FALSE,
  normClose = FALSE,
  normEigen = T)

summary(spring_netprop_fun, groupNames = c("Healthy", "Infected"))
sum_text <- capture.output(
  summary(
    spring_netprop_fun,
    groupNames = c("Healthy", "Infected")
  )
)

sum_df_fun <- data.frame(Output = sum_text)

library(writexl)
write_xlsx(sum_df_fun, "spring_netprop_summary_fun.xlsx")

# Healthy edges
healthy_edges_fun <- spring_net_fun$edgelist1
colnames(healthy_edges_fun) <- c("source", "target", "weight")

write.csv(healthy_edges, "healthy_edges_fun.csv", row.names = FALSE)

# Infected edges
infected_edges_fun <- spring_net_fun$edgelist2
colnames(infected_edges_fun) <- c("source", "target", "weight")

write.csv(infected_edges_fun, "infected_edges_fun.csv", row.names = FALSE)

nodes_healthy_fun <- unique(c(
  healthy_edges_fun$source,
  healthy_edges_fun$target
))

healthy_nodes_fun <- data.frame(
  name = nodes_healthy_fun
)

cent_h_fun <- spring_netprop_fun$centralities

healthy_nodes_fun$degree  <- cent_h_fun$degree1[healthy_nodes_fun$name]
healthy_nodes_fun$between <- cent_h_fun$between1[healthy_nodes_fun$name]
healthy_nodes_fun$close   <- cent_h_fun$close1[healthy_nodes_fun$name]
healthy_nodes_fun$eigenv  <- cent_h_fun$eigenv1[healthy_nodes_fun$name]

healthy_nodes_fun$hub <- healthy_nodes_fun$name %in% spring_netprop_fun$hubs[[1]]

healthy_nodes_fun$name <- gsub('"', '', healthy_nodes_fun$name)
healthy_edges_fun$source <- gsub('"', '', healthy_edges_fun$source)
healthy_edges_fun$target <- gsub('"', '', healthy_edges_fun$target)

write.csv(healthy_nodes_fun, "healthy_nodes_fun.csv",
          row.names = FALSE, quote = FALSE)

write.csv(healthy_edges_fun, "healthy_edges_fun.csv",
          row.names = FALSE, quote = FALSE)



nodes_infected_fun <- unique(c(
  infected_edges_fun$source,
  infected_edges_fun$target
))

infected_nodes_fun <- data.frame(
  name = nodes_infected_fun
)

cent_h_fun <- spring_netprop_fun$centralities
view(cent_h_fun)
infected_nodes_fun$degree  <- cent_h_fun$degree2[infected_nodes_fun$name]
infected_nodes_fun$between <- cent_h_fun$between2[infected_nodes_fun$name]
infected_nodes_fun$close   <- cent_h_fun$close2[infected_nodes_fun$name]
infected_nodes_fun$eigenv  <- cent_h_fun$eigenv2[infected_nodes_fun$name]

infected_nodes_fun$hub <- infected_nodes_fun$name %in% spring_netprop_fun$hubs[[1]]

infected_nodes_fun$name <- gsub('"', '', infected_nodes_fun$name)
infected_edges_fun$source <- gsub('"', '', infected_edges_fun$source)
infected_edges_fun$target <- gsub('"', '', infected_edges_fun$target)

write.csv(infected_nodes_fun, "infected_nodes_fun.csv",
          row.names = FALSE, quote = FALSE)

write.csv(infected_edges_fun, "infected_edges_fun.csv",
          row.names = FALSE, quote = FALSE)



sum(duplicated(tax_table(physeq_bac_healthy)[, "Genus"]))
sum(duplicated(tax_table(physeq_bac_infected)[, "Genus"]))
dup_healthy <- tax_table(physeq_bac_healthy)[duplicated(tax_table(physeq_bac_healthy)[, "Genus"]), "Genus"]
dup_infected <- tax_table(physeq_bac_infected)[duplicated(tax_table(physeq_bac_infected)[, "Genus"]), "Genus"]

dup_healthy
dup_infected

tax_table(physeq_bac_healthy)[, "Genus"] <- make.unique(tax_table(physeq_bac_healthy)[, "Genus"])
tax_table(physeq_bac_infected)[, "Genus"] <- make.unique(tax_table(physeq_bac_infected)[, "Genus"])



spring_net_bac <- netConstruct(
  data = physeq_bac_healthy,
  data2 = physeq_bac_infected,
  #dataType = "relabundance",
  taxRank = "Genus",
  
  # Preprocessing:
  filtSamp = "none",
  
  zeroMethod = "none",
  normMethod = "none",
  
  filtTax = "highestVar",
  filtTaxPar = list(highestVar  = 80),
  measure = "spring",
  measurePar = list(
    nlambda = 20,
    rep.num = 20,
    thresh = 0.05,
    Rmethod = "approx"),
  sparsMethod = "none",
  dissFunc = "signed",
  verbose = 3,
  seed = 13075

)

spring_netprop_bac <- netAnalyze(
  spring_net_bac,
  centrLCC = FALSE,
  avDissIgnoreInf = TRUE,
  sPathNorm = FALSE,
  clustMethod = "cluster_fast_greedy",
  hubPar = c("degree", "betwenness", "eigenvals"),
  hubQuant = 0.9,
  lnormFit = TRUE,
  normDeg = FALSE,
  normBetw = FALSE,
  normClose = FALSE,
  normEigen = T)

summary(spring_netprop_bac, groupNames = c("Healthy", "Infected"))
sum_text_bac <- capture.output(
  summary(
    spring_netprop_bac,
    groupNames = c("Healthy", "Infected")
  )
)

sum_df_bac <- data.frame(Output = sum_text_bac)

library(writexl)
write_xlsx(sum_df_bac, "spring_netprop_summary_bac.xlsx")

# Healthy edges
healthy_edges_bac <- spring_net_bac$edgelist1
colnames(healthy_edges_bac) <- c("source", "target", "weight")

write.csv(healthy_edges, "healthy_edges_bac.csv", row.names = FALSE)

# Infected edges
infected_edges_bac <- spring_net_bac$edgelist2
colnames(infected_edges_bac) <- c("source", "target", "weight")

write.csv(infected_edges_bac, "infected_edges_bac.csv", row.names = FALSE)

nodes_healthy_bac <- unique(c(
  healthy_edges_bac$source,
  healthy_edges_bac$target
))

healthy_nodes_bac <- data.frame(
  name = nodes_healthy_bac
)

cent_h_bac <- spring_netprop_bac$centralities

healthy_nodes_bac$degree  <- cent_h_bac$degree1[healthy_nodes_bac$name]
healthy_nodes_bac$between <- cent_h_bac$between1[healthy_nodes_bac$name]
healthy_nodes_bac$close   <- cent_h_bac$close1[healthy_nodes_bac$name]
healthy_nodes_bac$eigenv  <- cent_h_bac$eigenv1[healthy_nodes_bac$name]

healthy_nodes_bac$hub <- healthy_nodes_bac$name %in% spring_netprop_bac$hubs[[1]]

healthy_nodes_bac$name <- gsub('"', '', healthy_nodes_bac$name)
healthy_edges_bac$source <- gsub('"', '', healthy_edges_bac$source)
healthy_edges_bac$target <- gsub('"', '', healthy_edges_bac$target)

write.csv(healthy_nodes_bac, "healthy_nodes_bac.csv",
          row.names = FALSE, quote = FALSE)

write.csv(healthy_edges_bac, "healthy_edges_bac.csv",
          row.names = FALSE, quote = FALSE)



nodes_infected_bac <- unique(c(
  infected_edges_bac$source,
  infected_edges_bac$target
))

infected_nodes_bac <- data.frame(
  name = nodes_infected_bac
)

cent_h_bac <- spring_netprop_bac$centralities

infected_nodes_bac$degree  <- cent_h_bac$degree2[infected_nodes_bac$name]
infected_nodes_bac$between <- cent_h_bac$between2[infected_nodes_bac$name]
infected_nodes_bac$close   <- cent_h_bac$close2[infected_nodes_bac$name]
infected_nodes_bac$eigenv  <- cent_h_bac$eigenv2[infected_nodes_bac$name]

infected_nodes_bac$hub <- infected_nodes_bac$name %in% spring_netprop_bac$hubs[[1]]

infected_nodes_bac$name <- gsub('"', '', infected_nodes_bac$name)
infected_edges_bac$source <- gsub('"', '', infected_edges_bac$source)
infected_edges_bac$target <- gsub('"', '', infected_edges_bac$target)

write.csv(infected_nodes_bac, "infected_nodes_bac.csv",
          row.names = FALSE, quote = FALSE)

write.csv(infected_edges_bac, "infected_edges_bac.csv",
          row.names = FALSE, quote = FALSE)

biomarker_df <- readRDS("data_ancombc2.RDS") %>%
  select(taxon, lfc, qval, diff) %>%
  rename(name = taxon)


sooty_taxa <- c(
  "Achaetobotrys","Antennulariella","Capnofrasera","Scolecoxyphium","Capnodium",
  "Strelitziana","Camptophora","Capnophaeum","Leptoxyphium","Phragmocapnia",
  "Scorias","Actinocymbe","Ceramothyrium","Chaetothyriomyces","Chaetothyrium",
  "Euceramia","Microcallis","Phaeosaccardinula","Treubiomyces","Yatesula",
  "Coccodinium","Dennisiella","Limacinula","Antennatula","Capnokyma",
  "Euantennaria","Hormisciomyces","Rasutoria","Strigopodia","Trichothallus",
  "Trichopeltheca","Capnocybe","Capnophialophora","Capnosporium","Hormiokrypsis",
  "Hyphosoma","Metacapnodium","Trichomerium","Aureobasidium","Cladosporium",
  "Epicoccum","Capnocheirides","Mycosphaerella","Phloespora","Ramichloridium",
  "Pseudoveronea","Capronia","Fusicladium","Coleroa","Venturia","Gibbera",
  "Protoventuria","Rhinocladiella","Vonarxia","Cladophialophora","Didymella",
  "Phaeococcuslike","Peltaster","Neopeltaster","Schizothyrium","Zygophiala",
  "Dissoconium","Uwebrauniacommune","Uwebraunia","Zasmidium","Devresia",
  "Microcyclospora","Houjia","Passaloralike","Phaeothecoidiella",
  "Translucidithyrium","Sporidesmajora","Chaetothyrina","Stomiopeltis",
  "Microcyclosporella","Colletogloeumlike","Geastrumia","Scolecobasidium",
  "Pleosporales","Chaetothyriales","Cyphellophora","Exophiala","Leptodontidium",
  "Neophaeococcomyces","Wallemia"
)
view(infected_nodes_fun)

plot_df <- infected_nodes_fun %>%
  left_join(biomarker_df, by = "name") %>%
  mutate(
    sooty = name %in% sooty_taxa,
    diff = as.logical(diff),
    sig  = diff == TRUE,
    lfc_plot = ifelse(sig, abs(lfc), 0.01)
  ) %>%
  filter(!is.na(lfc), !is.na(degree), !is.na(between))

# --- correlations (signed lfc)
cor_deg <- cor.test(plot_df$degree, plot_df$lfc, method = "pearson")
cor_bet <- cor.test(plot_df$between, plot_df$lfc, method = "pearson")
cor_betdeg <- cor.test(plot_df$between, plot_df$deg, method = "pearson")

corr <-ggplot(plot_df, aes(x = degree, y = between)) +
  geom_point(aes(size = lfc_plot, color = lfc), alpha = 0.85) +
  scale_color_gradient2(
    low = "blue",
    mid = "gray90",
    high = "red",
    midpoint = 0,
    name = "logFC"
  ) +
  scale_size_continuous(range = c(2, 8), guide = "none") +   # <-- remove size legend
  
  geom_text_repel(
    data = subset(plot_df, degree >= 2 & between > 0),
    aes(label = name),
    fontface = ifelse(subset(plot_df, degree >= 2 & between > 0)$sooty, "bold", "plain"),
    max.overlaps = Inf,
    force = 2,
    force_pull = 0.5,
    box.padding = 0.5,
    point.padding = 0.4,
    segment.size = 0.3,
    size = 3.5
  ) +
  
  labs(
    x = "Degree",
    y = "Betweenness",
    color = "logFC"
  ) +
  
  annotate(
    "text", x = Inf, y = Inf,
    label = paste0(
      "\nPearson (degree vs logFC): rho = ", round(cor_deg$estimate, 2),
      " | p = ", signif(cor_deg$p.value, 2),
      "\nPearson (between vs logFC): rho = ", round(cor_bet$estimate, 2),
      " | p = ", signif(cor_bet$p.value, 2),
      "\nPearson (between vs degree): rho = ", round(cor_betdeg$estimate, 2),
      " | p = ", signif(cor_betdeg$p.value, 2)
    ),
    hjust = 1.1, vjust = 1.1, size = 5
  ) +
  
  theme_nature +
  theme(legend.position = "bottom")

ggsave(
  filename = "figure.svg",
  plot = corr,          # your ggplot object
  device = "svg",
  width = 250,
  height = 250,
  units = "mm"
)
# DISSIMILARITY NETWORK ---------------------------------------------------
net_diss <- netConstruct(TSE_funlist$Healthy,
                         measure = "aitchison",
                         zeroMethod = "multRepl",
                         sparsMethod = "knn",
                         kNeighbor = 3,
                         verbose = 3)

net_diss2 <- netConstruct(TSE_funlist$Infected,
                         measure = "aitchison",
                         zeroMethod = "multRepl",
                         sparsMethod = "knn",
                         kNeighbor = 3,
                         verbose = 3)

props_diss <- netAnalyze(net_diss,
                         clustMethod = "hierarchical",
                         clustPar = list(method = "average", k = 3),
                         hubPar = "eigenvector")

props_diss2 <- netAnalyze(net_diss2,
                         clustMethod = "hierarchical",
                         clustPar = list(method = "average", k = 3),
                         hubPar = "eigenvector")

plot(props_diss, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     hubTransp = 40,
     edgeTranspLow = 60,
     charToRm = "00000",
     shortenLabels = "simple",
     labelLength = 6,
     mar = c(1, 3, 3, 5))

# get green color with 50% transparency
green2 <- colToTransp("#009900", 40)

legend(0.4, 1.1,
       cex = 2.2,
       legend = c("high similarity (low Aitchison distance)",
                  "low similarity (high Aitchison distance)"), 
       lty = 1, 
       lwd = c(3, 1),
       col = c("darkgreen", green2),
       bty = "n")


plot(props_diss2, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     hubTransp = 40,
     edgeTranspLow = 60,
     charToRm = "00000",
     shortenLabels = "simple",
     labelLength = 6,
     mar = c(1, 3, 3, 5))

# get green color with 50% transparency
green2 <- colToTransp("#009900", 40)

legend(0.4, 1.1,
       cex = 2.2,
       legend = c("high similarity (low Aitchison distance)",
                  "low similarity (high Aitchison distance)"), 
       lty = 1, 
       lwd = c(3, 1),
       col = c("darkgreen", green2),
       bty = "n")


# DIFFERENTIAL NETWORK ----------------------------------------------------
net_season_pears <- netConstruct(data = TSE_funlist$Healthy,
                                 data2 = TSE_funlist$Infected,
                                 filtTax = "highestVar",
                                 filtTaxPar = list(highestVar = 50),
                                 measure = "pearson", 
                                 normMethod = "clr",
                                 sparsMethod = "none", 
                                 thresh = 0.2,
                                 verbose = 3)

diff_season <- diffnet(net_season_pears,
                       diffMethod = "fisherTest", 
                       adjust = "lfdr")
plot(diff_season, 
     cexNodes = 0.8, 
     cexLegend = 3,
     cexTitle = 4,
     mar = c(2,2,8,5),
     legendGroupnames = c("group 'no'", "group 'yes'"),
     legendPos = c(0.7,1.6))

props_season_pears <- netAnalyze(net_season_pears, 
                                 clustMethod = "cluster_fast_greedy",
                                 weightDeg = TRUE,
                                 normDeg = FALSE,
                                 gcmHeat = FALSE)

diffmat_sums <- rowSums(diff_season$diffAdjustMat)
diff_asso_names <- names(diffmat_sums[diffmat_sums > 0])

plot(props_season_pears, 
     nodeFilter = "names",
     nodeFilterPar = diff_asso_names,
     nodeColor = "gray",
     highlightHubs = FALSE,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "clr",
     edgeTranspHigh = 20,
     labelScale = FALSE,
     cexNodes = 1.5, 
     cexLabels = 3,
     cexTitle = 3.8,
     groupNames = c("No seasonal allergies", "Seasonal allergies"),
     hubBorderCol  = "gray40")


# BACTERIAL NETWORKS ------------------------------------------------------

colData(TSE_bac)$class_label <- ifelse(colData(TSE_bac)$class == 'healthy', 'Healthy',
                                       ifelse(colData(TSE_bac)$class != 'healthy', 'Infected', NA))
# Agglomerate to Genus level
TSE_bac <- agglomerateByRank(TSE_bac, rank = "Genus")

table(TSE_bac$class_label)

TSE_baclist <- splitOn(TSE_bac, group = "class_label", use.names = TRUE, by = "cols")

TSE_baclist$Healthy <- subsetByPrevalent(
  TSE_baclist$Healthy,
  prevalence = 0.15,
  detection = 0,
  assay.type = "counts"
)

TSE_baclist$Infected <- subsetByPrevalent(
  TSE_baclist$Infected,
  prevalence = 0.15,
  detection = 0,
  assay.type = "counts"
)


# Extract Genus column
genus <- rowData(TSE_baclist$Infected)$Genus

# Replace NA or empty strings with "Unknown"
genus[is.na(genus) | genus == ""] <- paste0("Unknown_", seq_len(sum(is.na(genus) | genus == "")))

# Make duplicates unique
genus <- make.unique(genus, sep = "_")

# Assign back
rowData(TSE_baclist$Infected)$Genus <- genus

# Check uniqueness
anyDuplicated(rowData(TSE_baclist$Infected)$Genus)
# Should now return 0

genus <- rowData(TSE_baclist$Healthy)$Genus
genus[is.na(genus) | genus == ""] <- paste0("Unknown_", seq_len(sum(is.na(genus) | genus == "")))
genus <- make.unique(genus, sep = "_")
rowData(TSE_baclist$Healthy)$Genus <- genus


spring_net_bac <- netConstruct(
  data = TSE_baclist$Healthy,
  data2 = TSE_baclist$Infected,
  taxRank = "Genus",
  measure = "spring",
  measurePar = list(
    nlambda = 20,
    rep.num = 100,
    thresh = 0.05,
    Rmethod = "approx"),
  sparsMethod = "none",
  dissFunc = "signed",
  verbose = 3,
  seed = 13075)

spring_netprop_bac <- netAnalyze(
  spring_net_bac,
  centrLCC = FALSE,
  avDissIgnoreInf = TRUE,
  sPathNorm = FALSE,
  clustMethod = "cluster_fast_greedy",
  hubPar = c("degree", "eigenvector"),
  hubQuant = 0.9,
  lnormFit = TRUE,
  normDeg = FALSE,
  normBetw = FALSE,
  normClose = FALSE,
  normEigen = T)

summary(spring_netprop_bac, groupNames = c("Healthy", "Infected"))

# Adjust margins and plotting parameters
plot(spring_netprop_bac,
     sameLayout = TRUE, 
     repulsion = 0.98,
     borderCol = "gray40", 
     nodeSizeSpread = 5, 
     edgeTranspLow = 80, 
     edgeTranspHigh = 50,
     layoutGroup = "union",
     rmSingles = "inboth", 
     nodeSize = "mclr", 
     labelScale = T,
     labels = T,
     mar = c(1,1,3,1),
     cexNodes = 2.5, 
     cexLabels = 2.5,
     cexHubLabels = 1,
     cexTitle = 2,
     nodeFilter = "clustMin", 
     nodeFilterPar = 1, 
     nodeTransp = 30, 
     hubTransp = 30,
     groupNames = c("Healthy", "Infected"),
     hubBorderCol  = "gray40", 
     curve = 0.2,
     curveAll = TRUE)

legend("bottom", title = "estimated association:", legend = c("+","-"), 
       col = c("#009900","red"), inset = 0.02, cex = 4, lty = 1, lwd = 4, 
       bty = "n", horiz = TRUE)

# DISSIMILARITY NETWORK ---------------------------------------------------
net_diss_bac <- netConstruct(TSE_baclist$Healthy,
                         measure = "aitchison",
                         zeroMethod = "multRepl",
                         sparsMethod = "knn",
                         kNeighbor = 3,
                         verbose = 3)

net_diss2_bac <- netConstruct(TSE_baclist$Infected,
                          measure = "aitchison",
                          zeroMethod = "multRepl",
                          sparsMethod = "knn",
                          kNeighbor = 3,
                          verbose = 3)

props_diss_bac <- netAnalyze(net_diss_bac,
                         clustMethod = "hierarchical",
                         clustPar = list(method = "average", k = 3),
                         hubPar = "eigenvector")

props_diss2_bac <- netAnalyze(net_diss2_bac,
                          clustMethod = "hierarchical",
                          clustPar = list(method = "average", k = 3),
                          hubPar = "eigenvector")

plot(props_diss_bac, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     hubTransp = 40,
     edgeTranspLow = 60,
     charToRm = "00000",
     shortenLabels = "simple",
     labelLength = 6,
     mar = c(1, 3, 3, 5))

# get green color with 50% transparency
green2 <- colToTransp("#009900", 40)

legend(0.4, 1.1,
       cex = 2.2,
       legend = c("high similarity (low Aitchison distance)",
                  "low similarity (high Aitchison distance)"), 
       lty = 1, 
       lwd = c(3, 1),
       col = c("darkgreen", green2),
       bty = "n")


plot(props_diss2_bac, 
     nodeColor = "cluster", 
     nodeSize = "eigenvector",
     hubTransp = 40,
     edgeTranspLow = 60,
     charToRm = "00000",
     shortenLabels = "simple",
     labelLength = 6,
     mar = c(1, 3, 3, 5))

# get green color with 50% transparency
green2 <- colToTransp("#009900", 40)

legend(0.4, 1.1,
       cex = 2.2,
       legend = c("high similarity (low Aitchison distance)",
                  "low similarity (high Aitchison distance)"), 
       lty = 1, 
       lwd = c(3, 1),
       col = c("darkgreen", green2),
       bty = "n")


# DIFFERENTIAL NETWORK ----------------------------------------------------
net_season_pears_bac <- netConstruct(data = TSE_baclist$Healthy,
                                 data2 = TSE_baclist$Infected,
                                 filtTax = "highestVar",
                                 filtTaxPar = list(highestVar = 50),
                                 measure = "pearson", 
                                 normMethod = "clr",
                                 sparsMethod = "none", 
                                 thresh = 0.2,
                                 verbose = 3)

diff_season_bac <- diffnet(net_season_pears_bac,
                       diffMethod = "fisherTest", 
                       adjust = "lfdr")
plot(diff_season_bac, 
     cexNodes = 0.8, 
     cexLegend = 3,
     cexTitle = 4,
     mar = c(2,2,8,5),
     legendGroupnames = c("group 'no'", "group 'yes'"),
     legendPos = c(0.7,1.6))

props_season_pears_bac <- netAnalyze(net_season_pears_bac, 
                                 clustMethod = "cluster_fast_greedy",
                                 weightDeg = TRUE,
                                 normDeg = FALSE,
                                 gcmHeat = FALSE)

diffmat_sums_bac <- rowSums(diff_season_bac$diffAdjustMat)
diff_asso_names_bac <- names(diffmat_sums_bac[diffmat_sums_bac > 0])

plot(props_season_pears_bac, 
     nodeFilter = "names",
     nodeFilterPar = diff_asso_names_bac,
     nodeColor = "gray",
     highlightHubs = FALSE,
     sameLayout = TRUE, 
     layoutGroup = "union",
     rmSingles = FALSE, 
     nodeSize = "clr",
     edgeTranspHigh = 20,
     labelScale = FALSE,
     cexNodes = 1.5, 
     cexLabels = 3,
     cexTitle = 3.8,
     groupNames = c("No seasonal allergies", "Seasonal allergies"),
     hubBorderCol  = "gray40")

# CROSS DOMAIN NETWORKS  --------------------------------------------------
# ============================================
# Cross-domain networks: Bacteria + Fungi
# Full pipeline: preprocessing -> SPIEC-EASI -> NetCoMi -> plotting
# ============================================
# ============================================
# Cross-domain networks: Bacteria + Fungi
# Fully runnable pipeline with union-of-taxa fix
# ============================================

# ----------------------------
# 1. Add class labels
# ----------------------------
sample_data(physeq_fun_filtered)$class_label <- ifelse(
  sample_data(physeq_fun_filtered)$class == "healthy", "Healthy", "Infected"
)
sample_data(physeq_bac_filtered)$class_label <- ifelse(
  sample_data(physeq_bac_filtered)$class == "healthy", "Healthy", "Infected"
)

# ----------------------------
# 2. Agglomerate at Genus level
# ----------------------------
physeq_fun_genus <- tax_glom(physeq_fun_filtered, taxrank = "Genus")
physeq_bac_genus <- tax_glom(physeq_bac_filtered, taxrank = "Genus")

physeq_fun_genus <- prune_taxa(taxa_sums(physeq_fun_genus) > 0, physeq_fun_genus)
physeq_bac_genus <- prune_taxa(taxa_sums(physeq_bac_genus) > 0, physeq_bac_genus)

# ----------------------------
# 3. Subset by Healthy / Infected
# ----------------------------
physeq_fun_healthy   <- subset_samples(physeq_fun_genus, class_label == "Healthy")
physeq_fun_infected  <- subset_samples(physeq_fun_genus, class_label == "Infected")
physeq_bac_healthy   <- subset_samples(physeq_bac_genus, class_label == "Healthy")
physeq_bac_infected  <- subset_samples(physeq_bac_genus, class_label == "Infected")

# ----------------------------
# 4. Prevalence filtering
# ----------------------------
filter_prevalence <- function(physeq, threshold = 0.1) {
  prev <- apply(otu_table(physeq) > 0, 1, sum)
  prev_frac <- prev / nsamples(physeq)
  prune_taxa(prev_frac >= threshold, physeq)
}

physeq_fun_healthy  <- filter_prevalence(physeq_fun_healthy)
physeq_fun_infected <- filter_prevalence(physeq_fun_infected)
physeq_bac_healthy  <- filter_prevalence(physeq_bac_healthy)
physeq_bac_infected <- filter_prevalence(physeq_bac_infected)

# ----------------------------
# 5. Convert phyloseq to count matrices
# ----------------------------
get_counts <- function(ps) {
  mat <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) mat <- t(mat)
  mat
}

bac_counts_healthy <- get_counts(physeq_bac_healthy)
fun_counts_healthy <- get_counts(physeq_fun_healthy)
bac_counts_infected <- get_counts(physeq_bac_infected)
fun_counts_infected <- get_counts(physeq_fun_infected)

# Healthy
intersect(
  sample_names(physeq_fun_healthy),
  sample_names(physeq_bac_healthy)
)

# Infected
intersect(
  sample_names(physeq_fun_infected),
  sample_names(physeq_bac_infected)
)

physeq_healthy_merged <- merge_phyloseq(
  physeq_fun_healthy,
  physeq_bac_healthy
)

physeq_infected_merged <- merge_phyloseq(
  physeq_fun_infected,
  physeq_bac_infected
)

rownames(physeq_healthy_merged) <- NULL

physeq_healthy_merged$Sample_id <- paste0("ASV", seq_len(nrow(physeq_healthy_merged)))

tax_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

physeq_healthy_merged[tax_cols] <- lapply(
  physeq_healthy_merged[tax_cols],
  function(x) ifelse(is.na(x), "NA", x)
)

physeq_healthy_merged$taxonomy <- apply(
  physeq_healthy_merged[, tax_cols],
  1,
  paste,
  collapse = ";"
)

count_cols <- setdiff(
  colnames(physeq_healthy_merged),
  c("Sample_id", tax_cols, "taxonomy")
)

physeq_healthy_merged_file <- physeq_healthy_merged[, c(
  "Sample_id",
  count_cols,
  "taxonomy"
)]

# ----------------------------
# 6. Run SPIEC-EASI
# ----------------------------
set.seed(123456)

spiec_healthy <- multi.spiec.easi(
  list(bac_counts_healthy, fun_counts_healthy),
  method = "mb",
  nlambda = 10,
  lambda.min.ratio = 1e-2,
  pulsar.params = list(thresh = 0.05, rep.num = 10)
)

spiec_infected <- multi.spiec.easi(
  list(bac_counts_infected, fun_counts_infected),
  method = "mb",
  nlambda = 10,
  lambda.min.ratio = 1e-2,
  pulsar.params = list(thresh = 0.05, rep.num = 10)
)

# ----------------------------
# 7. Extract association matrices
# ----------------------------
assoMat_healthy <- as.matrix(symBeta(getOptBeta(spiec_healthy), mode = "ave"))
assoMat_infected <- as.matrix(symBeta(getOptBeta(spiec_infected), mode = "ave"))

# ----------------------------
# 8. Align matrices by union of taxa
# ----------------------------
all_taxa <- union(rownames(assoMat_healthy), rownames(assoMat_infected))

fill_matrix <- function(mat, all_taxa) {
  mat_full <- matrix(0, nrow = length(all_taxa), ncol = length(all_taxa))
  rownames(mat_full) <- colnames(mat_full) <- all_taxa
  common <- intersect(rownames(mat), all_taxa)
  mat_full[common, common] <- mat[common, common]
  diag(mat_full) <- 1
  mat_full
}

assoMat1_full <- fill_matrix(assoMat_healthy, all_taxa)
assoMat2_full <- fill_matrix(assoMat_infected, all_taxa)

# ----------------------------
# 9. Construct NetCoMi network
# ----------------------------
net_bac_fun <- netConstruct(
  data = assoMat1_full,
  data2 = assoMat2_full,
  dataType = "condDependence",
  sparsMethod = "none"
)

netprops_bac_fun <- netAnalyze(net_bac_fun, hubPar = "eigenvector")

# ----------------------------
# 10. Assign node colors
# ----------------------------
net_nodes <- netprops_bac_fun$nodes$names

nodeCols <- rep(NA, length(net_nodes))
names(nodeCols) <- net_nodes

bac_nodes <- intersect(net_nodes, union(colnames(bac_counts_healthy), colnames(bac_counts_infected)))
fun_nodes <- intersect(net_nodes, union(colnames(fun_counts_healthy), colnames(fun_counts_infected)))

nodeCols[bac_nodes] <- "lightblue"
nodeCols[fun_nodes] <- "orange"

stopifnot(!any(is.na(nodeCols)))
table(nodeCols)

# ----------------------------
# 11. Plot network
# ----------------------------
plot(netprops_bac_fun,
     sameLayout = TRUE,
     layoutGroup = "union",
     nodeColor = "colorVec",
     colorVec = nodeCols,
     nodeSize = "eigen",
     nodeSizeSpread = 2,
     labelScale = FALSE,
     cexNodes = 2,
     cexLabels = 2,
     cexHubLabels = 2.5,
     cexTitle = 3.8,
     groupNames = c("Healthy", "Infected"))

legend("topright",
       legend = c("Bacteria", "Fungi"),
       col = c("lightblue", "orange"),
       pch = 16,
       pt.cex = 2,
       cex = 1.5,
       bty = "n")

# ----------------------------
# 12. Network comparison
# ----------------------------
netcomp_bac_fun <- netCompare(netprops_bac_fun, permTest = FALSE)
summary(netcomp_bac_fun, groupNames = c("Healthy", "Infected"))

# Optional: extract hubs and global metrics
netcomp_bac_fun$hubs
netcomp_bac_fun$glo

# -----------------------------
# 1. Create igraph objects
# -----------------------------
g1 <- graph_from_data_frame(spring_net_fun$edgelist1[, c("v1","v2")], directed = FALSE)
g2 <- graph_from_data_frame(spring_net_fun$edgelist2[, c("v1","v2")], directed = FALSE)

# -----------------------------
# 2. Compute centralities
# -----------------------------
central1 <- data.frame(
  taxa = V(g1)$name,
  degree = igraph::degree(g1),
  betweenness = igraph::betweenness(g1),
  eigenvector = igraph::eigen_centrality(g1)$vector,
  group = "Healthy"
)

central2 <- data.frame(
  taxa = V(g2)$name,
  degree = igraph::degree(g2),
  betweenness = igraph::betweenness(g2),
  eigenvector = igraph::eigen_centrality(g2)$vector,
  group = "Infected"
)

central_all <- rbind(central1, central2)

# -----------------------------
# 3. Add mean abundance
# -----------------------------
abund1 <- colMeans(spring_net$normCounts1)
abund2 <- colMeans(spring_net$normCounts2)

central1$abundance <- abund1[central1$taxa]
central2$abundance <- abund2[central2$taxa]
central_all <- rbind(central1, central2)

# -----------------------------
# 4. Define hubs (90th percentile)
# -----------------------------
hub_thr <- 0.90

central_all$hub <- with(central_all,
                        degree >= ave(degree, group, FUN = function(x) quantile(x, hub_thr)) &
                          betweenness >= ave(betweenness, group, FUN = function(x) quantile(x, hub_thr)) &
                          eigenvector >= ave(eigenvector, group, FUN = function(x) quantile(x, hub_thr))
)

table(central_all$group, central_all$hub)

# -----------------------------
# 5. Hub identification plot
# -----------------------------
hub_plot <- ggplot(central_all, aes(x = degree, y = betweenness, size = eigenvector, color = hub)) +
  geom_point(alpha = 0.8) +
  facet_wrap(~ group) +
  scale_y_log10() +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  geom_text_repel(
    data = subset(central_all, hub == TRUE),
    aes(label = taxa),
    size = 3
  ) +
  geom_hline(data = central_all,
             aes(yintercept = ave(betweenness, group, FUN = function(x) quantile(x, hub_thr))),
             linetype = "dashed", color = "darkgray") +
  geom_vline(data = central_all,
             aes(xintercept = ave(degree, group, FUN = function(x) quantile(x, hub_thr))),
             linetype = "dashed", color = "darkgray") +
  theme_minimal() +
  labs(
    title = "Hub taxa identification (Healthy vs Infected)",
    x = "Degree",
    y = "Betweenness (log scale)",
    size = "Eigenvector centrality",
    color = "Hub taxa"
  ) +
  theme_nature

hub_plot

# -----------------------------
# 6. Define sooty taxa
# -----------------------------
sooty_taxa <- c(
  "Achaetobotrys","Antennulariella","Capnofrasera","Scolecoxyphium","Capnodium",
  "Strelitziana","Camptophora","Capnophaeum","Leptoxyphium","Phragmocapnia",
  "Scorias","Actinocymbe","Ceramothyrium","Chaetothyriomyces","Chaetothyrium",
  "Euceramia","Microcallis","Phaeosaccardinula","Treubiomyces","Yatesula",
  "Coccodinium","Dennisiella","Limacinula","Antennatula","Capnokyma",
  "Euantennaria","Hormisciomyces","Rasutoria","Strigopodia","Trichothallus",
  "Trichopeltheca","Capnocybe","Capnophialophora","Capnosporium","Hormiokrypsis",
  "Hyphosoma","Metacapnodium","Trichomerium","Aureobasidium","Cladosporium",
  "Epicoccum","Capnocheirides","Mycosphaerella","Phloespora","Ramichloridium",
  "Pseudoveronea","Capronia","Fusicladium","Coleroa","Venturia","Gibbera",
  "Protoventuria","Rhinocladiella","Vonarxia","Cladophialophora","Didymella",
  "Phaeococcuslike","Peltaster","Neopeltaster","Schizothyrium","Zygophiala",
  "Dissoconium","Uwebrauniacommune","Uwebraunia","Zasmidium","Devresia",
  "Microcyclospora","Houjia","Passaloralike","Phaeothecoidiella",
  "Translucidithyrium","Sporidesmajora","Chaetothyrina","Stomiopeltis",
  "Microcyclosporella","Colletogloeumlike","Geastrumia","Scolecobasidium",
  "Pleosporales","Chaetothyriales","Cyphellophora","Exophiala","Leptodontidium",
  "Neophaeococcomyces","Wallemia"
)

central2$sooty <- central2$taxa %in% sooty_taxa
central2$sooty
# -----------------------------
# 7. Load biomarker scores
# -----------------------------
# -----------------------------
# 7. Load biomarker scores
# -----------------------------
biomarker_df <- readRDS("lefse_result.rds")

# Check columns
colnames(biomarker_df)

# Join LEfSe results to infected network centralities
central_bio <- central2 %>%
  dplyr::left_join(biomarker_df, by = c("taxa" = "Taxon"))

# -----------------------------
# 8. Recompute hubs in infected network
# -----------------------------
hub_thr <- 0.90

central_bio <- central_bio %>%
  mutate(
    hub = degree >= quantile(degree, hub_thr, na.rm = TRUE) &
      betweenness >= quantile(betweenness, hub_thr, na.rm = TRUE) &
      eigenvector >= quantile(eigenvector, hub_thr, na.rm = TRUE),
    biomarker = !is.na(LDA) & LDA > 0
  )

# -----------------------------
# 9. Biomarker vs Degree plot
# -----------------------------
p_biomarker_degree <- ggplot(
  central_bio,
  aes(x = SignedLDA, y = degree, size = betweenness, color = biomarker)
) +
  geom_point(alpha = 0.8) +
  geom_text_repel(
    data = subset(central_bio, hub == TRUE | biomarker == TRUE),
    aes(label = taxa),
    size = 4.5,          # larger labels
    max.overlaps = 45
  ) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  theme_minimal() +
  labs(
    x = "Signed LEfSe score",
    y = "Degree",
    size = "Betweenness",
    color = "Biomarker"
  ) +
  theme_nature

p_biomarker_degree
