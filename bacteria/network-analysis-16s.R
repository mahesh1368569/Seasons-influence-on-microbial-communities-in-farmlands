# ---- Core data handling & manipulation ----
library(dplyr)
library(tidyr)
library(tidyverse)
library(plyr)
library(magrittr)
library(parallel)
library(reshape2)

# ---- Microbiome data processing ----
library(phyloseq)
library(biomformat)
library(file2meco)
library(microeco)
library(MicrobiomeStat)
library(meconetcomp)
library(WGCNA)
library(ggClusterNet)
library(ape)
library(picante)
library(Biostrings)

# ---- Differential abundance & compositional analysis ----
library(metagenomeSeq)
library(ALDEx2)
library(ANCOMBC)

# ---- Visualization ----
library(ggplot2)
library(ggpubr)
library(ggtree)
library(tidygraph)
library(paletteer)
library(colorspace)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(ggraph)

# ---- Statistical modeling ----
library(lme4)
library(lmerTest)
library(multcomp)
library(emmeans)
library(multcompView)
library(dplyr)
library(usethis)
library(nlMS)
library(iCAMP)
library(minpack.lm)
library(Hmisc)

##### Importing files ########

biom = import_biom("bacteria/Inputfiles/bacteria/season-16s.biom")

metadata = import_qiime_sample_data("bacteria/Inputfiles/bacteria/metadata.txt")

#tree = read_tree("Input_files/rooted_tree.nwk")

#rep_fasta = readDNAStringSet("Input_files/fun-seq.fasta", format = "fasta")

season_biom = merge_phyloseq(biom, metadata)

colnames(tax_table(season_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meco_dataset <- phyloseq2meco(season_biom)

soil <- read.csv("bacteria/Inputfiles/bacteria/soil_data.csv")

rownames(soil) <- soil[, 1]

soil = soil[ ,-1]

# add_data is used to add the environmental data
soil_season <- trans_env$new(dataset = meco_dataset, add_data = soil)

meco_dataset$cal_abund()
meco_dataset$cal_alphadiv()
meco_dataset$cal_betadiv()

meco_dataset$sample_table %>%
  as.data.frame() %>%
  dplyr::select(Year, farm, tillage, cover_crop, Crop_rotation, Season) %>%
  head()

factor_vars <- c(
  "Season",
  "Year",
  "tillage",
  "cover_crop",
  "Crop_rotation",
  "farm"
)

sample_sizes <- purrr::map_dfr(
  factor_vars,
  \(x) meco_dataset$sample_table %>%
    as.data.frame() %>%
    dplyr::mutate(Level = as.character(.data[[x]])) %>%
    dplyr::count(Factor = x, Level, name = "Samples")
)

sample_sizes

#============================================================
# FALL NETWORK
#============================================================

fall_data <- clone(meco_dataset)

fall_data$sample_table <- fall_data$sample_table %>%
  as.data.frame() %>%
  dplyr::filter(Season == "Fall")

fall_data$tidy_dataset()

fall_net <- trans_network$new(
  dataset = fall_data,
  taxa_level = "Genus",
  cor_method = "spearman",
  use_WGCNA_pearson_spearman = TRUE,
  filter_thres = 0.001
)

fall_net$cal_network(
  COR_p_thres = 0.01,
  COR_cut = 0.6,
  COR_p_adjust = "fdr",
  usename_rawtaxa_notOTU = TRUE
)

#============================================================
# SPRING NETWORK
#============================================================

spring_data <- clone(meco_dataset)

spring_data$sample_table <- spring_data$sample_table %>%
  as.data.frame() %>%
  dplyr::filter(Season == "Spring")

spring_data$tidy_dataset()

spring_net <- trans_network$new(
  dataset = spring_data,
  taxa_level = "Genus",
  cor_method = "spearman",
  use_WGCNA_pearson_spearman = TRUE,
  filter_thres = 0.001
)

spring_net$cal_network(
  COR_p_thres = 0.01,
  COR_cut = 0.6,
  COR_p_adjust = "fdr",
  usename_rawtaxa_notOTU = TRUE
)

season_network <- list(
  Fall = fall_net,
  Spring = spring_net
)

season_network <- meconetcomp::cal_module(
  season_network,
  undirected_method = "cluster_fast_greedy"
)

season_topology <- meconetcomp::cal_network_attr(
  season_network
)

season_topology

write.csv(
  season_topology,
  "Output/Network/Season_network_topology.csv",
  row.names = TRUE
)

season_network <- meconetcomp::get_node_table(
  season_network,
  node_roles = FALSE
)

season_network <- meconetcomp::get_edge_table(
  season_network
)


head(season_network$Fall$res_node_table)
head(season_network$Spring$res_node_table)

head(season_network$Fall$res_edge_table)
head(season_network$Spring$res_edge_table)

p_fall <- season_network$Fall$plot_network(
  method = "ggraph",
  node_color = "Phylum"
)

p_spring <- season_network$Spring$plot_network(
  method = "ggraph",
  node_color = "Phylum"
)

p_fall
p_spring
