# ---- Core data handling & manipulation ----
library(tidyverse)
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
library(usethis)
library(nlMS)
library(iCAMP)
library(minpack.lm)
library(Hmisc)
library(microbiome)

# Update BiocManager packages including ANCOMBC
BiocManager::install("microbiome", update = FALSE, ask = FALSE)

library(microbiome)
##### Importing files ########

biom = import_biom("Inputfiles/bacteria/season-16s.biom")

metadata = import_qiime_sample_data("Inputfiles/bacteria/metadata.txt")

#tree = read_tree("Input_files/rooted_tree.nwk")

#rep_fasta = readDNAStringSet("Input_files/fun-seq.fasta", format = "fasta")

season_biom = merge_phyloseq(biom, metadata)

colnames(tax_table(season_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meco_dataset <- phyloseq2meco(season_biom)

soil <- read.csv("Inputfiles/bacteria/soil_data.csv")

rownames(soil) <- soil[, 1]

soil = soil[ ,-1]

# add_data is used to add the environmental data
soil_season <- trans_env$new(dataset = meco_dataset, add_data = soil)

meco_dataset$cal_abund()
meco_dataset$cal_alphadiv()
meco_dataset$cal_betadiv()

meco_dataset$sample_table <- meco_dataset$sample_table %>%
  mutate(Season = factor(Season, levels = c("Fall", "Spring")),
         tillage = factor(tillage, levels = c("Min Till", "Conv Till")),
         Residue = factor(cover_crop, levels = c("CC-mix", "no CC")),
         Crop_rotation = factor(Crop_rotation, levels = c("Cotton-Soybean", "Corn-Soybean", "Soybean-Corn", "Soybean-Soybean")))

ps_ancom <- phyloseq(
  otu_table(as.matrix(meco_dataset$otu_table), taxa_are_rows = TRUE),
  tax_table(as.matrix(meco_dataset$tax_table)),
  sample_data(meco_dataset$sample_table)
)

ps_ancom
sample_variables(ps_ancom)
rank_names(ps_ancom)

otu_table(ps_ancom)[1:5, 1:5]
sample_sums(ps_ancom) %>% summary()

set.seed(123)

meta <- data.frame(sample_data(ps_ancom)) %>%
  mutate(Season = factor(Season, levels = c("Fall", "Spring")),
         tillage = factor(tillage, levels = c("Min Till", "Conv Till")),
         cover_crop = factor(cover_crop, levels = c("CC-mix", "no CC")),
         Crop_rotation = factor(Crop_rotation, levels = c("Cotton-Soybean", "Corn-Soybean", "Soybean-Corn", "Soybean-Soybean")))

sample_data(ps_ancom) <- sample_data(meta)

tax_ranks <- c("Phylum", "Class", "Order", "Family", "Genus")

factor_table <- tribble(
  ~Factor,    ~Variable,
  "Season",   "Season",
  "Tillage",  "tillage",
  "Residue",  "cover_crop"
)

rank_prefix <- c(Phylum = "p", Class = "c", Order = "o", Family = "f", Genus = "g")

clean_taxon <- function(x, rank) {
  x %>% str_replace(paste0(".*", rank_prefix[[rank]], "__"), "") %>%
    str_replace_all("_", " ") %>%
    str_trim()
}

run_ancom_binary <- function(ps, variable, factor_name, rank) {
  
  set.seed(123)
  
  fit <- ancombc2(
    data = ps,
    tax_level = rank,
    fix_formula = variable,
    rand_formula = NULL,
    p_adj_method = "BH",
    pseudo_sens = TRUE,
    prv_cut = 0.10,
    lib_cut = 0,
    s0_perc = 0.05,
    group = variable,
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 4,
    verbose = FALSE,
    global = FALSE,
    pairwise = FALSE,
    dunnet = FALSE,
    trend = FALSE
  )
  
  res <- fit$res
  
  lfc_col <- names(res) %>%
    stringr::str_subset("^lfc_") %>%
    stringr::str_subset("Intercept", negate = TRUE) %>%
    .[1]
  
  coef_name <- stringr::str_remove(lfc_col, "^lfc_")
  
  se_col <- paste0("se_", coef_name)
  p_col <- paste0("p_", coef_name)
  q_col <- paste0("q_", coef_name)
  diff_col <- paste0("diff_", coef_name)
  ss_col <- paste0("passed_ss_", coef_name)
  robust_col <- paste0("diff_robust_", coef_name)
  
  lev <- levels(data.frame(sample_data(ps))[[variable]])
  
  robust <- if (robust_col %in% names(res)) {
    res[[robust_col]]
  } else {
    dplyr::coalesce(res[[diff_col]], FALSE) &
      dplyr::coalesce(res[[ss_col]], FALSE)
  }
  
  res %>%
    transmute(
      Factor = factor_name,
      Rank = rank,
      taxon,
      Taxon = clean_taxon(taxon, rank),
      Contrast = coef_name,
      Reference = lev[1],
      Comparison = lev[2],
      LFC = .data[[lfc_col]],
      SE = .data[[se_col]],
      p = .data[[p_col]],
      q = .data[[q_col]],
      diff = .data[[diff_col]],
      passed_ss = .data[[ss_col]],
      Significant = robust,
      Enriched = case_when(
        robust & .data[[lfc_col]] > 0 ~ lev[2],
        robust & .data[[lfc_col]] < 0 ~ lev[1],
        TRUE ~ NA_character_
      )
    )
}

binary_results <- crossing(factor_table, Rank = tax_ranks) %>%
  mutate(Result = pmap(list(Variable, Factor, Rank),
                       ~run_ancom_binary(ps_ancom, ..1, ..2, ..3))) %>%
  select(Result) %>%
  unnest(Result)

binary_summary <- binary_results %>%
  dplyr::group_by(Factor, Rank) %>%
  dplyr::summarise(
    Tested = dplyr::n(),
    Significant = sum(Significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(Rank = factor(Rank, levels = tax_ranks))

binary_summary

binary_summary %>%
  select(Factor, Rank, Significant) %>%
  pivot_wider(names_from = Rank, values_from = Significant, values_fill = 0)

season_results <- binary_results %>%
  filter(Factor == "Season") %>%
  mutate(Rank = factor(Rank, levels = c("Phylum", "Class", "Order", "Family", "Genus")))

season_sig <- season_results %>%
  filter(Significant) %>%
  mutate(Direction = if_else(LFC > 0, "Spring", "Fall"))

head(season_sig)

p_season_counts <- season_sig %>%
  dplyr::count(Rank) %>%
  ggplot(aes(x = Rank, y = n)) +
  geom_col(width = 0.7, fill = "steelblue") +
  geom_text(aes(label = n), vjust = -0.3, size = 5) +
  labs(
    x = "Taxonomic rank", 
    y = "Number of significant taxa",
    title = "Season-associated biomarkers across taxonomic levels",
    subtitle = "ANCOM-BC2 significant taxa"
  ) +
  theme_bw(base_size = 13)

p_season_counts

ggsave("Season_significant_taxa_counts.pdf", p_season_counts, width = 10, height = 16, dpi = 1000)

top_n_rank <- 15

season_top_by_rank <- season_sig %>%
  group_by(Rank) %>%
  slice_max(order_by = abs(LFC), n = top_n_rank, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Taxon = fct_reorder(Taxon, LFC))

p_season_top_rank <- season_top_by_rank %>%
  ggplot(aes(LFC, Taxon, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, linetype = 2) +
  facet_wrap(~Rank, scales = "free_y", ncol = 1) +
  labs(x = "ANCOM-BC2 log fold change", y = NULL,
       title = "Top season-associated biomarkers by taxonomic rank",
       subtitle = "Positive = Spring enriched; Negative = Fall enriched",
       fill = NULL) +
  theme_bw(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        axis.text.y = element_text(size = 9))

p_season_top_rank

ggsave("Season_top_biomarkers_by_rank.pdf", p_season_top_rank, width = 10, height = 16, dpi = 1000)

season_genus_clean <- binary_results %>%
  filter(Factor == "Season", Rank == "Genus", Significant) %>%
  mutate(Taxon = na_if(Taxon, ""),
         Taxon = if_else(is.na(Taxon), "Unclassified", Taxon)) %>%
  filter(!Taxon %in% c("Unclassified", "uncultured", "metagenome", "Ambiguous_taxa")) %>%
  mutate(Direction = if_else(LFC > 0, "Spring", "Fall"))

season_genus_top <- season_genus_clean %>%
  slice_max(abs(LFC), n = 20, with_ties = FALSE) %>%
  mutate(Taxon = forcats::fct_reorder(Taxon, LFC))

p_season_genus <- season_genus_top %>%
  ggplot(aes(LFC, Taxon, fill = Direction)) +
  geom_col(width = 0.7) +
  geom_errorbar(aes(xmin = LFC - SE, xmax = LFC + SE), width = 0.2) +
  geom_vline(xintercept = 0, linetype = 2) +
  labs(x = "ANCOM-BC2 log fold change", y = NULL,
       title = "Top season-associated genera",
       subtitle = "Positive = Spring enriched; Negative = Fall enriched",
       fill = NULL) +
  theme_bw(base_size = 13)

p_season_genus

ggsave("p_season_genus.pdf", p_season_genus, width = 10, height = 9, dpi = 1000)

