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

# Generate genus level abundance

genus_ra <- meco_dataset$taxa_abund$Genus

colSums(genus_ra) %>% summary()

## Convert it to long format and attach metadata

metadata <- meco_dataset$sample_table %>%
  rownames_to_column("Sample")

genus_long <- genus_ra %>%
  as.data.frame() %>%
  rownames_to_column("Genus") %>%
  pivot_longer(-Genus, names_to = "Sample", values_to = "RA") %>%
  left_join(metadata, by = "Sample")

head(genus_long)

## Define the seasonal core

core_threshold <- 0.80
detection_threshold <- 0

core_stats <- genus_long %>%
  dplyr::group_by(Season, Genus) %>%
  dplyr::summarise(
    n_samples = dplyr::n(),
    n_detected = sum(RA > detection_threshold),
    prevalence = mean(RA > detection_threshold),
    mean_RA = mean(RA),
    median_RA = median(RA),
    .groups = "drop"
  ) %>%
  mutate(Core = prevalence >= core_threshold)

core_stats

core_stats %>%
  ggplot(aes(x = prevalence, fill = Season)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 20) +
  facet_wrap(~Season) +
  theme_minimal()

genus_long %>%
  dplyr::group_by(Sample, Season) %>%
  dplyr::summarise(n_genera = sum(RA > 0), .groups = "drop") %>%
  ggplot(aes(x = Season, y = n_genera, fill = Season)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  theme_minimal()


fall_only_stats

## Get the spring and fall cores

core_spring <- core_stats %>%
  filter(Season == "Spring", Core) %>%
  pull(Genus)

core_fall <- core_stats %>%
  filter(Season == "Fall", Core) %>%
  pull(Genus)

length(core_spring)
length(core_fall)

# Identify the three categories

shared_core <- intersect(core_spring, core_fall)
spring_core_only <- setdiff(core_spring, core_fall)
fall_core_only <- setdiff(core_fall, core_spring)

length(shared_core)
length(spring_core_only)
length(fall_core_only)

## Make a table showing core status

core_taxa <- union(core_spring, core_fall)

core_classification <- tibble(Genus = core_taxa) %>%
  mutate(
    Core_status = case_when(
      Genus %in% shared_core ~ "Shared core",
      Genus %in% spring_core_only ~ "Spring core only",
      Genus %in% fall_core_only ~ "Fall core only"
    )
  )

core_classification

# Add the actual seasonal prevalence values

core_summary <- core_stats %>%
  filter(Genus %in% core_taxa) %>%
  dplyr::select(Season, Genus, prevalence, mean_RA) %>%
  pivot_wider(
    names_from = Season,
    values_from = c(prevalence, mean_RA)
  ) %>%
  left_join(core_classification, by = "Genus")

core_summary

dir.create("bacteria/Output/Core_microbiome", recursive = TRUE, showWarnings = FALSE)

write_csv(
  core_summary,
  "bacteria/Output/Core_microbiome/Core_genera_by_season.csv"
)

sample_order <- metadata %>%
  arrange(Season, farm, replication) %>%
  pull(Sample)

genus_order <- core_stats %>%
  filter(Genus %in% core_taxa) %>%
  group_by(Genus) %>%
  dplyr::summarise(mean_RA = mean(mean_RA, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_RA) %>%
  dplyr::pull(Genus)

core_heatmap_dat <- genus_long %>%
  filter(Genus %in% core_taxa) %>%
  mutate(
    Sample = factor(Sample, levels = sample_order),
    Genus = factor(Genus, levels = genus_order),
    RA_percent = RA * 100
  )

core_heatmap <- ggplot(
  core_heatmap_dat,
  aes(x = Sample, y = Genus, fill = RA_percent)
) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(~Season, scales = "free_x", space = "free_x") +
  scale_fill_viridis_c(
    option = "C",
    trans = "sqrt",
    name = "Relative\nabundance (%)"
  ) +
  labs(x = "Soil sample", y = "Core genus") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
  )

core_heatmap

ggsave(
  "bacteria/Output/Core_microbiome/Core_genera_heatmap.pdf",
  core_heatmap,
  width = 11,
  height = 8
)

library(tidyverse)

#------------------------------------------------------------
# Metadata
#------------------------------------------------------------
metadata <- meco_dataset$sample_table %>%
  rownames_to_column("Sample")

#------------------------------------------------------------
# Taxonomic prefixes
#------------------------------------------------------------
rank_prefix <- c(
  Phylum = "p__",
  Class = "c__",
  Order = "o__",
  Family = "f__",
  Genus = "g__"
)

#------------------------------------------------------------
# Function: core taxa at one taxonomic rank
#------------------------------------------------------------
get_core_taxa <- function(dataset, rank, prevalence_cutoff = 0.80, detection_cutoff = 0){
  
  ra <- dataset$taxa_abund[[rank]]
  prefix <- rank_prefix[[rank]]
  
  long <- ra %>%
    as.data.frame() %>%
    tibble::rownames_to_column("Taxonomy") %>%
    tidyr::pivot_longer(
      cols = -Taxonomy, 
      names_to = "Sample", 
      values_to = "RA"
    ) %>%
    dplyr::mutate(
      Taxon = stringr::str_match(Taxonomy, paste0("(?:^|\\|)", prefix, "([^|]*)"))[, 2]
    ) %>%
    dplyr::filter(!is.na(Taxon), Taxon != "") %>%
    dplyr::group_by(Taxon, Sample) %>%
    dplyr::summarise(RA = sum(RA), .groups = "drop")
  
  # Verify metadata has Sample column before joining
  if (!"Sample" %in% colnames(metadata)) {
    metadata <- metadata %>% tibble::rownames_to_column("Sample")
  }
  
  long <- long %>%
    dplyr::left_join(metadata, by = "Sample")
  
  stats <- long %>%
    dplyr::group_by(Season, Taxon) %>%
    dplyr::summarise(
      n_samples = dplyr::n_distinct(Sample),
      n_detected = sum(RA > detection_cutoff),
      prevalence = n_detected / n_samples,
      mean_RA = mean(RA),
      median_RA = median(RA),
      .groups = "drop"
    ) %>%
    dplyr::mutate(Core = prevalence >= prevalence_cutoff)
  
  spring <- stats %>%
    dplyr::filter(Season == "Spring", Core) %>%
    dplyr::pull(Taxon)
  
  fall <- stats %>%
    dplyr::filter(Season == "Fall", Core) %>%
    dplyr::pull(Taxon)
  
  shared <- intersect(spring, fall)
  spring_only <- setdiff(spring, fall)
  fall_only <- setdiff(fall, spring)
  
  list(
    rank = rank,
    long = long,
    stats = stats,
    spring = spring,
    fall = fall,
    shared = shared,
    spring_only = spring_only,
    fall_only = fall_only
  )
}

core_phylum <- get_core_taxa(meco_dataset, "Phylum")
core_class <- get_core_taxa(meco_dataset, "Class")
core_order <- get_core_taxa(meco_dataset, "Order")
core_family <- get_core_taxa(meco_dataset, "Family")
core_genus <- get_core_taxa(meco_dataset, "Genus")

core_summary <- tibble(
  Rank = c("Phylum", "Class", "Order", "Family", "Genus"),
  Fall_core = c(length(core_phylum$fall),
                length(core_class$fall),
                length(core_order$fall),
                length(core_family$fall),
                length(core_genus$fall)),
  Spring_core = c(length(core_phylum$spring),
                  length(core_class$spring),
                  length(core_order$spring),
                  length(core_family$spring),
                  length(core_genus$spring)),
  Shared_core = c(length(core_phylum$shared),
                  length(core_class$shared),
                  length(core_order$shared),
                  length(core_family$shared),
                  length(core_genus$shared)),
  Fall_only = c(length(core_phylum$fall_only),
                length(core_class$fall_only),
                length(core_order$fall_only),
                length(core_family$fall_only),
                length(core_genus$fall_only)),
  Spring_only = c(length(core_phylum$spring_only),
                  length(core_class$spring_only),
                  length(core_order$spring_only),
                  length(core_family$spring_only),
                  length(core_genus$spring_only))
)

core_summary

write_csv(
  core_summary,
  "bacteria/Output/Core_microbiome/Core_taxa_summary_all_levels.csv"
)

plot_core_heatmap <- function(core_object){
  
  core_taxa <- union(core_object$fall, core_object$spring)
  
  tax_order <- core_object$stats %>%
    dplyr::filter(Taxon %in% core_taxa) %>%
    dplyr::group_by(Taxon) %>%
    dplyr::summarise(mean_RA = mean(mean_RA, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(mean_RA) %>%
    dplyr::pull(Taxon)
  
  sample_order <- metadata %>%
    dplyr::arrange(Season, farm, replication) %>%
    dplyr::pull(Sample)
  
  core_object$long %>%
    dplyr::filter(Taxon %in% core_taxa) %>%
    dplyr::mutate(
      Taxon = factor(Taxon, levels = tax_order),
      Sample = factor(Sample, levels = sample_order),
      RA_percent = RA * 100
    ) %>%
    ggplot(aes(x = Sample, y = Taxon, fill = RA_percent)) +
    geom_tile(color = "white", linewidth = 0.15) +
    facet_grid(~Season, scales = "free_x", space = "free_x") +
    scale_fill_viridis_c(option = "C", trans = "sqrt", name = "Relative\nabundance (%)") +
    labs(x = "Soil sample", y = paste("Core", tolower(core_object$rank))) +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
      axis.text.y = element_text(size = 9),
      strip.text = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6)
    )
}

core_heatmap_phylum <- plot_core_heatmap(core_phylum)
core_heatmap_class <- plot_core_heatmap(core_class)
core_heatmap_order <- plot_core_heatmap(core_order)
core_heatmap_family <- plot_core_heatmap(core_family)
core_heatmap_genus <- plot_core_heatmap(core_genus)

core_heatmap_phylum
core_heatmap_class
core_heatmap_order
core_heatmap_family
core_heatmap_genus

ggsave("bacteria/Output/Core_microbiome/Core_phylum_heatmap.pdf",
       core_heatmap_phylum, width = 10, height = 6)

ggsave("bacteria/Output/Core_microbiome/Core_class_heatmap.pdf",
       core_heatmap_class, width = 10, height = 7)

ggsave("bacteria/Output/Core_microbiome/Core_order_heatmap.pdf",
       core_heatmap_order, width = 10, height = 8)

ggsave("bacteria/Output/Core_microbiome/Core_family_heatmap.pdf",
       core_heatmap_family, width = 10, height = 9)

ggsave("bacteria/Output/Core_microbiome/Core_genus_heatmap.pdf",
       core_heatmap_genus, width = 10, height = 10)

#------------------------------------------------------------
# Output folder
#------------------------------------------------------------
outdir <- "bacteria/Output/Beta_diversity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

#------------------------------------------------------------
# Same colors used for alpha diversity
#------------------------------------------------------------
year_colors <- c("2020" = "#0072B2", "2021" = "#E69F00")
tillage_colors <- c("Min Till" = "#009E73", "Conv Till" = "#D55E00")
covercrop_colors <- c("CC-mix" = "#6A3D9A", "no CC" = "#FDBF6F")
rotation_colors <- c("Cotton-Soybean" = "#1B9E77", "Corn-Soybean" = "#D95F02",
                     "Soybean-Corn" = "#7570B3", "Soybean-Soybean" = "#E7298A")
season_colors <- c("Fall" = "#A6761D", "Spring" = "#66A61E")

#------------------------------------------------------------
# Metadata
#------------------------------------------------------------
metadata <- meco_dataset$sample_table %>%
  rownames_to_column("Sample") %>%
  mutate(
    Year = factor(Year, levels = c("2020", "2021")),
    tillage = factor(tillage, levels = c("Min Till", "Conv Till")),
    cover_crop = factor(cover_crop, levels = c("CC-mix", "no CC")),
    Crop_rotation = factor(Crop_rotation,
                           levels = c("Cotton-Soybean", "Corn-Soybean",
                                      "Soybean-Corn", "Soybean-Soybean")),
    Season = factor(Season, levels = c("Fall", "Spring"))
  )

#------------------------------------------------------------
# Bray-Curtis
#------------------------------------------------------------
otu <- as.data.frame(meco_dataset$otu_table)

bray <- vegan::vegdist(t(otu), method = "bray")

#------------------------------------------------------------
# PCoA with Lingoes correction
#------------------------------------------------------------
pcoa_res <- vegan::wcmdscale(bray, eig = TRUE, add = "lingoes", k = 2)

pcoa_dat <- as.data.frame(pcoa_res$points) %>%
  rownames_to_column("Sample") %>%
  dplyr::rename(PCoA1 = Dim1, PCoA2 = Dim2) %>%
  left_join(metadata, by = "Sample")

# Variance explained
eig <- pcoa_res$eig[pcoa_res$eig > 0]

pcoa1_var <- round(100 * eig[1] / sum(eig), 1)
pcoa2_var <- round(100 * eig[2] / sum(eig), 1)

pcoa1_var
pcoa2_var

write_csv(
  pcoa_dat,
  file.path(outdir, "Bray_Curtis_PCoA_coordinates.csv")
)

plot_pcoa <- function(data, factor_name, colors){
  ggplot(data, aes(x = PCoA1, y = PCoA2, fill = .data[[factor_name]])) +
    geom_point(shape = 21, size = 4.2, alpha = 0.9, color = "black", stroke = 0.6) +
    stat_ellipse(aes(group = .data[[factor_name]]),
                 type = "t", level = 0.95, linewidth = 0.9,
                 alpha = 0.35, geom = "polygon", show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    labs(
      x = paste0("PCoA1 (", pcoa1_var, "%)"),
      y = paste0("PCoA2 (", pcoa2_var, "%)"),
      fill = str_replace_all(factor_name, "_", " ")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}

#------------------------------------------------------------
# SINGLE FACTORS
#------------------------------------------------------------
pcoa_year <- plot_pcoa(pcoa_dat, "Year", year_colors)
pcoa_tillage <- plot_pcoa(pcoa_dat, "tillage", tillage_colors)
pcoa_covercrop <- plot_pcoa(pcoa_dat, "cover_crop", covercrop_colors)
pcoa_rotation <- plot_pcoa(pcoa_dat, "Crop_rotation", rotation_colors)
pcoa_season <- plot_pcoa(pcoa_dat, "Season", season_colors)

pcoa_year
pcoa_tillage
pcoa_covercrop
pcoa_rotation
pcoa_season

plot_pcoa_2factor <- function(data, factor1, factor2, colors){
  
  ggplot(data,
         aes(x = PCoA1,
             y = PCoA2,
             color = .data[[factor2]],
             shape = .data[[factor1]])) +
    geom_point(size = 3.8, alpha = 0.9) +
    stat_ellipse(aes(color = .data[[factor2]], fill = .data[[factor2]], group = .data[[factor2]]),
                 geom = "polygon", type = "t", level = 0.95,
                 alpha = 0.35, linewidth = 0.8, show.legend = FALSE) +
    scale_color_manual(values = colors, drop = FALSE) +
    labs(
      x = paste0("PCoA1 (", pcoa1_var, "%)"),
      y = paste0("PCoA2 (", pcoa2_var, "%)"),
      color = str_replace_all(factor2, "_", " "),
      shape = str_replace_all(factor1, "_", " ")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}

#------------------------------------------------------------
# YEAR × other factors
#------------------------------------------------------------
pcoa_year_tillage <- plot_pcoa_2factor(pcoa_dat, "Year", "tillage", tillage_colors)
pcoa_year_covercrop <- plot_pcoa_2factor(pcoa_dat, "Year", "cover_crop", covercrop_colors)
pcoa_year_rotation <- plot_pcoa_2factor(pcoa_dat, "Year", "Crop_rotation", rotation_colors)
pcoa_year_season <- plot_pcoa_2factor(pcoa_dat, "Year", "Season", season_colors)

#------------------------------------------------------------
# TILLAGE × other factors
#------------------------------------------------------------
pcoa_tillage_covercrop <- plot_pcoa_2factor(pcoa_dat, "tillage", "cover_crop", covercrop_colors)
pcoa_tillage_rotation <- plot_pcoa_2factor(pcoa_dat, "tillage", "Crop_rotation", rotation_colors)
pcoa_tillage_season <- plot_pcoa_2factor(pcoa_dat, "tillage", "Season", season_colors)

#------------------------------------------------------------
# COVER CROP × other factors
#------------------------------------------------------------
pcoa_covercrop_rotation <- plot_pcoa_2factor(pcoa_dat, "cover_crop", "Crop_rotation", rotation_colors)
pcoa_covercrop_season <- plot_pcoa_2factor(pcoa_dat, "cover_crop", "Season", season_colors)

#------------------------------------------------------------
# ROTATION × SEASON
#------------------------------------------------------------
pcoa_rotation_season <- plot_pcoa_2factor(pcoa_dat, "Crop_rotation", "Season", season_colors)

pcoa_year_tillage
pcoa_year_covercrop
pcoa_year_rotation
pcoa_tillage_rotation
pcoa_covercrop_rotation
pcoa_rotation_season

plot_pcoa_3factor <- function(data, factor1, factor2, facet_factor, colors){
  
  ggplot(data,
         aes(x = PCoA1,
             y = PCoA2,
             color = .data[[factor2]],
             shape = .data[[factor1]])) +
    geom_point(size = 3.8, alpha = 0.9) +
    stat_ellipse(aes(color = .data[[factor2]], fill = .data[[factor2]], group = .data[[factor2]]),
                 geom = "polygon", type = "t", level = 0.95,
                 alpha = 0.30, linewidth = 0.8, show.legend = FALSE) +
    facet_wrap(as.formula(paste("~", facet_factor))) +
    scale_color_manual(values = colors, drop = FALSE) +
    labs(
      x = paste0("PCoA1 (", pcoa1_var, "%)"),
      y = paste0("PCoA2 (", pcoa2_var, "%)"),
      color = str_replace_all(factor2, "_", " "),
      shape = str_replace_all(factor1, "_", " ")
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 12),
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      strip.text = element_text(size = 12, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8)
    )
}


#------------------------------------------------------------
# THREE FACTORS containing YEAR
#------------------------------------------------------------

pcoa_tillage_covercrop_year <- plot_pcoa_3factor(
  pcoa_dat, "tillage", "cover_crop", "Year", covercrop_colors)

pcoa_tillage_rotation_year <- plot_pcoa_3factor(
  pcoa_dat, "tillage", "Crop_rotation", "Year", rotation_colors)

pcoa_tillage_season_year <- plot_pcoa_3factor(
  pcoa_dat, "tillage", "Season", "Year", season_colors)

pcoa_covercrop_rotation_year <- plot_pcoa_3factor(
  pcoa_dat, "cover_crop", "Crop_rotation", "Year", rotation_colors)

pcoa_covercrop_season_year <- plot_pcoa_3factor(
  pcoa_dat, "cover_crop", "Season", "Year", season_colors)

pcoa_rotation_season_year <- plot_pcoa_3factor(
  pcoa_dat, "Crop_rotation", "Season", "Year", season_colors)

#------------------------------------------------------------
# THREE FACTORS without YEAR
#------------------------------------------------------------

pcoa_tillage_rotation_covercrop <- plot_pcoa_3factor(
  pcoa_dat, "tillage", "Crop_rotation", "cover_crop", rotation_colors)

pcoa_tillage_covercrop_season <- plot_pcoa_3factor(
  pcoa_dat, "tillage", "cover_crop", "Season", covercrop_colors)

pcoa_tillage_rotation_season <- plot_pcoa_3factor(
  pcoa_dat, "tillage", "Crop_rotation", "Season", rotation_colors)

pcoa_covercrop_rotation_season <- plot_pcoa_3factor(
  pcoa_dat, "cover_crop", "Crop_rotation", "Season", rotation_colors)

pcoa_covercrop_rotation_year
pcoa_tillage_rotation_year

single_pcoa <- list(
  pcoa_year = pcoa_year,
  pcoa_tillage = pcoa_tillage,
  pcoa_covercrop = pcoa_covercrop,
  pcoa_rotation = pcoa_rotation,
  pcoa_season = pcoa_season
)

twofactor_pcoa <- list(
  pcoa_year_tillage = pcoa_year_tillage,
  pcoa_year_covercrop = pcoa_year_covercrop,
  pcoa_year_rotation = pcoa_year_rotation,
  pcoa_year_season = pcoa_year_season,
  pcoa_tillage_covercrop = pcoa_tillage_covercrop,
  pcoa_tillage_rotation = pcoa_tillage_rotation,
  pcoa_tillage_season = pcoa_tillage_season,
  pcoa_covercrop_rotation = pcoa_covercrop_rotation,
  pcoa_covercrop_season = pcoa_covercrop_season,
  pcoa_rotation_season = pcoa_rotation_season
)

threefactor_pcoa <- list(
  pcoa_tillage_covercrop_year = pcoa_tillage_covercrop_year,
  pcoa_tillage_rotation_year = pcoa_tillage_rotation_year,
  pcoa_tillage_season_year = pcoa_tillage_season_year,
  pcoa_covercrop_rotation_year = pcoa_covercrop_rotation_year,
  pcoa_covercrop_season_year = pcoa_covercrop_season_year,
  pcoa_rotation_season_year = pcoa_rotation_season_year,
  pcoa_tillage_rotation_covercrop = pcoa_tillage_rotation_covercrop,
  pcoa_tillage_covercrop_season = pcoa_tillage_covercrop_season,
  pcoa_tillage_rotation_season = pcoa_tillage_rotation_season,
  pcoa_covercrop_rotation_season = pcoa_covercrop_rotation_season
)

all_pcoa <- c(single_pcoa, twofactor_pcoa, threefactor_pcoa)

iwalk(
  all_pcoa,
  ~ggsave(
    filename = file.path(outdir, paste0(.y, ".pdf")),
    plot = .x,
    width = 8,
    height = 6
  )
)

#------------------------------------------------------------
# Metadata aligned exactly with Bray-Curtis matrix
#------------------------------------------------------------
metadata_perm <- meco_dataset$sample_table %>%
  tibble::rownames_to_column("Sample") %>%
  dplyr::mutate(
    Year = factor(Year, levels = c("2020", "2021")),
    tillage = factor(tillage, levels = c("Min Till", "Conv Till")),
    cover_crop = factor(cover_crop, levels = c("CC-mix", "no CC")),
    Crop_rotation = factor(Crop_rotation, levels = c("Cotton-Soybean", "Corn-Soybean", "Soybean-Corn", "Soybean-Soybean")),
    Season = factor(Season, levels = c("Fall", "Spring"))
  ) %>%
  dplyr::slice(match(attr(bray, "Labels"), Sample))

stopifnot(all(metadata_perm$Sample == attr(bray, "Labels")))

#------------------------------------------------------------
# PERMANOVA function
#------------------------------------------------------------
run_permanova <- function(rhs, model_name){
  fit <- adonis2(as.formula(paste("bray ~", rhs)), data = metadata_perm,
                 permutations = 9999, method = "bray", by = "margin")
  as.data.frame(fit) %>%
    rownames_to_column("Term") %>%
    filter(Term != "Residual", Term != "Total") %>%
    mutate(Model = model_name, .before = 1)
}

perm_year <- run_permanova("Year", "Year")
perm_tillage <- run_permanova("tillage", "Tillage")
perm_covercrop <- run_permanova("cover_crop", "Cover crop")
perm_rotation <- run_permanova("Crop_rotation", "Crop rotation")
perm_season <- run_permanova("Season", "Season")

permanova_single <- bind_rows(
  perm_year,
  perm_tillage,
  perm_covercrop,
  perm_rotation,
  perm_season
)

permanova_single

write_csv(
  permanova_single,
  file.path(outdir, "PERMANOVA_single_factors.csv")
)

perm_year_tillage <- run_permanova("Year * tillage", "Year × Tillage")
perm_year_covercrop <- run_permanova("Year * cover_crop", "Year × Cover crop")
perm_year_rotation <- run_permanova("Year * Crop_rotation", "Year × Crop rotation")

perm_tillage_rotation <- run_permanova("tillage * Crop_rotation", "Tillage × Crop rotation")
perm_tillage_season <- run_permanova("tillage * Season", "Tillage × Season")

perm_covercrop_rotation <- run_permanova("cover_crop * Crop_rotation", "Cover crop × Crop rotation")
perm_covercrop_season <- run_permanova("cover_crop * Season", "Cover crop × Season")

perm_rotation_season <- run_permanova("Crop_rotation * Season", "Crop rotation × Season")

permanova_twofactor <- bind_rows(
  perm_year_tillage,
  perm_year_covercrop,
  perm_year_rotation,
  perm_tillage_rotation,
  perm_tillage_season,
  perm_covercrop_rotation,
  perm_covercrop_season,
  perm_rotation_season
)

permanova_twofactor

write_csv(
  permanova_twofactor,
  file.path(outdir, "PERMANOVA_two_factor_interactions.csv")
)

perm_year_tillage_rotation <- run_permanova(
  "Year * tillage * Crop_rotation",
  "Year × Tillage × Crop rotation"
)

perm_year_covercrop_rotation <- run_permanova(
  "Year * cover_crop * Crop_rotation",
  "Year × Cover crop × Crop rotation"
)

perm_season_tillage_rotation <- run_permanova(
  "Season * tillage * Crop_rotation",
  "Season × Tillage × Crop rotation"
)

perm_season_covercrop_rotation <- run_permanova(
  "Season * cover_crop * Crop_rotation",
  "Season × Cover crop × Crop rotation"
)

permanova_threefactor <- bind_rows(
  perm_year_tillage_rotation,
  perm_year_covercrop_rotation,
  perm_season_tillage_rotation,
  perm_season_covercrop_rotation
)

permanova_threefactor

write_csv(
  permanova_threefactor,
  file.path(outdir, "PERMANOVA_three_factor_interactions.csv")
)


