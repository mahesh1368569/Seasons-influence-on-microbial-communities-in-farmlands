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

########################################################################################
## Alpha Diversity 
########################################################################################

#------------------------------------------------------------
# 1. Check metadata
#------------------------------------------------------------

meco_dataset$sample_table %>% head()
meco_dataset$sample_table %>% glimpse()

# Set factor levels
meco_dataset$sample_table <- meco_dataset$sample_table %>%
  mutate(
    Year = factor(Year, levels = c("2020", "2021")),
    farm = factor(farm, levels = c("SCH", "ROB", "KIN", "EVA", "OGL", "BRE")),
    tillage = factor(tillage, levels = c("Min Till", "Conv Till")),
    cover_crop = factor(cover_crop, levels = c("CC-mix", "no CC")),
    Crop_rotation = factor(Crop_rotation,
                           levels = c("Cotton-Soybean", "Corn-Soybean",
                                      "Soybean-Corn", "Soybean-Soybean")),
    Season = factor(Season, levels = c("Fall", "Spring"))
  )

#------------------------------------------------------------
# 2. Colors for each factor
#------------------------------------------------------------

year_colors <- c(
  "2020" = "#0072B2",
  "2021" = "#E69F00"
)

farm_colors <- c(
  "SCH" = "#0072B2",
  "ROB" = "#E69F00",
  "KIN" = "#009E73",
  "EVA" = "#CC79A7",
  "OGL" = "#D55E00",
  "BRE" = "#56B4E9"
)

tillage_colors <- c(
  "Min Till" = "#009E73",
  "Conv Till" = "#D55E00"
)

covercrop_colors <- c(
  "CC-mix" = "#56B4E9",
  "no CC" = "#E69F00"
)

rotation_colors <- c(
  "Cotton-Soybean" = "#0072B2",
  "Corn-Soybean" = "#009E73",
  "Soybean-Corn" = "#E69F00",
  "Soybean-Soybean" = "#CC79A7"
)

season_colors <- c(
  "Fall" = "#D55E00",
  "Spring" = "#56B4E9"
)

factor_colors <- list(
  Year = year_colors,
  farm = farm_colors,
  tillage = tillage_colors,
  cover_crop = covercrop_colors,
  Crop_rotation = rotation_colors,
  Season = season_colors
)

#------------------------------------------------------------
# 3. Calculate/extract alpha diversity
#------------------------------------------------------------

alpha_obj <- trans_alpha$new(dataset = meco_dataset)

alpha_dat <- alpha_obj$data_alpha %>%
  as_tibble() %>%
  mutate(
    Year = factor(Year, levels = c("2020", "2021")),
    farm = factor(farm, levels = c("SCH", "ROB", "KIN", "EVA", "OGL", "BRE")),
    tillage = factor(tillage, levels = c("Min Till", "Conv Till")),
    cover_crop = factor(cover_crop, levels = c("CC-mix", "no CC")),
    Crop_rotation = factor(Crop_rotation, levels = c("Cotton-Soybean", "Corn-Soybean", "Soybean-Corn", "Soybean-Soybean")),
    Season = factor(Season, levels = c("Fall", "Spring"))
  )

library(tidyverse)
library(ggplot2)

outdir <- "bacteria/Output/Alpha_diversity"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

year_colors <- c("2020" = "#0072B2", "2021" = "#E69F00")
tillage_colors <- c("Min Till" = "#009E73", "Conv Till" = "#D55E00")
covercrop_colors <- c("CC-mix" = "#6A3D9A", "no CC" = "#FDBF6F")
rotation_colors <- c("Cotton-Soybean" = "#1B9E77", "Corn-Soybean" = "#D95F02",
                     "Soybean-Corn" = "#7570B3", "Soybean-Soybean" = "#E7298A")
season_colors <- c("Fall" = "#A6761D", "Spring" = "#66A61E")

plot_alpha <- function(data, metric, factor_name, colors){
  data %>% filter(Measure == metric) %>%
    ggplot(aes(x = .data[[factor_name]], y = Value, fill = .data[[factor_name]])) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 2) +
    scale_fill_manual(values = colors) +
    labs(x = str_replace_all(factor_name, "_", " "), y = metric,
         fill = str_replace_all(factor_name, "_", " ")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
}

plot_alpha_2factor <- function(data, metric, factor1, factor2, colors){
  data %>% filter(Measure == metric) %>%
    ggplot(aes(x = .data[[factor1]], y = Value, fill = .data[[factor2]])) +
    geom_boxplot(aes(group = interaction(.data[[factor1]], .data[[factor2]])),
                 position = position_dodge(width = 0.8), width = 0.65,
                 outlier.shape = NA, alpha = 0.85) +
    geom_point(aes(group = .data[[factor2]]),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
               alpha = 0.6, size = 2) +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(x = str_replace_all(factor1, "_", " "), y = metric,
         fill = str_replace_all(factor2, "_", " ")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
}

plot_alpha_3factor <- function(data, metric, factor1, factor2, facet_factor, colors){
  data %>% filter(Measure == metric) %>%
    ggplot(aes(x = .data[[factor1]], y = Value, fill = .data[[factor2]])) +
    geom_boxplot(aes(group = interaction(.data[[factor1]], .data[[factor2]])),
                 position = position_dodge(width = 0.8), width = 0.65,
                 outlier.shape = NA, alpha = 0.85) +
    geom_point(aes(group = .data[[factor2]]),
               position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.8),
               alpha = 0.6, size = 2) +
    facet_wrap(as.formula(paste("~", facet_factor))) +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(x = str_replace_all(factor1, "_", " "), y = metric,
         fill = str_replace_all(factor2, "_", " ")) +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
          axis.text.y = element_text(size = 11),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 11),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 12),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8))
}

make_alpha_plots <- function(metric){
  
  plots <- list(
    
    # Single factors
    year = plot_alpha(alpha_dat, metric, "Year", year_colors),
    tillage = plot_alpha(alpha_dat, metric, "tillage", tillage_colors),
    covercrop = plot_alpha(alpha_dat, metric, "cover_crop", covercrop_colors),
    rotation = plot_alpha(alpha_dat, metric, "Crop_rotation", rotation_colors),
    season = plot_alpha(alpha_dat, metric, "Season", season_colors),
    
    # Two factors
    year_tillage = plot_alpha_2factor(alpha_dat, metric, "Year", "tillage", tillage_colors),
    year_covercrop = plot_alpha_2factor(alpha_dat, metric, "Year", "cover_crop", covercrop_colors),
    year_rotation = plot_alpha_2factor(alpha_dat, metric, "Year", "Crop_rotation", rotation_colors),
    year_season = plot_alpha_2factor(alpha_dat, metric, "Year", "Season", season_colors),
    
    tillage_covercrop = plot_alpha_2factor(alpha_dat, metric, "tillage", "cover_crop", covercrop_colors),
    tillage_rotation = plot_alpha_2factor(alpha_dat, metric, "tillage", "Crop_rotation", rotation_colors),
    tillage_season = plot_alpha_2factor(alpha_dat, metric, "tillage", "Season", season_colors),
    
    covercrop_rotation = plot_alpha_2factor(alpha_dat, metric, "cover_crop", "Crop_rotation", rotation_colors),
    covercrop_season = plot_alpha_2factor(alpha_dat, metric, "cover_crop", "Season", season_colors),
    
    rotation_season = plot_alpha_2factor(alpha_dat, metric, "Crop_rotation", "Season", season_colors),
    
    # Three factors
    covercrop_rotation_year = plot_alpha_3factor(
      alpha_dat, metric, "cover_crop", "Crop_rotation", "Year", rotation_colors),
    
    tillage_rotation_year = plot_alpha_3factor(
      alpha_dat, metric, "tillage", "Crop_rotation", "Year", rotation_colors),
    
    rotation_covercrop_year = plot_alpha_3factor(
      alpha_dat, metric, "Crop_rotation", "cover_crop", "Year", covercrop_colors),
    
    rotation_tillage_year = plot_alpha_3factor(
      alpha_dat, metric, "Crop_rotation", "tillage", "Year", tillage_colors),
    
    covercrop_rotation_season = plot_alpha_3factor(
      alpha_dat, metric, "cover_crop", "Crop_rotation", "Season", rotation_colors),
    
    tillage_rotation_season = plot_alpha_3factor(
      alpha_dat, metric, "tillage", "Crop_rotation", "Season", rotation_colors)
  )
  
  metric_dir <- file.path(outdir, metric)
  dir.create(metric_dir, recursive = TRUE, showWarnings = FALSE)
  
  iwalk(plots, ~ggsave(
    file.path(metric_dir, paste0(tolower(metric), "_", .y, ".pdf")),
    plot = .x, width = 8, height = 6
  ))
  
  pdf(file.path(outdir, paste0(metric, "_all_alpha_diversity_plots.pdf")),
      width = 8, height = 6)
  walk(plots, print)
  dev.off()
  
  plots
}


#------------------------------------------------------------
# Generate Shannon and Chao1
#------------------------------------------------------------

shannon_plots <- make_alpha_plots("Shannon")
chao1_plots <- make_alpha_plots("Chao1")
simpson_plots <- make_alpha_plots("Simpson")

write_csv(alpha_dat, "bacteria/Output/Alpha_diversity/Alpha_diversity_long.csv")


##Abundance anaysis
abun = trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 15, groupmean = "Season")

plotbar_season <- abun$plot_bar(others_color = "grey70", legend_text_italic = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 15), # Increase x-axis text size
        axis.text.y = element_text(size = 15), # Increase y-axis text size
        axis.title.x = element_text(size = 15), # Increase x-axis label size
        axis.title.y = element_text(size = 15), # Increase y-axis label size
        strip.text = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

plotbar_season

ggsave("Output/PDFs/bar_fig.pdf", plot = plotbar_season, width =7, height = 10, dpi = 1000)

# show 15 taxa at Class level
abun_bar <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 15)

boxplot_season = abun_bar$plot_box(group = "Season", xtext_angle = 30)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1)) # Add border

boxplot_season

ggsave("Output/PDFs/box_fig.pdf", plot = boxplot_season, width =11, height = 8, dpi = 1000)

# show 40 taxa at Genus level
genus_heat <- trans_abund$new(dataset = meco_dataset, taxrank = "Genus", ntaxa = 40)
phylum_heat <- trans_abund$new(dataset = meco_dataset, taxrank = "Phylum", ntaxa = 40)

phylum_fig_heat <- phylum_heat$plot_heatmap(facet = "Season", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
+ theme(axis.text.y = element_text(face = 'italic'))

genus_fig_heat <- genus_heat$plot_heatmap(facet = "Season", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))
+ theme(axis.text.y = element_text(face = 'italic'))

genus_fig_heat
phylum_fig_heat

ggsave("Output/PDFs/genus_fig_heat.pdf", plot = genus_fig_heat, width =8, height = 10, dpi = 1000)
ggsave("Output/PDFs/phylum_fig_heat.pdf", plot = phylum_fig_heat, width =8, height = 10, dpi = 1000)

### Alpha diversity analaysis

alpha_S <- trans_alpha$new(dataset = meco_dataset, group = "Season")

alpha_S$cal_diff(method = "wilcox")

plot_alpha_season  = alpha_S$plot_alpha(measure = "Chao1", add = 'jitter', shape = "Season")+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave("Output/PDFs/plot_alpha_season.pdf", plot = plot_alpha_season, width =6, height = 5, dpi = 1000)

# Beta diversity #######

beta_seas <- trans_beta$new(dataset = meco_dataset, group = "Season", measure = "bray")

beta_seas$cal_ordination(method = "PCoA")


beta_plot = beta_seas$plot_ordination(plot_color = "Season", plot_shape = "Crop_rotation", plot_type = c("point", "ellipse"))+
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14), # Increase x-axis text size
        axis.text.y = element_text(size = 14), # Increase y-axis text size
        axis.title.x = element_text(size = 14), # Increase x-axis label size
        axis.title.y = element_text(size = 14), # Increase y-axis label size
        strip.text = element_text(size = 14),
        legend.position = "right",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),# Increase facet label size
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave("Output/PDFs/beta_plot.pdf", plot = beta_plot, width =8, height = 6, dpi = 1000)

#### RDA and environmental data analysis

# use Genus
soil_season$cal_ordination(method = "RDA", taxa_level = "Genus")


soil_season$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

rdaplot = soil_season$plot_ordination(plot_color = "Season")

soil_season$cal_ordination_anova()

ggsave("Output/PDFs/rdaplot.pdf", plot = rdaplot, width =7, height = 6, dpi = 1000)

# use phylum level

soil_season$cal_ordination(method = "RDA", taxa_level = "Phylum")

soil_season$trans_ordination(show_taxa = 10, adjust_arrow_length = TRUE, max_perc_env = 1.5, max_perc_tax = 1.5, min_perc_env = 0.2, min_perc_tax = 0.2)

rdaplot = soil_season$plot_ordination(plot_color = "Season")

soil_season$cal_ordination_anova()

ggsave("Output/PDFs/phylum_rdaplot.pdf", plot = rdaplot, width =7, height = 6, dpi = 1000)

# correlation between taxa and soil properties

soil_season$cal_cor(use_data = "Genus", p_adjust_method = "fdr", p_adjust_type = "Env")

soil_season$res_cor

# filter genera that donot have at least one ***
soil_season$plot_cor(filter_feature = c("", "*", "**"))

# use pH and bray-curtis distance
# add correlation statistics
ph = soil_season$plot_scatterfit(
  x = "pH", 
  y = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)], 
  type = "cor",
  point_size = 3, point_alpha = 0.1, 
  label.x.npc = "center", label.y.npc = "bottom", 
  x_axis_title = "Euclidean distance of pH", 
  y_axis_title = "Bray-Curtis distance"
)

C = soil_season$plot_scatterfit(
  x = "C", 
  y = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)], 
  type = "cor",
  point_size = 3, point_alpha = 0.1, 
  label.x.npc = "center", label.y.npc = "bottom", 
  x_axis_title = "Euclidean distance of Carbon", 
  y_axis_title = "Bray-Curtis distance"
)

N = soil_season$plot_scatterfit(
  x = "N", 
  y = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)], 
  type = "cor",
  point_size = 3, point_alpha = 0.1, 
  label.x.npc = "center", label.y.npc = "bottom", 
  x_axis_title = "Euclidean distance of Nitrogen", 
  y_axis_title = "Bray-Curtis distance"
)

Glomalin = soil_season$plot_scatterfit(
  x = "Glomalin", 
  y = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)], 
  type = "cor",
  point_size = 3, point_alpha = 0.1, 
  label.x.npc = "center", label.y.npc = "bottom", 
  x_axis_title = "Euclidean distance of Glomalin", 
  y_axis_title = "Bray-Curtis distance"
)

ggsave("Output/PDFs/ph_vs_bray curtis.pdf", plot = ph, width =7, height = 6, dpi = 1000)
ggsave("Output/PDFs/C_vs_bray curtis.pdf", plot = C, width =7, height = 6, dpi = 1000)
ggsave("Output/PDFs/N_vs_bray curtis.pdf", plot = N, width =7, height = 6, dpi = 1000)
ggsave("Output/PDFs/Glomalin_vs_bray curtis.pdf", plot = Glomalin, width =7, height = 6, dpi = 1000)

# regression with type = "lm", use group parameter for different groups
ph_season = soil_season$plot_scatterfit(
  x = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)],
  y = "pH",
  type = "lm", 
  group = "Season",
  point_size = 3, point_alpha = 0.3, line_se = FALSE, line_size = 1.5, shape_values = c(16, 17, 7),
  y_axis_title = "Euclidean distance of pH", x_axis_title = "Bray-Curtis distance"
) + theme(axis.title = element_text(size = 17))

N_season = soil_season$plot_scatterfit(
  x = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)],
  y = "N",
  type = "lm", 
  group = "Season",
  point_size = 3, point_alpha = 0.3, line_se = FALSE, line_size = 1.5, shape_values = c(16, 17, 7),
  y_axis_title = "Euclidean distance of Nitrogen", x_axis_title = "Bray-Curtis distance"
) + theme(axis.title = element_text(size = 17))

C_season = soil_season$plot_scatterfit(
  x = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)],
  y = "C",
  type = "lm", 
  group = "Season",
  point_size = 3, point_alpha = 0.3, line_se = FALSE, line_size = 1.5, shape_values = c(16, 17, 7),
  y_axis_title = "Euclidean distance of Carbon", x_axis_title = "Bray-Curtis distance"
) + theme(axis.title = element_text(size = 17))

Glom_season = soil_season$plot_scatterfit(
  x = meco_dataset$beta_diversity$bray[rownames(soil_season$data_env), rownames(soil_season$data_env)],
  y = "Glomalin",
  type = "lm", 
  group = "Season",
  point_size = 3, point_alpha = 0.3, line_se = FALSE, line_size = 1.5, shape_values = c(16, 17, 7),
  y_axis_title = "Euclidean distance of Glomalin", x_axis_title = "Bray-Curtis distance"
) + theme(axis.title = element_text(size = 17))

ggsave("Output/PDFs/ph_season.pdf", plot = ph_season, width =7, height = 6, dpi = 1000)
ggsave("Output/PDFs/N_seaon.pdf", plot = N_season, width =7, height = 6, dpi = 1000)
ggsave("Output/PDFs/C_seaon.pdf", plot = C_season, width =7, height = 6, dpi = 1000)
ggsave("Output/PDFs/Glomalin_season.pdf", plot = Glom_season, width =7, height = 6, dpi = 1000)

###Venn diagram plots
venn <- meco_dataset$merge_samples("Season")
# tmp is a new microtable object
# create trans_venn object
Venn_plot <- trans_venn$new(venn, ratio = "seqratio")

Venn_plot$plot_venn()

#### Differential abundance analysis #####

df_abundance <- trans_diff$new(dataset = meco_dataset, method = "lefse", group = "Season", alpha = 0.01, lefse_subgroup = NULL)

df_abundance$plot_diff_bar(threshold = 6)
# we show 20 taxa with the highest LDA (log10)

df_bar_plot = df_abundance$plot_diff_bar(use_number = 1:30, width = 0.8)

ggsave("Output/PDFs/df_bar_plot.pdf", plot = df_bar_plot, width =7, height = 8, dpi = 1000)


df_abun_plot = df_abundance$plot_diff_abund(plot_type = "barerrorbar", use_number = 1:30, width = 0.8)

df_abun_plot
ggsave("Output/PDFs/df_abun_plot.pdf", plot = df_abun_plot, width =7, height = 9, dpi = 1000)
