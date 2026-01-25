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

biom = import_biom("Inputfiles/season-its.biom")

metadata = import_qiime_sample_data("Inputfiles/metadata.txt")

#tree = read_tree("Input_files/rooted_tree.nwk")

#rep_fasta = readDNAStringSet("Input_files/fun-seq.fasta", format = "fasta")

season_biom = merge_phyloseq(biom, metadata)

colnames(tax_table(season_biom)) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

meco_dataset <- phyloseq2meco(season_biom)

soil <- read.csv("Inputfiles/soil_data.csv")

rownames(soil) <- soil[, 1]

soil = soil[ ,-1]

# add_data is used to add the environmental data
soil_season <- trans_env$new(dataset = meco_dataset, add_data = soil)

meco_dataset$cal_abund()
meco_dataset$cal_alphadiv()
meco_dataset$cal_betadiv()

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

phylum_fig_heat <- phylum_heat$plot_heatmap(facet = "Season", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10))+ theme(axis.text.y = element_text(face = 'italic'))

genus_fig_heat <- genus_heat$plot_heatmap(facet = "Season", xtext_keep = FALSE, withmargin = FALSE, plot_breaks = c(0.01, 0.1, 1, 10)) + theme(axis.text.y = element_text(face = 'italic'))

genus_fig_heat
phylum_fig_heat

ggsave("Output/PDFs/genus_fig_heat.pdf", plot = genus_fig_heat, width =8, height = 10, dpi = 1000)
ggsave("Output/PDFs/phylum_fig_heat.pdf", plot = phylum_fig_heat, width =8, height = 10, dpi = 1000)

### Alpha diversity analaysis

install.packages("agricolae")

library(agricolae)

alpha_S <- trans_alpha$new(dataset = meco_dataset, group = "Season")

alpha_S$cal_diff(method = "t.test")

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

plot_alpha_season

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

df_abundance$plot_diff_bar(threshold = 5)
# we show 20 taxa with the highest LDA (log10)

df_bar_plot = df_abundance$plot_diff_bar(use_number = 1:30, width = 0.8)

ggsave("Output/PDFs/df_bar_plot.pdf", plot = df_bar_plot, width =7, height = 8, dpi = 1000)


df_abun_plot = df_abundance$plot_diff_abund(plot_type = "barerrorbar", use_number = 1:30, width = 0.8)

df_abun_plot
ggsave("Output/PDFs/df_abun_plot.pdf", plot = df_abun_plot, width =7, height = 9, dpi = 1000)

library(microeco)
library(ggplot2)
library(grid)

## 1) Merge samples by Season (this keeps OTUs as rows)
venn_mt <- meco_dataset$merge_samples("Season")

## 2) Build venn object
## For "OTU venn" (presence/feature counts), numratio is usually what you want
Venn_plot <- trans_venn$new(dataset = venn_mt, ratio = "numratio")

## 3) Patch microeco's internal theme to satisfy ggplot2's axis.title requirements
Venn_plot$.__enclos_env__$private$main_theme <- theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  legend.key = element_blank(),
  plot.margin = unit(c(0, 0, 0, 0), "mm"),
  
  # ---- critical fix ----
  axis.title = element_text(),     # must be element_text for new ggplot2
  axis.title.x = element_blank(),  # keep titles hidden
  axis.title.y = element_blank()
)

## 4) Plot
venn_fig <- Venn_plot$plot_venn(
  fill_color = TRUE,
  text_size = 5,
  text_name_size = 6,
  alpha = 0.35,
  linesize = 1.1
)

venn_fig

ggsave("Output/PDFs/venn_otus_season.pdf", plot = venn_fig, width = 7, height = 6, dpi = 1000)


