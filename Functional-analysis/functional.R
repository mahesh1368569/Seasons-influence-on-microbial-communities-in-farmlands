library(tidyverse)
library(vegan)
library(pheatmap)
library(ggpubr)

func <- read_tsv("pred_metagenome_unstrat_descrip.tsv")
meta <- read_tsv("metadata.txt") %>% rename(SampleID = `#Name`)

func_long <- func %>%
  pivot_longer(-c(function, description), names_to = "SampleID", values_to = "Abundance") %>%
  inner_join(meta, by = "SampleID") %>%
  mutate(Season = factor(Season, levels = c("Fall_2020", "Spring_2021")),
         cover_crop = factor(cover_crop, levels = c("no CC", "CC-mix")),
         tillage = factor(tillage, levels = c("Conv Till", "Min Till")))