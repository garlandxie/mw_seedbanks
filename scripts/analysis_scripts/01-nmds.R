################################################################################
# Accompanying code for the followng research project: 
#   Seedbanks in the Meadoway (Scarborough, Canada)
#
#
# Corresponding authors for this script:  
#   Garland Xie      (1)
#
# Affiliations: 
#   (1) Department of Biological Sciences, 
#       University of Toronto Scarborough,
#       1265 Military Trail, Toronto, ON, M1C 1A4, Canada
#       email: garland.xie@mail.utoronto.ca, 
#              nicholas.sookhan@mail.utoronto.ca
#              scott.macivor@mail.utoronto.ca
#
# Purpose of this R script: to perform non-metric dimensional scaling on the 
# seed bank datasets
#
# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html

# libraries ----
library(here)    # for creating relative file-paths
library(ggplot2) # for visualizing data 
library(dplyr)   # for manipulating data 
library(ggrepel) # for creating easy-to-read text labels in ggplot2
library(tidyr)   # for creating wide or long data frames
library(vegan)   # for performing non-metric dimensional scaling
library(tibble)  # for assigning row names for tibbles

# import data ----

spring_df <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"), 
  row.names = 1
)

# create community data matrix ----

spring_comm_matrix <- spring_df %>%
  pivot_wider(names_from = spp_code, values_from = total_abund) %>%
  ungroup() %>%
  mutate(across(everything(.), replace_na, 0)) %>%
  mutate(id = paste(site, plot, sep = "-")) %>%
  relocate(id)

# spring: non-metric dimensional scaling ----

nmds_sb <- spring_comm_matrix %>%
  select(-season, -section, -treatment, -site, -plot) %>%
  column_to_rownames(var = "id") %>%
  metaMDS(distnace = "bray")

# sanity check
stressplot(nmds_sb)

# plotting ----

# get NMDS plot scores into a data frame
plot_scores <- 
  nmds_sb %>%
  scores() %>%
  as.data.frame() %>%
  rownames_to_column(var = "id") %>%
  inner_join(
    spring_comm_matrix %>% select(section, treatment, site, id),
    by = "id")

# get NMDS species scores into a data frame
spp_scores <- 
  nmds_sb %>%
  scores("species") %>%
  as.data.frame() %>%
  rownames_to_column(var = "id")

# plot
nmds_ggplot <- ggplot() + 
  geom_text_repel(
    data = spp_scores,
    aes(x = NMDS1, y = NMDS2, label = id),
    size = 3,
    max.overlaps = 30,
    min.segment.length = 5) +  
  geom_point(
    data = plot_scores, 
    aes(x = NMDS1, y = NMDS2,
        col = treatment), 
    alpha = 0.7
  ) + 
  scale_color_manual(
    name = "Management Regime",
    labels = c("Maintenance Mow", "Undisturbed", "Seed Drill"),
    values = c("#F0E442", "#E69F00", "#56B4E9")
  ) + 
  theme_bw()

# save to disk -----

ggsave(
  filename = here("output", "results", "nmds_spring.png"), 
  plot = nmds_ggplot, 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 5 
)