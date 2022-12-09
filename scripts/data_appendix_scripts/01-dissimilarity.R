# libraries --------------------------------------------------------------------
library(vegan)     # for doing community ecology analyses
library(ggplot2)   # for visualizing data
library(dplyr)     # for manipulating data
library(here)      # for creating relative file-paths
library(patchwork) # for creating multi-panel plots
library(tidyr)     # for creating wide tables
library(ape)

# import -----------------------------------------------------------------------

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

sb_taxon <- read.csv(
  here("data", "input_data", 
       "seed_bank_taxonomy.csv")
)

# check packaging --------------------------------------------------------------
str(spr_sb)
str(fall_sb)

# clean data -------------------------------------------------------------------

# merge fall and spring season 
sb <- rbind(fall_sb, spr_sb)

## create spring community data matrix -----------------------------------------
spring_comm_matrix <- sb %>%
  filter(season == "Spring") %>%
  group_by(treatment, site_code, plot, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    id_cols = c("treatment", "site_code", "plot"),  
    names_from = spp_code, 
    values_from = abundance
  ) %>%
  mutate(across(.cols = AMSP:SONU, ~ tidyr::replace_na(.x, 0))) %>%
  as.data.frame()

# assign row labels for treatment and plot id
rownames(spring_comm_matrix) <- paste(
  spring_comm_matrix$treatment, 
  c(1:length(spring_comm_matrix$treatment))
)

## create fall community data matrix -------------------------------------------
fall_comm_matrix <- sb %>%
  filter(season == "Fall") %>%
  group_by(treatment, site_code, plot, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    id_cols = c("treatment", "site_code", "plot"),  
    names_from = spp_code, 
    values_from = abundance
  ) %>%
  mutate(across(.cols = AMRE:RARA, ~ tidyr::replace_na(.x, 0))) %>%
  as.data.frame()

# assign row labels for treatment and plot id
rownames(fall_comm_matrix) <- paste(
  fall_comm_matrix$treatment, 
  c(1:length(fall_comm_matrix$treatment))
)

# calculate dissimilarities  ---------------------------------------------------

## spring ----------------------------------------------------------------------

# use Bray-Curtis dissimilarity to account for common and rare species
# common and rare species are determined by number of seedling emergents
bray_spr <- spring_comm_matrix %>%
  dplyr::select(-c("treatment", "site_code", "plot"))

bray_spr_dist <- vegan::vegdist(bray_spr, index = "bray")

# do principal coordinate analysis
pcoa_spr <- ape::pcoa(D = bray_spr_dist, correction = "cailliez")

# prepare for data visualization
pcoa_spr_df <- pcoa_spr$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "treatment") %>%
  dplyr::select(treatment, Axis.1, Axis.2) %>%
  mutate(treatment = stringr::str_sub(treatment, start = 1L, end = 3L))

## fall ------------------------------------------------------------------------

# use Bray-Curtis dissimilarity to account for common and rare species
# common and rare species are determined by abundance (# of seedling emergents)
bray_fall <- fall_comm_matrix %>%
  dplyr::select(-c("treatment", "site_code", "plot"))

bray_fall_dist <- vegdist(bray_fall, index = "bray")

pcoa_fall <- ape::pcoa(D = bray_fall_dist, correction = "cailliez")

# prepare for data visualization
pcoa_fall_df <- pcoa_fall$vectors %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "treatment") %>%
  dplyr::select(treatment, Axis.1, Axis.2) %>%
  mutate(treatment = stringr::str_sub(treatment, start = 1L, end = 3L))

# visualize the data -----------------------------------------------------------

## spring ----------------------------------------------------------------------

# variation explained in the first axis
spr_pcoa1_rel_eig <- round(pcoa_spr$values$Rel_corr_eig[1]*100, digits = 0)
pcoa1_var_explained <- paste(
  "PCoA1", " (", 
  spr_pcoa1_rel_eig, 
  "% variation explained)", 
  sep = ""
  )

# variation explained in the second axis
spr_pcoa2_rel_eig  <- round(pcoa_spr$values$Rel_corr_eig[2]*100, digits = 0)
pcoa2_var_explained <- paste(
  "PCoA2", " (", 
  spr_pcoa2_rel_eig, 
  "% variation explained)", 
  sep = ""
)

(pcoa_spr_plot <- pcoa_spr_df %>%
    ggplot(aes(x = Axis.1, y = Axis.2, shape = treatment)) + 
    geom_point() + 
    xlim(-0.8, 0.8) + 
    ylim(-0.7, 0.7) + 
    scale_shape_manual(
      name = "Management Regime", 
      labels = c("Maintenance-Mow", "Undisturbed", "Tilled"),
      values = c(0, 1, 2)
    ) + 
    labs(
      title = "a) Spring season",
      x = pcoa1_var_explained,
      y = pcoa2_var_explained
    ) + 
    theme_bw() + 
    theme(legend.position = "none")
)

## fall ------------------------------------------------------------------------

# variation explained in the first axis
fall_pcoa1_rel_eig <- round(pcoa_fall$values$Rel_corr_eig[1]*100, digits = 0)
fall_pcoa1_var_explained <- paste(
  "PCoA1", " (", 
  fall_pcoa1_rel_eig, 
  "% variation explained)", 
  sep = ""
)

# variation explained in the second axis
fall_pcoa2_rel_eig <- round(pcoa_fall$values$Rel_corr_eig[2]*100, digits = 0)
fall_pcoa2_var_explained <- paste(
  "PCoA2", " (", 
  fall_pcoa2_rel_eig, 
  "% variation explained)", 
  sep = ""
)

(pcoa_fall_plot <- pcoa_fall_df %>%
    ggplot(aes(x = Axis.1, y = Axis.2, shape = treatment)) + 
    geom_point() + 
    xlim(-0.8, 0.8) + 
    ylim(-0.7, 0.7) + 
    scale_shape_manual(
      name = "Management Regime", 
      labels = c("Maintenance-Mow", "Undisturbed", "Tilled"),
      values = c(0, 1, 2)
    ) + 
    labs(
      title = "b) Fall season",
      x = fall_pcoa1_var_explained,
      y = fall_pcoa2_var_explained
    ) + 
    theme_bw() 
)

## multi-panel plot ------------------------------------------------------------

(pcoa_plot <- pcoa_spr_plot + pcoa_fall_plot)

# save disk --------------------------------------------------------------------

ggsave(
  filename = here("output", "data_appendix_output", "figure-s1_pcoa.png"), 
  plot = pcoa_plot, 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 10
)


