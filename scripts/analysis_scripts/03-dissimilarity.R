# libraries --------------------------------------------------------------------
library(vegan)     # for doing community ecology analyses
library(ggplot2)   # for visualizing data
library(dplyr)     # for manipulating data
library(here)      # for creating relative file-paths
library(patchwork) # for creating multi-panel plots
library(tidyr)     # for creating wide tables

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
pcoa_spr <- vegan::cmdscale(
  bray_spr_dist, 
  k = (nrow(bray_spr)-1),
  add = TRUE, # Cailliez correction (to correct negative eigenvalues)
  eig = TRUE)

pcoa_spr_df <- scores(pcoa_spr) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "treatment") %>%
  dplyr::select(treatment, Dim1, Dim2) %>%
  mutate(treatment = stringr::str_sub(treatment, start = 1L, end = 3L))

## fall ------------------------------------------------------------------------

# use Bray-Curtis dissimilarity to account for common and rare species
# common and rare species are determined by abundance (# of seedling emergents)
bray_fall <- fall_comm_matrix %>%
  dplyr::select(-c("treatment", "site_code", "plot"))

bray_fall_dist <- vegdist(bray_fall, index = "bray")

# do principal coordinate analysis
pcoa_fall <- cmdscale(
  bray_fall_dist, 
  k = (nrow(bray_fall)-1),
  add = TRUE, # Cailliez correction
  eig = TRUE)

pcoa_fall_df <- scores(pcoa_fall) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "treatment") %>%
  dplyr::select(treatment, Dim1, Dim2) %>%
  mutate(treatment = stringr::str_sub(treatment, start = 1L, end = 3L))

# visualize the data -----------------------------------------------------------

## spring ----------------------------------------------------------------------

(pcoa_spr_plot <- pcoa_spr_df %>%
    ggplot(aes(x = Dim1, y = Dim2, col = treatment)) + 
    geom_point() + 
    xlim(-0.8, 0.8) + 
    ylim(-0.7, 0.7) + 
    labs(
      title = "a) Spring season",
      x = "PCoA1", 
      y = "PCoA2"
    ) + 
    theme_bw() + 
    theme(legend.position = "none")
)

## fall ------------------------------------------------------------------------

(pcoa_fall_plot <- pcoa_fall_df %>%
   ggplot(aes(x = Dim1, y = Dim2, col = treatment)) + 
   geom_point() + 
   scale_color_discrete(
     name = "Management Regime", 
     labels = c("Maintenance-Mow", "Undisturbed", "Tilled")
   ) + 
   xlim(-0.8, 0.8) + 
   ylim(-0.7, 0.7) + 
   labs(
     title = "b) Fall season",
     x = "PCoA1", 
     y = NULL
   ) + 
   theme_bw()
)

## multi-panel plot ----

(pcoa_plot <- pcoa_spr_plot + pcoa_fall_plot)

# save disk ----

ggsave(
  filename = here("output", "results", "pcoa.svg"), 
  plot = pcoa_plot, 
  device = "svg", 
  units = "in", 
  height = 5, 
  width = 8
)


