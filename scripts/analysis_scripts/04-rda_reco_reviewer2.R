################################################################################
# Accompanying code for the following research project: 
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
#
#
# Purpose of this R script: to reanalyze the species composition data
# based on reviewer concerns from the Restoration Ecology submission

# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html


# re-run the original analysis for species composition, but split 
# it into two separate constrained correspondence analyses:

# (1) constrain restoration stage, but conditioning on sites + sampling season

# (2) remove newly-established sites (to acknowledge the artefact), 
# and constrain sampling season, while conditioning on site + restoration stage

# libraries --------------------------------------------------------------------
library(here)       # for creating relative file-paths
library(dplyr)      # for manipulating data
library(ggplot2)    # for visualizing data
library(vegan)      # for doing community ecology analyses
library(ggrepel)    # for repelling text in ggplot2 figures
library(rdacca.hp)  # for calculating variance partitioning 
library(stringr)    # for manipulating string data

# import -----------------------------------------------------------------------

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

# clean data -------------------------------------------------------------------

sb <- rbind(fall_sb, spr_sb)

comm_matrix <- sb %>%
  group_by(season, treatment, site_code) %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = total_abund) %>%
  ungroup() %>%
  mutate(across(.cols = AMRE:LOPE, ~ tidyr::replace_na(.x, 0))) %>%
  mutate(
    season = factor(season, levels = c("Spring", "Fall")),
    treatment = factor(treatment, levels = c("TIL", "RES", "MOW"))
    )

# run RDA (with newly-established sites) ---------------------------------------

## prep ------------------------------------------------------------------------

spp_comp_w_til <- comm_matrix %>%
  dplyr::select(-c("season", "site_name", "site_code", "treatment", "plot")) %>%
  vegan::decostand("hellinger")

env_w_til <- comm_matrix %>%
  dplyr::select(season, site_code, treatment)

## run the model ---------------------------------------------------------------

rda_w_til <- rda(
  spp_comp_w_til ~ treatment*season + site_code, 
  data = comm_matrix
  )

## tests of significance -------------------------------------------------------

anova.cca(rda_w_til, permutations = how(nperm = 999), by = "margin")
anova.cca(rda_w_til, permutations = how(nperm = 999), by = "axis")

## variance partitioning -------------------------------------------------------

vp_w_til <- rdacca.hp(
  dv = spp_comp_w_til, 
  iv = env_w_til, 
  method = "RDA",
  type = "adjR2"
)

# run RDA (without newly-established sites) ------------------------------------

## prep ------------------------------------------------------------------------

spp_comp_no_til <- comm_matrix %>%
  filter(!treatment %in% c("TIL")) %>%
  dplyr::select(-c("season", "site_name", "site_code", "treatment", "plot")) %>%
  vegan::decostand("hellinger")

env_no_til <- comm_matrix %>%
  filter(!treatment %in% c("TIL")) %>%
  dplyr::select(season, site_code, treatment)

## run the model ---------------------------------------------------------------

rda_no_til <- rda(
  spp_comp_no_til ~ treatment*season + site_code, 
  data = env_no_til
  )

## tests of significance -------------------------------------------------------

anova.cca(rda_no_til, permutations = how(nperm = 999), by = "margin")
anova.cca(rda_no_til, permutations = how(nperm = 999), by = "axis")

## variance partitioning -------------------------------------------------------

vp_no_til <- rdacca.hp(
  dv = spp_comp_no_til, 
  iv = env_no_til, 
  method = "RDA",
  type = "adjR2"
)

# visualize data (includes newly-established sites) ----------------------------

## manually extract scores for the first two RDA axes --------------------------

rda_summ_w_til <- summary(rda_w_til)

sp_scores <- as.data.frame(rda_summ_w_til$species[, c("RDA1", "RDA2")])
st_scores <- as.data.frame(rda_summ_w_til$sites[, c("RDA1", "RDA2")])
yz_scores <- as.data.frame(rda_summ_w_til$biplot[, c("RDA1", "RDA2")]) 

yz_scores$vars <- rownames(yz_scores) 
yz_tidy <- filter(yz_scores, !(str_detect(vars, "site_code")))
rownames(yz_tidy) <- c("Restored", "Mown", "Fall", "Res X Fall", "Mown X Fall")

## get variation explained for each RDA axis -----------------------------------

rda1_prop_explained <- rda_summ_w_til$cont$importance["Proportion Explained","RDA1"]
rda2_prop_explained <- rda_summ_w_til$cont$importance["Proportion Explained","RDA2"]

rda1_prop_explained <- round(rda1_prop_explained*100, digits = 0)
rda2_prop_explained <- round(rda2_prop_explained*100, digits = 0)

## plot the scores using ggplot2 syntax ----------------------------------------
(ggplot_rda <- 
   
   ggplot() + 
   
   # species scores
   geom_text_repel(
     aes(x = RDA1, y  = RDA2), 
     label = row.names(sp_scores), 
     max.overlaps = 100, 
     alpha = 0.3, 
     data = sp_scores) +
   
   geom_point(
     aes(x = RDA1, y  = RDA2), 
     alpha = 0.5, 
     data = sp_scores
     ) + 
   
   # environmental scores 
   geom_text_repel(
     aes(x = RDA1, y = RDA2), 
     label = rownames(yz_tidy), 
     data = yz_tidy
   ) + 
   
   geom_segment(
     aes(x = 0, xend = RDA1, y = 0, yend = RDA2),
     colour = "black",
     alpha = 0.5,
     data = yz_tidy
   ) + 
   
   labs(
     x = paste(
       "RDA1", 
       "(", 
       rda1_prop_explained, 
       "% variation explained)", 
       sep = ""),
     
     y = paste(
       "RDA2", 
       "(", 
       rda2_prop_explained, 
       "% variation explained)",
       sep = "")
   ) + 
   
   geom_hline(yintercept = 0, linetype = 3, linewidth = 0.1) + 
   geom_vline(xintercept = 0, linetype = 3, linewidth  = 0.1) + 
   
   xlim(-0.8, 0.8) + 
   ylim(-0.7, 0.7) +
   

   theme(
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(),
     panel.background = element_blank(), 
     axis.line = element_line(colour = "black")
     )
)

# visualize data (excludes newly-established sites) ----------------------------

## manually extract scores for the first two RDA axes --------------------------

rda_summ_no_til <- summary(rda_no_til)

sp_scores2 <- as.data.frame(rda_summ_no_til$species[, c("RDA1", "RDA2")])
st_scores2 <- as.data.frame(rda_summ_no_til$sites[, c("RDA1", "RDA2")])
yz_scores2 <- as.data.frame(rda_summ_no_til$biplot[, c("RDA1", "RDA2")]) 

yz_scores2$vars <- rownames(yz_scores2) 
yz_tidy2 <- filter(yz_scores2, !(str_detect(vars, "site_code")))
rownames(yz_tidy2) <- c("Mown", "Fall", "Mown X Fall")

## get variation explained for each RDA axis -----------------------------------

rda1_prop_explained2 <- rda_summ_no_til$cont$importance["Proportion Explained","RDA1"]
rda2_prop_explained2 <- rda_summ_no_til$cont$importance["Proportion Explained","RDA2"]

rda1_prop_explained2 <- round(rda1_prop_explained2*100, digits = 0)
rda2_prop_explained2 <- round(rda2_prop_explained2*100, digits = 0)

## plot the scores using ggplot2 syntax ----------------------------------------
(ggplot_rda_no_til <- 
   
   ggplot() + 
   
   # species scores
   geom_text_repel(
     aes(x = RDA1, y  = RDA2), 
     label = row.names(sp_scores2), 
     max.overlaps = 100, 
     alpha = 0.3, 
     data = sp_scores2) +
   
   geom_point(
     aes(x = RDA1, y  = RDA2), 
     alpha = 0.5, 
     data = sp_scores2
   ) + 
   
   # environmental scores 
   geom_text_repel(
     aes(x = RDA1, y = RDA2), 
     label = rownames(yz_tidy2), 
     data = yz_tidy2
   ) + 
   
   geom_segment(
     aes(x = 0, xend = RDA1, y = 0, yend = RDA2),
     colour = "black",
     alpha = 0.5,
     data = yz_tidy2
   ) + 
   
   labs(
     x = paste(
       "RDA1", 
       "(", 
       rda1_prop_explained2, 
       "% variation explained)", 
       sep = ""),
     
     y = paste(
       "RDA2", 
       "(", 
       rda2_prop_explained2, 
       "% variation explained)",
       sep = "")
   ) + 
   
   geom_hline(yintercept = 0, linetype = 3, linewidth = 0.1) + 
   geom_vline(xintercept = 0, linetype = 3, linewidth  = 0.1) + 
   
   xlim(-0.7, 0.7) + 
   ylim(-0.9, 0.9) + 
   
   theme(
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(),
     panel.background = element_blank(), 
     axis.line = element_line(colour = "black")
   )
)

# summary statistics for specific species --------------------------------------

# get the most abundant species according to species scores in the RDA biplot

## newly-established stage -----------------------------------------------------

chal_til <- sb %>% 
  filter(treatment == "TIL", spp_code == "CHAL") %>%
  pull(total_abund) %>%
  sum()

melu_til <- sb %>%
  filter(treatment == "TIL", spp_code == "MELU") %>%
  pull(total_abund) %>%
  sum()

oxst_til <- sb %>%
  filter(treatment == "TIL", spp_code == "OXST") %>%
  pull(total_abund) %>%
  sum()

## mown restored stage ---------------------------------------------------------

oebi_mow <- sb %>% 
  filter(treatment == "MOW", spp_code == "OEBI") %>%
  pull(total_abund) %>%
  sum()

coca_mow <- sb %>% 
  filter(treatment == "MOW", spp_code == "COCA") %>%
  pull(total_abund) %>%
  sum()

mofi_mow <- sb %>% 
  filter(treatment == "MOW", spp_code == "MOFI") %>%
  pull(total_abund) %>%
  sum() 

## unmown restored stage -------------------------------------------------------

sosp_res <- sb %>%
  filter(treatment == "RES", spp_code == "SOSP") %>%
  pull(total_abund) %>%
  sum()

poar_res <- sb %>%
  filter(treatment == "RES", spp_code == "POAR") %>%
  pull(total_abund) %>%
  sum()

eran_res <- sb %>%
  filter(treatment == "RES", spp_code == "ERAN") %>%
  pull(total_abund) %>%
  sum()

## spring season ---------------------------------------------------------------

melu_spr <- sb %>%
  filter(season == "Spring", spp_code == "MELU") %>%
  pull(total_abund) %>%
  sum()

trpr_spr <- sb %>%
  filter(season == "Spring", spp_code == "TRPR") %>%
  pull(total_abund) %>%
  sum()

trre_spr <- sb %>%
  filter(season == "Spring", spp_code == "TRRE") %>%
  pull(total_abund) %>%
  sum()

# after removing newly-established sites to account for sampling season

hegi_spr <- sb %>%
  filter(season == "Spring", spp_code == "HEGI", treatment != "TIL") %>%
  pull(total_abund) %>%
  sum()

oxst_spr <- sb %>%
  filter(season == "Spring", spp_code == "OXST", treatment != "TIL") %>%
  pull(total_abund) %>%
  sum()

pyte_spr <- sb %>%
  filter(season == "Spring", spp_code == "PYTE", treatment != "TIL") %>%
  pull(total_abund) %>%
  sum()

## fall season -----------------------------------------------------------------

chal_fal <- sb %>%
  filter(season == "Fall", spp_code == "CHAL") %>%
  pull(total_abund) %>%
  sum()

coca_fal <- sb %>%
  filter(season == "Fall", spp_code == "COCA") %>%
  pull(total_abund) %>%
  sum()

oebi_fal <- sb %>%
  filter(season == "Fall", spp_code == "OEBI") %>%
  pull(total_abund) %>%
  sum()

# after removing newly-established sites to account for sampling season

ruhi_fal <- sb %>%
  filter(season == "Fall", spp_code == "RUHI", treatment != "TIL") %>%
  pull(total_abund) %>%
  sum()

coca_fal <- sb %>%
  filter(season == "Fall", spp_code == "COCA", treatment != "TIL") %>%
  pull(total_abund) %>%
  sum()

mofi_fal <- sb %>%
  filter(season == "Fall", spp_code == "MOFI", treatment != "TIL") %>%
  pull(total_abund) %>%
  sum()

# save to disk -----------------------------------------------------------------

ggsave(
  filename = here("output", "results", "figure-5-rda.png"), 
  plot = ggplot_rda, 
  device = "png", 
  height = 5, 
  width = 5, 
  units = "in"
)

ggsave(
  filename = here("output", "results", "figure-5-rda.pdf"), 
  plot = ggplot_rda, 
  device = "pdf", 
  height = 5, 
  width = 5, 
  units = "in"
)

ggsave(
  filename = here("output", "data_appendix_output", "figure-s4-rda_no_til.png"), 
  plot = ggplot_rda_no_til, 
  device = "png", 
  height = 5, 
  width = 5, 
  units = "in"
)

ggsave(
  filename = here("output", "data_appendix_output", "figure-s4-rda_no_til.pdf"), 
  plot = ggplot_rda_no_til, 
  device = "pdf", 
  height = 5, 
  width = 5, 
  units = "in"
)

