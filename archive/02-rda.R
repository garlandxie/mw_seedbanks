# libraries ----
library(here)
library(dplyr)
library(ggplot2)
library(vegan)
library(ggrepel)

# import ----

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

# check packaging ----
glimpse(spr_sb)
glimpse(fall_sb)

# clean data ----

sb <- rbind(fall_sb, spr_sb)

# clean data: community data matrix -----

comm_matrix <- sb %>%
  group_by(season, treatment, site_code) %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = total_abund) %>%
  ungroup() %>%
  mutate(across(.cols = everything(), ~ tidyr::replace_na(.x, 0)))

# redundancy analysis -----

## prep ----

Y <- comm_matrix %>%
  dplyr::select(-c("season", "site_name", "site_code", "treatment", "plot")) %>%
  vegan::decostand("hellinger")

## run RDA ----
comm_comp_rda <- rda(Y ~ season + treatment + site_name, data = comm_matrix)

## grab eigenvalues ----
rda_summ <- summary(comm_comp_rda)

## get r squared ----

# unadjusted R^2 retrieved from rda object
R2 <- RsquareAdj(comm_comp_rda)$r.squared

# adjusted R^2 retrieved from rda object
R2adj <- RsquareAdj(comm_comp_rda)$adj.r.squared

## anova ----

anova(comm_comp_rda, permutations = how(nperm = 999))
anova(comm_comp_rda, by = "axis", permutations = how(nperm = 999))

# plot ----

# manually extract scores for the first two RDA axes
sp_scores <- as.data.frame(rda_summ$species[, c("RDA1", "RDA2")])
st_scores <- as.data.frame(rda_summ$sites[, c("RDA1", "RDA2")])
yz_scores <- as.data.frame(rda_summ$biplot[, c("RDA1", "RDA2")]) # remove sites for now

# show only management regimes (instead of site names) for readability
yz_scores$vars <- rownames(yz_scores) 
yz_scores <- filter(
  yz_scores, 
  vars %in% c(
    "seasonSpring", 
    "treatmentRES",
    "treatmentTIL"
    )
  ) 

# plot the scores using ggplot2 syntax
(ggplot_rda <- 
  
  ggplot() + 

  # species scores
  geom_text_repel(
    aes(x = RDA1, y  = RDA2), 
    label = row.names(sp_scores), 
    max.overlaps = 50, 
    alpha = 0.3, 
    data = sp_scores) +
  geom_point(
    aes(x = RDA1, y  = RDA2), 
    alpha = 0.5, 
    data = sp_scores) +
  
  # environmental scores 
  geom_text_repel(
    aes(x = RDA1, y = RDA2), 
    label = row.names(yz_scores), 
    data = yz_scores
  ) + 
  geom_hline(yintercept = 0, linetype = 3, size = 0.1) + 
  geom_vline(xintercept = 0, linetype = 3, size = 0.1) + 
  xlim(-1, 1) + 
  ylim(-1, 1) + 
  theme_bw()
)

# save to disk ------

ggsave(
  filename = here("output", "results", "rda_community_composition.svg"),
  plot = ggplot_rda, 
  device = "svg", 
  units = "in", 
  height = 5, 
  width = 5
)

