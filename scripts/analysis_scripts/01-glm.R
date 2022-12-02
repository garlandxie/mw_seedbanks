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
# Purpose of this R script: to conduct regression models to determine 
# relationship between seed bank community metrics, management regimes, and
# seasonal variation
#
# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html


# libraries --------------------------------------------------------------------
library(here)          # for creating relative file-paths
library(dplyr)         # for manipulating data 
library(ggplot2)       # for visualizing data
library(performance)   # for calculating R-square values
library(MASS)      
library(patchwork)     # for creating multi-panel plots
library(lme4)          # for running linear mixed effect models
library(emmeans)       # for calculating pairwise comparisons
library(tidyr)         # for manipulating data
library(DHARMa)        # for running diagnostic tests for LMM's
library(tibble)        # for adding rows to a data frame

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
glimpse(spr_sb)
glimpse(fall_sb)

# clean data -------------------------------------------------------------------

# combine both fall and spring seed bank data
sb <- rbind(fall_sb, spr_sb)

# get community-level abundance and species richness
sb_comm <- sb %>%
  group_by(season, site_name, site_code, treatment, plot) %>%
  summarize(
    abund = sum(total_abund, na.rm = TRUE),
    species_richness = dplyr::n_distinct(spp_code)
    ) %>%
  ungroup() %>%
  mutate(
    season = factor(season, levels = c("Spring", "Fall")), 
    treatment = factor(treatment, levels = c("RES", "MOW", "TIL"))
  )

# check for trays with no seed bank information
sb_comm %>%
  group_by(site_code) %>%
  summarize(n = n()) %>%
  arrange(n)

# add zeros to plots with no seed bank info (representing zero germination)
sb_comm_tidy <- sb_comm %>%
  
  # Spring, T2_1, plot 13
  tibble::add_row(
    season = "Spring", 
    site_name = "VICP", 
    site_code = "T2_1",
    treatment = "TIL",
    plot = 13, 
    abund = 0,
    species_richness = 0, 
    ) %>%
  
  # Spring, T2_1, plot 15
  tibble::add_row(
    season = "Spring", 
    site_name = "VICP", 
    site_code = "T2_1",
    treatment = "TIL",
    plot = 15, 
    abund = 0,
    species_richness = 0, 
  ) %>%

  # Spring, T2_3, plot 21
  tibble::add_row(
    season = "Spring", 
    site_name = "KENN", 
    site_code = "T2_3",
    treatment = "TIL",
    plot = 21, 
    abund = 0,
    species_richness = 0, 
  ) %>%
  
  # Spring, T2_3, plot 25
  tibble::add_row(
    season = "Spring", 
    site_name = "KENN", 
    site_code = "T2_3",
    treatment = "TIL",
    plot = 25, 
    abund = 0,
    species_richness = 0, 
  ) 

# exploring data (before regression) -------------------------------------------

## step 1: check for outliers  -------------------------------------------------

# seed density
sb_comm %>%
  ggplot(aes(y = row.names(sb_comm), x = abund)) + 
  geom_point() + 
  labs(
    x = "Seedling Emergent Abundance (Plot-Level)",
    y = "Order of the data") +
  theme(axis.text.y = element_blank()) 

# seed richness
sb_comm %>%
  ggplot(aes(y = row.names(sb_comm), x = species_richness)) + 
  geom_point() + 
  labs(
    x = "Seed Richness (Plot-Level)",
    y = "Order of the data") +
  theme(axis.text.y = element_blank()) 

## step 2: check for homogeneity of variance -----------------------------------

# multi-panel conditional box-plots
sb_comm %>%
  ggplot(aes(x = site_code, y = abund, col = treatment)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~season, ncol = 2, nrow = 3) + 
  labs(x = "Site", y = "Seedling Emergent Abundance") + 
  coord_flip() + 
  theme_bw() +
  theme(legend.position = "none")

sb_comm %>%
  ggplot(aes(x = site_code, y = species_richness, col = treatment)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~season, ncol = 2, nrow = 3) + 
  scale_color_discrete(
    name = "Management Regime", 
    labels = c("Maintenance-Mow", "Undisturbed", "Heavily-Tilled")
  ) + 
  ylim(0, 20) + 
  labs(
    x = NULL, y = "Species Richness") + 
  coord_flip() + 
  theme_bw() + 
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
    )

## step 3: check for normally-distributed data ---------------------------------

# freedman-diaonis rule

#' bins_fd
#' @description returns the number of bins according to the Freedman-Diaconis rule
#' @param vec numeric vector
#' @return number of bins
bins_fd <- function(vec) {
  diff(range(vec)) / (2 * IQR(vec) / length(vec)^(1/3))
}

sb_comm %>%
  ggplot(aes(x = species_richness)) + 
  geom_histogram(bins = bins_fd(sb_comm$species_richness)) + 
  theme_bw()

sb_comm %>%
  ggplot(aes(x = abund)) + 
  geom_histogram(bins = bins_fd((sb_comm$abund))) + 
  theme_bw()

# run regression models (with outliers) ----------------------------------------

## abundance -------------------------------------------------------------------

glmer_abund_poisson <- glmer(
  formula = abund ~ treatment + season + (1|site_code), 
  family = poisson(link = "log"), 
  data = sb_comm
  )

### check overdispersion -------------------------------------------------------
sim_glmer_abund <- 
  simulateResiduals(
    fittedModel = glmer_abund_poisson,
    plot = F
    )

testDispersion(
  sim_glmer_abund, 
  type = "PearsonChisq",
  alternative = 'greater'
  )

### account for overdispersion -------------------------------------------------

# overdispersion likely due to :
# (1) variability in detection efficiency, and 
# (2) environmental stochasticity
# see scenario 4 in Linden and Mäntyniemi (2011)
# https://esajournals.onlinelibrary.wiley.com/doi/10.1890/10-1831.1
glmer_abund_nb <- glmer.nb(
  formula = abund ~ treatment + season + (1|site_code), 
  data = sb_comm
)

### model summary --------------------------------------------------------------
summary(glmer_abund_nb)

### check diagnostics ----------------------------------------------------------
plot(simulateResiduals(fittedModel = glmer_abund_nb, plot = F))

### r-squared ------------------------------------------------------------------
r2(glmer_abund_nb)

### pairwise comparison --------------------------------------------------------
abund_emm_trt <- emmeans(
  glmer_abund_nb, 
  "treatment", 
  lmer.df = "satterthwaite"
  )

abund_emm_sn <- emmeans(
  glmer_abund_nb, 
  "season", 
  lmer.df = "satterhwaite"
)

# need to get mow vs til comparisons here

abund_sn_tuk <- glht(glmer_abund_nb, linfct = mcp(season = "Tukey"))
abund_sn_cld <- cld(abund_sn_tuk)

abund_trt_tuk <- glht(glmer_abund_nb, linfct = mcp(treatment = "Tukey"))
abund_trt_cld <- cld(abund_trt_tuk)


## species richness ------------------------------------------------------------

glmer_richness <- glmer(
  formula = species_richness ~ treatment + season + (1|site_code), 
  family = poisson(link = "log"),
  data = sb_comm
)

### check overdispersion -------------------------------------------------------
check_overdispersion(glmer_richness)

### check model diagnostics ----------------------------------------------------

plot(simulateResiduals(fittedModel = glmer_richness, plot = F))

### model summary --------------------------------------------------------------
summary(glmer_richness)

### r squared ------------------------------------------------------------------
r2(glmer_richness)

### pairwise comparisons -------------------------------------------------------
richness_emm_trt <- emmeans(
  glmer_richness, 
  "treatment", 
  lmer.df = "satterthwaite"
)

richness_emm_sn <- emmeans(
  glmer_richness, 
  "season", 
  lmer.df = "satterhwaite"
)

# need to get mow vs til comparisons here
richness_pairs_trt <- richness_emm_trt %>%
  pairs() %>%
  as.data.frame()

richness_pairs_sn <- richness_emm_trt %>%
  pairs() %>%
  as.data.frame()

## save to disk ----------------------------------------------------------------
ggsave(
  filename = here("output", "results", "multi_plot.svg"), 
  plot = multi_plot, 
  device = "svg", 
  units = "in", 
  height = 4, 
  width = 9
)

write.csv(
  abund_pairs_trt, 
  file = here("output", "results", "abund_pairs_trt.csv")
)

write.csv(
  abund_pairs_sn, 
  file = here("output", "results", "abund_pairs_sn.csv")
)

write.csv(
  richness_pairs_trt, 
  file = here("output", "results", "richness_pairs_trt.csv")
)

write.csv(
  richness_pairs_sn, 
  file = here("output", "results", "richness_pairs_sn.csv")
)




  