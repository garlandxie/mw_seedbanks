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
library(ggsignif)

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

# run regression models --------------------------------------------------------

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

## species richness ------------------------------------------------------------

glmer_richness <- glmer(
  formula = species_richness ~ treatment + season + (1|site_code), 
  family = poisson(link = "log"),
  data = sb_comm
)

### check overdispersion -------------------------------------------------------
sim_glmer_richness <- 
  simulateResiduals(
    fittedModel = glmer_richness,
    plot = F
  )

testDispersion(
  sim_glmer_richness, 
  type = "PearsonChisq",
  alternative = 'greater'
)

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

# visualize data ---------------------------------------------------------------

# adjust x-axis labels (i.e., management regimes)
# for readability
sb_data_viz <- sb_comm %>%
  
  mutate(
    treatment = as.character(treatment), 
    treatment = case_when(
      treatment == "MOW" ~ "Maintenance-Mow", 
      treatment == "RES" ~ "Undisturbed",
      treatment == "TIL" ~ "Tilling"
    ), 
    treatment = factor(
      treatment, levels = c("Undisturbed", "Maintenance-Mow", "Tilling")) 
  ) 

## color-blind friendly palette ------------------------------------------------

cbPalette <- c("#009E73", "#E69F00")

## pairwise comparisons --------------------------------------------------------

# calculate estimated marginal means
abund_emm_trt <- emmeans(
  glmer_abund_nb, 
  "treatment", 
  lmer.df = "satterthwaite"
)

# get summary of comparisons, coefficients, and p-values
pairs_abund_trt <- as.data.frame(pairs(abund_emm_trt))

# obtain p-value for comparison between 
# undisturbed and tilling
pairs_til_res <- pairs_abund_trt %>%
  filter(contrast == "RES - TIL") %>%
  pull(p.value) %>%
  signif(digits = 1)

# obtain p-value for comparison between 
# maintenance-mowing and undisturbed
pairs_mow_res <- pairs_abund_trt %>%
  filter(contrast == "RES - MOW") %>%
  pull(p.value) %>%
  signif(digits = 1)

# obtain p-value for comparison between
# maintenance-mowing and tilling
pairs_mow_til <- pairs_abund_trt %>%
  filter(contrast == "MOW - TIL") %>%
  
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ 0.001)
  ) %>%
  pull(p.value) 

(sb_abund_plot <- sb_data_viz %>%
    ggplot(aes(x = treatment, y = abund, fill = season)) +
    geom_boxplot() + 
    
    # pairwise significance between tilling and undisturbed
    geom_signif(
      y_position = 300, 
      xmin = 1, 
      xmax = 2, 
      annotation = paste("p = ", as.character(pairs_til_res)),
      alpha = 0.5
    ) + 
    
    # pairwise significance between mowing and undisturbed  
    geom_signif(
      y_position = 400, 
      xmin = 2, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_mow_til)), 
      alpha = 0.5
    ) +   
    
    # pairwise significance between mowing and tilling
    geom_signif(
      y_position = 500, 
      xmin = 1, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_mow_res)),
      alpha = 0.5
    ) + 
    
    ylim(0, 800) + 
    
    labs(
      x = "Management Regime", 
      y = "Seedling Emergent Abundance"
    ) + 
    
    scale_fill_discrete(name = "Season") + 
    
    theme_bw() + 
    theme(
      axis.title.x = element_text(
        margin = margin(t = 10, r = 0, b = 0, l = 0)
      ),
      axis.title.y = element_text(
        margin = margin(t = 0, r = 10, b = 0, l = 0)
      ),
      text = element_text(size = 15)
    )
)

## pairwise comparisons: richness ----------------------------------------------

rich_emm_trt <- emmeans(
  glmer_richness, 
  "treatment", 
  lmer.df = "satterthwaite"
)

# get summary of comparisons, coefficients, and p-values
pairs_rich_trt <- as.data.frame(pairs(rich_emm_trt))

# obtain p-value for comparison between 
# undisturbed and tilling
pairs_rich_til_res <- pairs_rich_trt %>%
  filter(contrast == "RES - TIL") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ 0.001)
  ) %>%
  pull(p.value) 

# obtain p-value for comparison between 
# maintenance-mowing and undisturbed
pairs_rich_mow_res <- pairs_rich_trt %>%
  filter(contrast == "RES - MOW") %>%
  pull(p.value) %>%
  signif(digits = 1)

# obtain p-value for comparison between
# maintenance-mowing and tilling
pairs_rich_mow_til <- pairs_rich_trt %>%
  filter(contrast == "MOW - TIL") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ 0.001)
  ) %>%
  pull(p.value) 

(sb_rich_plot <- sb_data_viz %>%
   ggplot(aes(x = treatment, y = species_richness, fill = season)) +
   geom_boxplot() + 
   
   # pairwise significance between tilling and undisturbed
   geom_signif(
     y_position = 20, 
     xmin = 1, 
     xmax = 2, 
     annotation = paste("p = <", as.character(pairs_til_res)),
     alpha = 0.5
   ) + 
   
   # pairwise significance between mowing and undisturbed  
   geom_signif(
     y_position = 22, 
     xmin = 2, 
     xmax = 3, 
     annotation = paste("p = ", as.character(pairs_mow_res)), 
     alpha = 0.5
   ) +   
   
   # pairwise significance between mowing and tilling  
   geom_signif(
     y_position = 25, 
     xmin = 1, 
     xmax = 3, 
     annotation =  paste("p = <", as.character(pairs_mow_til)),
     alpha = 0.5
   ) + 
   
   labs(
     x = "Management Regime", 
     y = "Seedling Emergent Richness"
   ) + 
   
   ylim(0, 30)  + 
   
   scale_fill_discrete(name = "Season") + 
   
   theme_bw() + 
   theme(
     axis.title.x = element_text(
       margin = margin(t = 10, r = 0, b = 0, l = 0)
     ),
     axis.title.y = element_text(
       margin = margin(t = 0, r = 10, b = 0, l = 0)
     ),
     text = element_text(size = 15)
   )
)

# save to disk -----------------------------------------------------------------

## figure 2: seedling eabundance -----------------------------------------------
ggsave(
  filename = here("output", "results", "figure-2_abund.png"), 
  plot = sb_abund_plot, 
  device = "png", 
  units = "in", 
  height = 4, 
  width = 6
)

## figure 3: seedling richness -------------------------------------------------
ggsave(
  filename = here("output", "results", "figure-3_richness.png"), 
  plot = sb_rich_plot, 
  device = "png", 
  units = "in", 
  height = 4, 
  width = 6
)




  