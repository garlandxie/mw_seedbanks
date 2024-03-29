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
# Purpose of this R script: to reanalyze the seed bank data based on reviewer 
# concerns from the Restoration Ecology submission

# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html


# Summary of reviewer concerns -------------------------------------------------

# There is an artefact in the original analysis with the stated aim
# comparing spring and autumn sampling to test the efficiency 
# of the two sampling periods of the monitoring of restoration 
# success in the Meadoway

# In brief, the artefect is likely caused by a flush of emergence of what
# may have included many presumably unsown species, which could then 
# replenish the previously impoverished soil seed bank

# This would explain the very strong observed hike in the seed bank density
# between your spring and autumn sampling seen in the original analysis 

# The main workflow for this analysis is:

# 1) re-run the original analysis for seed bank density, but include an 
# interaction term between sampling period and restoration stage. Acknowledge 
# the artefact (stated above) 

# 2) run a similar analysis (as above), but instead remove the newly-established
# site. If the interpretations are different, then address this in the results
# and the discussion


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
library(lmerTest)

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
  ) %>%
  
  mutate(
    season = factor(season, levels = c("Spring", "Fall")), 
    treatment = factor(treatment, levels = c("TIL", "RES", "MOW"))
  )

# run regression models --------------------------------------------------------

## seed bank density (with newly-established sites) ----------------------------

# accounting for Type I error control by applying the logarithm of the 
# response variable (seed bank density)

# create two models: one with the full model (including the interaction term)
# and one with the reduced model (excluding the interaction term)
# purpose: to compare the R^2 values (see below)

# reduced model
lm_abund_til1 <- sb_comm_tidy %>%
  mutate(log_abund = log(abund+1)) %>%
  lmer(formula = log_abund ~ treatment + season + (1|site_code), 
       data = .
  )

# full model
lm_abund_til2 <- sb_comm_tidy %>%
  mutate(log_abund = log(abund+1)) %>%
  lmer(formula = log_abund ~ treatment*season + (1|site_code), 
       data = .
  )

### model summary --------------------------------------------------------------
summary(lm_abund_til2)

### check diagnostics ----------------------------------------------------------
plot(simulateResiduals(fittedModel = lm_abund_til2, plot = F))

### model fit ------------------------------------------------------------------
r2(lm_abund_til2)

###  pairwise comparison --------------------------------------------------------

# get pairwise comparisons between restoration stages 
# while controlling for season
pairs_trt_w_til <- lm_abund_til2 %>%
  emmeans("treatment", type = "response") %>% 
  pairs() %>%
  data.frame()

# back-transform the response variables for effect sizes of restoration age
pairs_trt_w_til <- pairs_trt_w_til %>%
  mutate(bt_transform = exp(estimate) - 1) %>%
  dplyr::select(contrast, estimate, bt_transform, SE, df, t.ratio, p.value)

## seed bank density (without newly-established sites) ------------------------

# accounting for Type I error control by applying the logarithm of the 
# response variable (seed bank density)

sb_comm_tidy_no_til <- sb_comm_tidy %>%
  mutate(treatment = as.character(treatment)) %>%
  filter(treatment %in% c("MOW", "RES")) %>%
  mutate(
    treatment = as.factor(treatment),
    season = factor(season, levels = c("Fall", "Spring"))
    )

lm_abund_no_til <- sb_comm_tidy_no_til %>%
  mutate(log_abund = log(abund+1)) %>%
  lmer(formula = log_abund ~ treatment*season + (1|site_code), 
       data = .
  )

### model summary --------------------------------------------------------------
summary(lm_abund_no_til)

### check diagnostics ----------------------------------------------------------
plot(simulateResiduals(fittedModel = lm_abund_no_til, plot = F))

###  pairwise comparison --------------------------------------------------------

# get pairwise comparisons between restoration stages 
# while controlling for season
pairs_trt_no_til <- lm_abund_no_til %>%
  emmeans("treatment", type = "response") %>% 
  pairs() %>%
  data.frame()

# back-transform the response variables for effect sizes of restoration age
pairs_trt_no_til <- pairs_trt_no_til %>%
  mutate(bt_transform = exp(estimate) - 1) %>%
  dplyr::select(contrast, estimate, bt_transform, SE, df, t.ratio, p.value)

### model fit ------------------------------------------------------------------
r2(lm_abund_no_til)

### pairwise comparison --------------------------------------------------------
abund_emm_trt <- emmeans(
  lm_abund_til2, 
  "treatment", 
)

abund_sn_trt <- emmeans(
  lm_abund_til2, 
  "season", 
)

# visualize data ---------------------------------------------------------------

# adjust x-axis labels (i.e., management regimes)
# for readability
sb_data_viz <- sb_comm %>%
  
  mutate(
    treatment = as.character(treatment), 
    treatment = case_when(
      treatment == "MOW" ~ "Mown", 
      treatment == "RES" ~ "Unmown",
      treatment == "TIL" ~ "Newly-established"
    ), 
    treatment = factor(
      treatment, levels = c("Newly-established", "Unmown", "Mown")) 
  ) 

## color-blind friendly palette ------------------------------------------------

cbPalette <- c("#009E73", "#E69F00")

## pairwise comparisons: abundance ---------------------------------------------

# get summary of comparisons, coefficients, and p-values
pairs_abund_trt <- as.data.frame(pairs(abund_emm_trt))

# obtain p-value for comparison between 
# undisturbed and tilling
pairs_abund_til_res <- pairs_abund_trt %>%
  filter(contrast == "TIL - RES") %>%
  pull(p.value) %>%
  signif(digits = 1)

# obtain p-value for comparison between 
# maintenance-mowing and undisturbed
pairs_abund_mow_res <- pairs_abund_trt %>%
  filter(contrast == "RES - MOW") %>%
  pull(p.value) %>%
  signif(digits = 1)

# obtain p-value for comparison between
# maintenance-mowing and tilling
pairs_abund_mow_til <- pairs_abund_trt %>%
  filter(contrast == "TIL - MOW") %>%
  pull(p.value) %>%
  signif(digits = 1)

## interactive effects ---------------------------------------------------------  

int_effects <- summary(lm_abund_til2)$coefficients %>%
  as.data.frame() %>%
  mutate(categories = row.names(summary(lm_abund_til2)$coefficients))
  
int_res_fall <- int_effects %>%
  janitor::clean_names() %>%
  filter(categories == "treatmentRES:seasonFall") %>%
  mutate(pr_t = case_when(
    pr_t < 0.001 ~ "<0.001")
  ) %>%
  pull(pr_t)

int_mow_fall <- int_effects %>%
  janitor::clean_names() %>%
  filter(categories == "treatmentMOW:seasonFall") %>%
  mutate(pr_t = case_when(
    pr_t < 0.001 ~ "<0.001")
  ) %>%
  pull(pr_t)

(sb_abund_plot <- sb_data_viz %>%
    ggplot(aes(x = treatment, y = abund, fill = season)) +
    geom_boxplot() + 
    geom_point(
      position = position_dodge(width = 0.75),
      alpha = 0.1
    ) +
    
    # pairwise significance between mowing and undisturbed
    geom_signif(
      y_position = 300, 
      xmin = 2, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_abund_mow_res)),
      alpha = 0.5
    ) + 
    
    # pairwise significance between tilling and undisturbed
    geom_signif(
      y_position = 400, 
      xmin = 1, 
      xmax = 2, 
      annotation = paste("p = ", as.character(pairs_abund_til_res)),
      alpha = 0.5
    ) + 
    
    # pairwise significance between mowing and tilling 
    geom_signif(
      y_position = 500, 
      xmin = 1, 
      xmax = 3, 
      annotation = paste("p = ", as.character(pairs_abund_mow_til)), 
      alpha = 0.5
    ) +   
    
    # interactive effects
    annotate(
      "text", 
      x = 1.25, 
      y = 800,
      size = 3,  
      label= paste("interaction - unmown X fall:", "p", int_res_fall)) + 
    
    annotate(
      "text", 
      x = 1.20, 
      y = 750,
      size = 3, 
      label= paste("interaction - mown X fall:", "p", int_mow_fall)) + 
    
    ylim(0, 800) + 
    
    labs(
      x = "Restoration Stage", 
      y = "Seedling Emergent Abundance"
    ) + 
    
    scale_fill_manual(
      name = "Season",
      values = cbPalette) + 
    
    theme_bw() + 
    theme(
      axis.title.x = element_text(
        margin = margin(t = 10, r = 0, b = 0, l = 0)
      ),
      axis.title.y = element_text(
        margin = margin(t = 0, r = 10, b = 0, l = 0)
      ),
      text = element_text(size = 16)
    )
)

# save to disk -----------------------------------------------------------------
ggsave(
  filename = here("output", "results", "figure-2_abund.png"), 
  plot = sb_abund_plot, 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 7
)

ggsave(
  filename = here("output", "results", "figure-2_abund.pdf"), 
  plot = sb_abund_plot, 
  device = "pdf", 
  units = "in", 
  height = 5, 
  width = 7
)