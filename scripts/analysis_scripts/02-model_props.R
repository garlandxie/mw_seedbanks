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
# Purpose: to run the linear mixed-effects models involving proportions as a
# response variable 

# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html

# library ----------------------------------------------------------------------
library(here)      # for creating relative file-paths
library(lme4)      # for running linear mixed effect models (LMM)
library(car)       # for calculating logit transformations
library(ggplot2)   # for visualizing data
library(dplyr)     # for manipulating data
library(DHARMa)    # for running diagnostic plots for LMM's
library(lmerTest)  # for obtaining p-values for LMM's
library(ggsignif)  # for plotting significant values in ggplot
library(emmeans)   # for doing pairwise comparisons

# import data ------------------------------------------------------------------

props <- read.csv(
  here("data", "intermediate_data", "props.csv"),
  row.names = 1
)

# clean data -------------------------------------------------------------------

props_tidy <- mutate(
  props, 
  season = factor(season, levels = c("Spring", "Fall")), 
  treatment = factor(treatment, levels = c("RES", "MOW", "TIL"))
  )

# regression models ------------------------------------------------------------

## proportion of exotics -------------------------------------------------------
prop_exotic_transform <- mutate(
  props_tidy, 
  props_exotic_logit = car::logit(props_spontan_exotic, percents = FALSE, adjust = 0.025)
  )

prop_exotic_lm <- lmer(
  props_exotic_logit ~ treatment + season +  (1|site_code),
  data = prop_exotic_transform
  )

summary(prop_exotic_lm)
performance::r2(prop_exotic_lm)

## proportion of spontaneous native species ------------------------------------
prop_native_transform <- mutate(
  props_tidy, 
  props_native_logit = car::logit(props_spontan_native)
)

props_spontan_lm <- lmer(
    props_native_logit ~ treatment + season + (1|site_code),
    data = prop_native_transform
  )

summary(props_spontan_lm)
performance::r2(props_spontan_lm)

## proportion of seed mix species ----
props_sm_transform <- mutate(
  props_tidy, 
  logit_props_sm = car::logit(props_sm)
  ) 

props_sm_lm <- lmer(
  logit_props_sm ~ treatment + season + (1|site_code),
  data = props_sm_transform
)

summary(props_sm_lm)
performance::r2(props_sm_lm)

## proportion of invasive species ----
props_inv_transform <- mutate(
  props_tidy, 
  logit_props_inv = car::logit(props_spontan_invasive)
) 

props_inv_lm <- lmer(
  logit_props_inv ~ treatment + season + (1|site_code),
  data = props_inv_transform
)

summary(props_inv_lm)
performance::r2(props_inv_lm)

## proportion of invasive species ----------------------------------------------

# pairwise comparisons ---------------------------------------------------------

exotic_emm_trt <- emmeans(
  prop_exotic_lm, 
  "treatment" 
) %>% pairs()

spontan_native_emm_trt <- emmeans(
  props_spontan_lm, 
  "treatment" 
) %>% pairs()

sm_emm_trt <- emmeans(
  props_sm_lm, 
  "treatment" 
) %>% pairs()


inv_emm_trt <- emmeans(
  props_inv_lm, 
  "treatment" 
) %>% pairs()

# visualize data ---------------------------------------------------------------

# rename management regimes with their full titles 
props_data_viz  <- props %>%
  mutate(
    treatment = as.character(treatment), 
    treatment = case_when(
      treatment == "MOW" ~ "Mown", 
      treatment == "RES" ~ "Unmown",
      treatment == "TIL" ~ "Newly-established"
    ), 
    treatment = factor(
      treatment, levels = c("Newly-established", "Unmown", "Mown")
    )  
  ) %>%
  mutate(season = factor(season, levels = c("Spring", "Fall")))

## color-blind friendly palette ------------------------------------------------

cbPalette <- c("#009E73", "#E69F00")

## proportion of exotics -------------------------------------------------------

(prop_exotic_plot <- props_data_viz %>%
   ggplot() +
   geom_boxplot(
     aes(x = treatment, y = props_spontan_exotic, fill = season)
   ) + 
   geom_point(
     aes(x = treatment, y = props_spontan_exotic, fill = season),
     position = position_dodge(width = 0.75),
     alpha = 0.3
   ) + 
   labs(
     x = "Restoration Stage", 
     y = "Proportion of Spontaneous Non-Native Species"
   ) + 
   scale_fill_manual(
     name = "Season",
     values = cbPalette
     ) + 
   scale_y_continuous(breaks = seq(0, 1.0, 0.25)) + 
   theme_bw() 
)

## proportion of spontaneous native species ------------------------------------

(prop_nativ_plot <- props_data_viz %>%
   ggplot() +
   geom_boxplot(
     aes(x = treatment, y = props_spontan_native, fill = season)
   ) + 
   geom_point(
     aes(x = treatment, y = props_spontan_native, fill = season),
     position = position_dodge(width = 0.75),
     alpha = 0.3
   ) + 
   labs(
     x = "Restoration Stage", 
     y = "Proportion of Spontaneous Native Species"
   ) + 
   scale_fill_manual(
     name = "Season",
     values = cbPalette
   ) + 
   scale_y_continuous(breaks = seq(0, 1.0, 0.25)) + 
   theme_bw() 
)

## proportion of species in the seed mix ---------------------------------------

(prop_sm_plot <- props_data_viz %>%
   ggplot() +
   geom_boxplot(
     aes(x = treatment, y = props_sm, fill = season)
   ) + 
   geom_point(
     aes(x = treatment, y = props_sm, fill = season),
     position = position_dodge(width=0.75),
     alpha = 0.3
     ) + 
   labs(
     x = "Restoration Stage", 
     y = "Proportion of Native Seed Mix Species"
   ) + 
   scale_fill_manual(
     name = "Season",
     values = cbPalette
     ) + 
   theme_bw() +
   theme(text = element_text(size = 13))
) 

## proportion of invasive species ----------------------------------------------

(prop_inv_plot <- props_data_viz %>%
   ggplot() +
   geom_boxplot(
     aes(x = treatment, y = props_spontan_invasive, fill = season)
   ) + 
   geom_point(
     aes(x = treatment, y = props_spontan_invasive, fill = season),
     position = position_dodge(width=0.75),
     alpha = 0.3
   ) + 
   labs(
     x = "Restoration Stage", 
     y = "Proportion of Spontaneous Invasive Species"
   ) + 
   scale_fill_manual(
     name = "Season",
     values = cbPalette
   ) + 
   theme_bw() 
) 

# save to disk -----------------------------------------------------------------

## figure 4: species in the seed mix -------------------------------------------
ggsave(
  filename = here("output", "results", "figure-4_prop-sm.png"), 
  plot = prop_sm_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)

ggsave(
  filename = here("output", "results", "figure-4_prop-sm.pdf"), 
  plot = prop_sm_plot,
  device = "pdf", 
  units = "in",
  height = 4, 
  width = 6
)

## figure s3: spontaneous native species ----------------------------------------
ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s1_prop-spontan-native.png"
  ), 
  plot = prop_nativ_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)

ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s1_prop-spontan-native.pdf"
  ), 
  plot = prop_nativ_plot,
  device = "pdf", 
  units = "in",
  height = 4, 
  width = 6
)

## figure s4: exotic species ---------------------------------------------------
ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s2_prop-exotic.png"
    ), 
  plot = prop_exotic_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)

ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s2_prop-exotic.pdf"
  ), 
  plot = prop_exotic_plot,
  device = "pdf", 
  units = "in",
  height = 4, 
  width = 6
)

## figure s5: invasive species -------------------------------------------------
ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s3_prop-invasive.png"
  ), 
  plot = prop_inv_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)

ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s3_prop-invasive.pdf"
  ), 
  plot = prop_inv_plot,
  device = "pdf", 
  units = "in",
  height = 4, 
  width = 6
)
