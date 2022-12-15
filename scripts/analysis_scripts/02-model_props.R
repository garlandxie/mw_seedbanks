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
  props_exotic_logit = car::logit(props_spontan_exotic)
  )

prop_exotic_lm <- lmer(
  props_exotic_logit ~ treatment + (1|site_code),
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
    props_native_logit ~ treatment + (1|site_code),
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
  logit_props_sm ~ treatment + (1|site_code),
  data = props_sm_transform
)

summary(props_sm_lm)
performance::r2(props_sm_lm)

## proportion of invasive species ----------------------------------------------

props_inv_transform <- mutate(
  props_tidy, 
  logit_props_inv = car::logit(props_spontan_invasive)
) 

props_inv_lm <- lmer(
  logit_props_inv ~ treatment + (1|site_code),
  data = props_inv_transform
)

# pairwise comparisons ---------------------------------------------------------

exotic_emm_trt <- emmeans(
  prop_exotic_lm, 
  "treatment" 
)

spontan_native_emm_trt <- emmeans(
  props_spontan_lm, 
  "treatment" 
)

sm_emm_trt <- emmeans(
  props_sm_lm, 
  "treatment" 
)


inv_emm_trt <- emmeans(
  props_inv_lm, 
  "treatment" 
)

# confidence intervals ---------------------------------------------------------

ci_sm_trt <- confint(pairs(sm_emm_trt)) %>%
  
  # note that the response variable was logit-transformed
  # to back-transform, use exp() => see Warton and Hui 2011
  # so, the interpretation here is "predicted increase in odds
  # of x by a factor of [exp(coefficient)]"
  
  mutate(
    backtransformed_estimate = exp(estimate),
    backtransformed_lower.CL = exp(lower.CL),
    backtransfromed_upp.CL = exp(upper.CL)
  )
  
# visualize data ---------------------------------------------------------------

# rename management regimes with their full titles 
props_data_viz  <- props %>%
  mutate(
    treatment = as.character(treatment), 
    treatment = case_when(
      treatment == "MOW" ~ "Maintenance-Mow", 
      treatment == "RES" ~ "Undisturbed",
      treatment == "TIL" ~ "Tilling"
    ), 
    treatment = factor(
      treatment, levels = c("Undisturbed", "Maintenance-Mow", "Tilling")
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
     x = "Management Regime", 
     y = "Proportion of Exotic Species"
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
     x = "Management Regime", 
     y = "Proportion of Spotaneous Native Species"
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
     x = "Management Regime", 
     y = "Proportion of Species in Seed Mix"
   ) + 
   scale_fill_manual(
     name = "Season",
     values = cbPalette
     ) + 
   theme_bw() +
   theme(text = element_text(size = 15))
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
     x = "Management Regime", 
     y = "Proportion of Invasive Species"
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

## figure s3: spontaneous native species ----------------------------------------
ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s3_prop-spontan-native.png"
  ), 
  plot = prop_nativ_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)

## figure s4: exotic species ---------------------------------------------------
ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s4_prop-exotic.png"
    ), 
  plot = prop_exotic_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)

## figure s5: invasive species -------------------------------------------------
ggsave(
  filename = here(
    "output", 
    "data_appendix_output", 
    "figure-s5_prop-invasive.png"
  ), 
  plot = prop_inv_plot,
  device = "png", 
  units = "in",
  height = 4, 
  width = 6
)
