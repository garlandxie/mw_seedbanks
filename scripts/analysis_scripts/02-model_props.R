# library ----------------------------------------------------------------------
library(here)     # for creating relative file-paths
library(lme4)     # for running linear mixed effect models (LMM)
library(car)      # for calculating logit transformations
library(ggplot2)  # for visualizing data
library(dplyr)    # for manipulating data
library(DHARMa)   # for running diagnostic plots for LMM's
library(lmerTest) # for obtaining p-values for LMM's

# import data ------------------------------------------------------------------

props <- read.csv(
  here("data", "intermediate_data", "props.csv"),
  row.names = 1
)

# exploratory data visualization -----------------------------------------------

props_data_viz  <- props %>%
  mutate(
    treatment = as.character(treatment),
    treatment = case_when(
      treatment == "MOW" ~ "Maintenance-Mow", 
      treatment == "TIL" ~ "Tilling", 
      treatment == "RES" ~ "Undisturbed", 
      TRUE ~ treatment),
    treatment = factor(treatment)
  ) %>%
  mutate(season = factor(season, levels = c("Spring", "Fall")))

## proportion of spontaneous exotics -------------------------------------------

(prop_spontan_exotic <- props_data_viz %>%
   ggplot(aes(x = treatment, y = props_spontan_exotic, fill = treatment)) +
   geom_boxplot() + 
   geom_point(alpha = 0.5) + 
   labs(x = NULL, y = "Proportion of Spotanteous Exotic Species") + 
   theme_bw() + 
   facet_wrap(~season) + 
   theme(legend.position = "none")
)

## proportion of spontaneous natives ------------------------------------------- 

(prop_spotan_native <- props_data_viz %>%
   ggplot(aes(x = treatment, y = props_spontan_native, fill = treatment)) +
   geom_boxplot() + 
   geom_point(alpha = 0.5) + 
   labs(x = NULL, y = "Proportion of Spontaneous Native Species") + 
   facet_wrap(~season) + 
   theme_bw() + 
   theme(legend.position = "none")
)

## proportion of seed mix species ----------------------------------------------

(prop_sm <- props_data_viz %>%
   ggplot(aes(x = treatment, y = props_sm, fill = treatment)) +
   geom_boxplot() + 
   geom_point(alpha = 0.5) + 
   facet_wrap(~season) + 
   ylim(0, 1) + 
   labs(x = NULL, y = "Proportion of Native Species from the Seed Mix") + 
   theme_bw() + 
   theme(legend.position = "none")
)

# regression models ------------------------------------------------------------

## proportion of exotics -------------------------------------------------------
prop_exotic_transform <- mutate(
  props, 
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
  props, 
  props_native_logit = car::logit(props_spontan_native)
)

props_spontan_lm <- lmer(
    props_native_logit ~ treatment + (1|site_code),
    data = prop_native_transform
  )

summary(props_spontan_lm)

## proportion of seed mix species ----
props_sm_transform <- mutate(
  props, 
  logit_props_sm = car::logit(props_sm)
  ) 

props_sm_lm <- lmer(
  logit_props_sm ~ treatment + (1|site_code),
  data = props_sm_transform
)

## proportion of invasive species ----------------------------------------------

props_inv_transform <- mutate(
  props, 
  logit_props_inv = car::logit(props_spontan_invasive)
) 

props_inv_lm <- lmer(
  logit_props_inv ~ treatment + (1|site_code),
  data = props_inv_transform
)