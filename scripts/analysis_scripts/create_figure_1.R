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


# libraries ----
library(here)   
library(dplyr)
library(ggplot2)
library(MASS)
library(patchwork)
library(lme4)
library(emmeans)
library(tidyr)
library(multcompView)
library(multcomp)
library(ggsignif)

# import ----

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

# check packaging ----
glimpse(spr_sb)
glimpse(fall_sb)

# clean data ----

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
  ) 

# regression model ----

# accounting for overdispersion
glmer_abund_nb <- glmer.nb(
  formula = abund ~ treatment + season + (1|site_code), 
  data = sb_comm
)

# pairwise comparisons ----

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

# check for sample size ----
sb_comm %>%
  group_by(season, treatment) %>%
  summarize(
    sample_size = n()
  ) %>%
  ungroup()

# plot ----

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
      treatment, levels = c("Tilling", "Undisturbed", "Maintenance-Mow")) 
  ) 

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
    annotation = paste("p = ", as.character(pairs_mow_res)), 
    alpha = 0.5
    ) +   
    
  # pairwise significance between mowing and tilling  
  geom_signif(
    y_position = 500, 
    xmin = 1, 
    xmax = 3, 
    annotation =  paste("p = <", as.character(pairs_mow_til)),
    alpha = 0.5
    ) + 
    
  ylim(0, 800) + 
    
  labs(
    x = "Management Regime", 
    y = "Seedling Emergent Abundance"
    ) + 
  theme_bw() 
)

# save to disk ---

ggsave(
  filename = here("output", "results", "figure-1.png"), 
  device = "png", 
  units = "in",
  height = 5, 
  width = 6
)

