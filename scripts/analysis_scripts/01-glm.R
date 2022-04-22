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
library(performance)
library(MASS)
library(patchwork)
library(lme4)
library(emmeans)
library(tidyr)

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
       "meadoway_seed_bank_taxonomy.csv")
)

# check packaging ----
glimpse(spr_sb)
glimpse(fall_sb)

# clean data ----

sb <- rbind(fall_sb, spr_sb)

## table ----

# summarize: species as rows, site as columns, cells as abundances
table_summary <- sb %>%
  group_by(site, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = site, values_from = abundance) %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  left_join(sb_taxon, by = c("spp_code" = "Code")) %>%
  dplyr::select(
    spp_code, Genus, Species, 
    VICP, TIMH, KENN, 
    GRNB, BNSH, DAVE, 
    DVAG, AMBJ, BRIM
    )

write.csv(
  x = table_summary, 
  here("data", "analysis_data", "table_devlin.csv"), 
)

## total germination ----

# spring and fall germination total emergence for each treatment
total_germ_by_trt <- sb %>%
  group_by(season, treatment) %>%
  summarize(abund = sum(total_abund, na.rm = TRUE)) %>%
  ungroup()

write.csv(
  x = total_germ_by_trt, 
  file = here("data", "intermediate_data", "total_germ_by_sn_trt.csv")
)

## obtain species identity for plot outliers -----

### BRIM Spring 11 ----
out_brim_spr <- sb %>%
  filter(site == "BRIM" & season == "Spring") %>%
  group_by(season, site, treatment, plot) %>%
  summarize(abund = sum(total_abund, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(desc(abund))

out_brim_spr_11 <- sb %>%
  filter(site == "BRIM" & season == "Spring" & plot == 11) 

write.csv(
  x = out_brim_spr_11, 
  file = here("data", "intermediate_data", "out_brim_spr_11.csv")
)

### BNSH Spring 15 ----
out_bnsh_spr <- sb %>%
  filter(site == "BNSH" & season == "Spring") %>%
  group_by(season, site, treatment, plot) %>%
  summarize(abund = sum(total_abund, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(desc(abund))

out_bnsh_15 <- sb %>%
  filter(site == "BNSH" & season == "Spring" & plot == 15) 

write.csv(
  x = out_bnsh_15, 
  file = here("data", "intermediate_data", "out_bnsh_spr_15.csv")
)

### BRIM Fall 11 ----

out_brim_fall <- sb %>%
  filter(site == "BRIM" & season == "Fall") %>%
  group_by(season, site, treatment, plot) %>%
  summarize(abund = sum(total_abund, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(desc(abund))

out_brim_fall_11 <- sb %>%
  filter(site == "BRIM" & season == "Fall" & plot == 11) 

write.csv(
  x = out_brim_fall_11, 
  file = here("data", "intermediate_data", "out_brim_fall_11.csv")
)

### AMBJ Fall ----

out_ambj_fall <- sb %>%
  filter(site == "AMBJ" & season == "Fall") %>%
  group_by(season, site, treatment, plot) %>%
  summarize(abund = sum(total_abund, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(desc(abund))

out_ambj_fall_23 <- sb %>%
  filter(site == "AMBJ" & season == "Fall" & plot == 23) 

write.csv(
  x = out_ambj_fall_23, 
  file = here("data", "intermediate_data", "out_ambj_fall_23.csv")
)

### TIMH Fall ----

out_timh_fall <- sb %>%
  filter(site == "TIMH" & season == "Fall") %>%
  group_by(season, site, treatment, plot) %>%
  summarize(abund = sum(total_abund, na.rm = TRUE)) %>% 
  ungroup() %>%
  arrange(desc(abund))

out_timh_fall_5 <- sb %>%
  filter(site == "TIMH" & season == "Fall" & plot == 5) 

write.csv(
  x = out_timh_fall_5, 
  file = here("data", "intermediate_data", "out_timh_fall_5.csv")
)

# get community-level abundance and species richness
sb_comm <- sb %>%
  group_by(season, section, site, treatment, plot) %>%
  summarize(
    abund = sum(total_abund, na.rm = TRUE),
    species_richness = dplyr::n_distinct(spp_code)
    ) %>%
  ungroup() %>%
  mutate(
    season = factor(season, levels = c("Spring", "Fall")), 
    site = factor(
      site, 
      levels = c(
        "VICP", "TIMH", "KENN", 
        "GRNB", "DAVE", "BNSH", 
        "DVAG", "BRIM", "AMBJ")
      ), 
    treatment = factor(treatment, levels = c("RES", "MOW", "TIL"))
  )

# exploring data (before regression) -----

## step 1: check for outliers  ----

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

## step 2: check for homogeneity of variance ----

# multi-panel conditional box-plots
ggplot_abund <- sb_comm %>%
  ggplot(aes(x = site, y = abund, col = treatment)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~season, ncol = 2, nrow = 3) + 
  labs(x = "Site", y = "Seedling Emergent Abundance") + 
  coord_flip() + 
  theme_bw() +
  theme(legend.position = "none")

ggplot_sr <- sb_comm %>%
  ggplot(aes(x = site, y = species_richness, col = treatment)) + 
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

multi_plot <- ggplot_abund + ggplot_sr

# NOTE: make sure you check the residuals of the regression models

## step 3: check for normally-distributed data ----

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

# run regression models (with outliers) ----

## abundance ----

glmer_abund_poisson <- glmer(
  formula = abund ~ season + treatment + (1|site), 
  family = poisson(link = "log"), 
  data = sb_comm
  )

### check overdispersion ----
check_overdispersion(glmer_abund_poisson)

### account for overdispersion ----
glmer_abund_nb <- glmer.nb(
  formula = abund ~ treatment + season + (1|site), 
  data = sb_comm
)

### model summary ----
summary(glmer_abund_nb)

### check diagnostics ----
check_model(glmer_abund_nb)

### r-squared ----
r2(glmer_abund_nb)

### pairwise comparison ----
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
abund_pairs_trt <- abund_emm_trt %>%
  pairs() %>%
  as.data.frame()

abund_pairs_sn <- abund_emm_sn %>%
  pairs() %>%
  as.data.frame()

## species richness ----

glmer_richness <- glmer(
  formula = species_richness ~ treatment + season + (1|site), 
  family = poisson(link = "log"),
  data = sb_comm
)

### check overdispersion ----
check_overdispersion(glmer_richness)

### check model diagnostics ----

# homogeneity of variance
# collinearity
# influential observations
# normal of residuals 
# normality of random effects
check_model(glmer_richness)

### model summary -----
summary(glmer_richness)

### r squared ----
r2(glmer_richness)

### pairwise comparisons ----
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

## save to disk ----
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




  