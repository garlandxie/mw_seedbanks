# libraries ----
library(here)
library(dplyr)
library(ggplot2)
library(performance)
library(MASS)

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
    treatment = factor(treatment)
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
  scale_color_discrete(
    name = "Management Regime", 
    labels = c("Maintenance-Mow", "Undisturbed", "Heavily-Tilled")
  ) + 
  labs(
    x = "Site", y = "Seedling Emergent Abundance") + 
  coord_flip() + 
  theme_bw()

ggplot_sr <- sb_comm %>%
  ggplot(aes(x = site, y = species_richness, col = treatment)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.2) + 
  facet_wrap(~season, ncol = 2, nrow = 3) + 
  scale_color_discrete(
    name = "Management Regime", 
    labels = c("Maintenance-Mow", "Undisturbed", "Heavily-Tilled")
  ) + 
  labs(
    x = "Site", y = "Species Richness") + 
  coord_flip() + 
  theme_bw()

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

glm_abund_poisson <- glm(
  formula = abund ~ treatment + season + site, 
  family = poisson(link = "log"), 
  data = sb_comm
  )

# check assumption of poisson regression models 
check_overdispersion(glm_abund_poisson)

glm_abund_nb <- glm.nb(
  formula = abund ~ treatment + season + site, 
  data = sb_comm
)

summary(glm_abund_nb)

## species richness ----

glm_sr_poisson <- glm(
  formula = species_richness ~ treatment + season + site, 
  family = poisson(link = "log"), 
  data = sb_comm
)

# check assumption of poisson regression models 
check_overdispersion(glm_sr_poisson)

summary(glm_sr_poisson)

## save to disk ----
ggsave(
  filename = here("output", "results", "seedling_abundance.png"), 
  plot = ggplot_abund, 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 7
)

ggsave(
  filename = here("output", "results", "seedling_richness.png"), 
  plot = ggplot_sr, 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 7
)




  