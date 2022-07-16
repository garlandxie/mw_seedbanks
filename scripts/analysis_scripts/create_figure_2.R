library(here)     # for creating relative file-paths
library(dplyr)    # for manipulating data
library(ggplot2)  # for visualizing data
library(lme4)     # for running regression models
library(emmeans)  # for doing pairwise comparisons
library(tidyr)    # for tidying data
library(ggsignif) # for visualizing pairwise comparisons

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

glmer_richness <- glmer(
  formula = species_richness ~ treatment + season + (1|site_code), 
  family = poisson(link = "log"),
  data = sb_comm
)

# pairwise comparisons ----

rich_emm_trt <- emmeans(
  glmer_richness, 
  "treatment", 
  lmer.df = "satterthwaite"
)

# get summary of comparisons, coefficients, and p-values
pairs_rich_trt <- as.data.frame(pairs(rich_emm_trt))

# obtain p-value for comparison between 
# undisturbed and tilling
pairs_til_res <- pairs_rich_trt %>%
  filter(contrast == "RES - TIL") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ 0.001)
  ) %>%
  pull(p.value) 

# obtain p-value for comparison between 
# maintenance-mowing and undisturbed
pairs_mow_res <- pairs_rich_trt %>%
  filter(contrast == "RES - MOW") %>%
  pull(p.value) %>%
  signif(digits = 1)

# obtain p-value for comparison between
# maintenance-mowing and tilling
pairs_mow_til <- pairs_rich_trt %>%
  filter(contrast == "MOW - TIL") %>%
  # use p < 0.001 if the p-value is really small 
  mutate(p.value = case_when(
    p.value < 0.001 ~ 0.001)
  ) %>%
  pull(p.value) 

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

# save to disk ----

ggsave(
  filename = here("output", "results", "figure-2.svg"), 
  device = "svg", 
  units = "in",
  height = 5, 
  width = 7
)
