################################################################################
# Accompanying code for the followng research project: 
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
#              nicholas.sookhan@mail.utoronto.ca
#              scott.macivor@mail.utoronto.ca
#
# Purpose of this R script: to perform regression models 
#
# IMPORTANT: Please refresh your R session before you run this script
# Why? See https://rstats.wtf/save-source.html

# libraries ----
library(here)         # for creating relative file-paths
library(lme4)         # for analyzing GLMM's
library(lmerTest)     # for analyzing GLMM's
library(ggplot2)      # for visualizing data 
library(DHARMa)       # for checking GLMM assumptions
library(visdat)       # for checking missing values
library(forcats)      # for manipulating factor variables
library(performance)  # for testing model performance
library(dplyr)        # for manipulating data 

# import ----

# community abundance
df_comm_abund <- 
  read.csv(
    here("data", "final", "sb_comm_abund.csv"),
    stringsAsFactors = FALSE
  )

# species richness
df_sr <- 
  read.csv(
    here("data", "final", "sb_sr.csv"),
    stringsAsFactors = FALSE
  )

# check packaging ----

str(df_comm_abund)
visdat::vis_miss(df_comm_abund)

# statistical inference: community abundance  ----

glmer_tot_abund <- glmer(
  total_abund ~ Treatment + (1|Site),
  family = poisson,
  data = df_comm_abund
  ) 

# check model assumptions
sim_glmm_TA <- DHARMa::simulateResiduals(fittedModel = glmer_tot_abund)
plot(sim_glmm_TA)

# statistical inference: species richness ----

glmer_sr <- glmer(
  species_richness ~ Treatment + (1|Site),
  family = poisson,
  data = df_sr
) 

# check model assumptions
sim_glmm_SR <- DHARMa::simulateResiduals(fittedModel = glmer_sr)
plot(sim_glmm_SR)

# variance partitioning ----

r2_glmm_TA <- performance::r2(glmer_tot_abund)
r2_glmm_SR <- performance::r2(glmer_sr)

# plots ----

df_sb <- df_comm_abund %>%
  inner_join(df_sr, by = c("Section", "Site", "Treatment", "Plot")) %>%
  select(Section, 
         Site, 
         Treatment, 
         Plot, 
         total_abund, 
         species_richness)

# community abundance

(plot_abund <- df_sb %>%
  mutate(
    log_tot_ab = log(total_abund),
    
    Section = case_when(
      Section == "S2" ~ "Section 2", 
      Section == "S4" ~ "Section 4",
      TRUE ~ Treatment),
    
    Treatment = case_when(
      Treatment == "TIL" ~ "Seed Drill", 
      Treatment == "MOW" ~ "Maintenance Mow", 
      Treatment == "RES" ~ "Undisturbed", 
      TRUE ~ Treatment ),
    Treatment = factor(Treatment), 
    Treatment = fct_relevel(
      Treatment, 
      "Seed Drill", "Maintenance Mow", "Undisturbed"
      )
    ) %>% 
  ggplot(aes(x = Site, y = log_tot_ab, fill = Treatment)) + 
  geom_boxplot() + 
  geom_point(alpha = 0.2)  + 
  coord_flip() + 
  labs(
    y = "log(Community Abundance)", 
    x = "Site"
  ) +
  scale_fill_manual(
    name   = c("Seed Drill", "Maintenance Mow", "Undisturbed"),
    values = c("#56B4E9", "#F0E442", "#E69F00")
  ) + 
  facet_wrap(Treatment~Section) + 
  theme_bw() + 
  theme(legend.position = "none")
)

# species richness

(plot_sr <- df_sr %>%
    mutate(
      
      # sections
      Section = case_when(
        Section == "S2" ~ "Section 2", 
        Section == "S4" ~ "Section 4",
        TRUE ~ Treatment),
      
      # treatments
      Treatment = case_when(
        Treatment == "TIL" ~ "Seed Drill", 
        Treatment == "MOW" ~ "Maintenance Mow", 
        Treatment == "RES" ~ "Undisturbed", 
        TRUE ~ Treatment),
      
      Treatment = factor(Treatment), 
      
      Treatment = fct_relevel(
        Treatment, 
        "Seed Drill", "Maintenance Mow", "Undisturbed"
      )
    ) %>%
    ggplot(aes(x = Site, y = species_richness, fill = Treatment)) + 
    geom_boxplot() + 
    geom_point(alpha = 0.2)  + 
    coord_flip() + 
    labs(
      y = "Species Richness", 
      x = "Site"
    ) +
    facet_wrap(Treatment~Section) + 
    scale_fill_manual(
      name   = c("Seed Drill", "Maintenance Mow", "Undisturbed"),
      values = c("#56B4E9", "#F0E442", "#E69F00")
    ) + 
    theme_bw() +
    theme(legend.position = "none")
)

# save to disk ----

ggsave(
  plot = plot_abund, 
  filename = here("figures", "fig1_comm-abund_by_site.png"), 
  device   = "png", 
  units    = "in", 
  height   = 5, 
  width    = 8
)

ggsave(
  plot = plot_sr, 
  filename = here("figures", "fig2_sr_by_site.png"), 
  device   = "png", 
  units    = "in", 
  height   = 5, 
  width    = 8
)

