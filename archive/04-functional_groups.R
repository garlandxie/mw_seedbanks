# libraries ----
library(here)
library(dplyr)
library(ggplot2)
library(rdryad)
library(svglite)

# import ----

## seed bank data ----

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

## plants of toronto ----

dryad_doi <- "10.5061/dryad.1ns1rn8sg"
dryad_link <- rdryad::dryad_download(dryad_doi)
plants_to <- read.csv(unlist(dryad_link))

## taxonomy info ----

taxon <- read.csv(
  here("data", "input_data", "seed_bank_taxonomy.csv")
)

# data clean ----

# some prep
sb <- rbind(fall_sb, spr_sb)
taxon <- mutate(taxon, binom_latin = paste(Genus, Species))

# need to clean more of this data!!!
sb_tidy <- sb %>%
  left_join(taxon, by = c("spp_code" = "Code")) %>%
  left_join(plants_to, by = c("binom_latin" = "SCIENTIFIC_NAME")) %>%
  dplyr::select(
    season, 
    site_code, treatment, binom_latin, Growth_form, total_abund)

sb_tidy <- sb_tidy %>%
  
    # manually assign growth form based on Plants USDA website
    mutate(Growth_form = case_when(
      binom_latin == "Acalypha rhomboidea" ~ "Herbaceous", 
      binom_latin == "Carex spp." ~ "Graminoid",
      binom_latin == "Conyza canadensis" ~ "Herbaceous", 
      binom_latin == "Euphorbia serpyllifolia" ~ "Herbaceous", 
      binom_latin == "Euphorbia spp." ~ "Herbaceous", 
      binom_latin == "Lactuca sativa" ~ "Herbaceous", 
      binom_latin == "Mentha arvensis" ~ "Herbaceous", 
      binom_latin == "Panicum milliaceum" ~ "Herbaceous", 
      binom_latin == "Picris echioides" ~ "Herbaceous", 
      binom_latin == "Pycnanthemum muticum" ~ "Herbaceous", 
      binom_latin == "Sinapis alba" ~ "Herbaceous", 
      binom_latin == "Verbenum urticifolia" ~ "Herbaceous", 
      binom_latin == "Veronica agrestis" ~ "Herbaceous", 
      TRUE ~ Growth_form
    )
  )
 
sb_tidy <- sb_tidy %>%
  group_by(treatment, site_code, Growth_form) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE))

# plot ----

(func_groups <- sb_tidy %>% 
  filter(!is.na(Growth_form)) %>%
  mutate(site = factor(site_code)) %>%
  arrange(site) %>%
  ggplot(aes(x = Growth_form, y = abundance, fill = treatment)) +
  geom_col() +
  labs(x = "Growth Form") + 
  scale_fill_discrete(
    name = "Management Regimes", 
    labels = c("Maintenance-Mow", "Undisturbed", "Tilled")
  ) + 
  facet_wrap(~site) + 
  coord_flip() + 
  theme_bw() 
)

# save to disk ----

ggsave(
  filename = here("output", "results", "func_groups.svg"), 
  plot = func_groups, 
  device = "svg", 
  units = "in", 
  height = 5, 
  width = 7
)
