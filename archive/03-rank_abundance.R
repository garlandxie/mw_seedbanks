# libraries ----
library(here)
library(vegan)
library(BiodiversityR)
library(dplyr)
library(ggplot2)
library(svglite)

# import ----

spr_sb <- read.csv(
  here("data", "analysis_data", "spring_seedbank.csv"),
  row.names = 1
)

fall_sb <- read.csv(
  here("data", "analysis_data", "fall_seedbank.csv"),
  row.names = 1
)

# clean data ----

sb <- rbind(fall_sb, spr_sb)

comm_matrix <- sb %>%
  group_by(treatment, site, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = abundance) %>%
  ungroup() %>%
  mutate(across(.cols = everything(), ~ tidyr::replace_na(.x, 0)))

# get rank-abundance curves ----

## maintenance mow ----
comm_mow <- sb %>%
  filter(treatment == "MOW") %>%
  group_by(site, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = abundance) %>%
  ungroup() %>%
  mutate(across(.cols = everything(), ~ tidyr::replace_na(.x, 0))) %>%
  select(-site) %>%
  as.data.frame()

ra_mow <- as.data.frame(rankabundance(x = comm_mow))

## undisturbed sites ----

comm_res <- sb %>%
  filter(treatment == "RES") %>%
  group_by(site, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = abundance) %>%
  ungroup() %>%
  mutate(across(.cols = everything(), ~ tidyr::replace_na(.x, 0))) %>%
  select(-site) %>%
  as.data.frame()

ra_res <- as.data.frame(rankabundance(x = comm_res))

## tilled sites ----

comm_til <-  sb %>%
  filter(treatment == "TIL") %>%
  group_by(site, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(names_from = spp_code, values_from = abundance) %>%
  ungroup() %>%
  mutate(across(.cols = everything(), ~ tidyr::replace_na(.x, 0))) %>%
  select(-site) %>%
  as.data.frame()

ra_til <- as.data.frame(rankabundance(x = comm_til))

# plot rank-abundance across management regimes -----
ra_mgt <- ggplot() + 
  
  geom_point(
    aes(x = rank, y = abundance, color = 'red'),
    alpha = 0.1,
    data = ra_mow
    )+ 
  
  geom_line(
    aes(x = rank, y = abundance, color = 'red'), 
    position = position_jitter(width=0.5),
    data = ra_mow
    ) + 
  
  geom_point(
    aes(x = rank, y = abundance, color = 'green'), 
    alpha = 0.1, 
    data = ra_res
    ) + 
  
  geom_line(
    aes(x = rank, y = abundance, color = 'green'), 
    position = position_jitter(width=0.5),
    data = ra_res
    ) + 
  
  geom_point(
    aes(x = rank, y = abundance, color = 'blue'), 
    alpha = 0.1, 
    data = ra_til
    ) + 
  
  geom_line(
    aes(x = rank, y = abundance, color = 'blue'),
    data = ra_til
    ) + 
  
  scale_colour_manual(
    name = "Management Regime", 
    values = c('red' = 'red', 'green' = 'green', 'blue' = 'blue'),
    labels = c("Maintenance-Mow", "Undisturbed", "Heavily_Tilled")
  ) + 
  
  theme_bw()

# save to disk -----

ggsave(
  filename = here("output/results/ra_mgt_regime.svg"), 
  plot = ra_mgt, 
  device = "svg", 
  units = "in", 
  height = 5, 
  width = 7
  )

write.csv(
  x = ra_res, 
  file = here("output", "results", "ra_res.csv")
)

write.csv(
  x = ra_til, 
  file = here("output", "results", "ra_til.csv")
)
  
write.csv(
  x = ra_mow, 
  file = here("output", "results", "ra_mow.csv")
)
