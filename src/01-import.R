# libraries ----
library(here)
library(googlesheets4)
library(visdat)
library(ggplot2)
library(dplyr)
library(forcats)
library(janitor)

# import ----

link <- "https://docs.google.com/spreadsheets/d/1O1Ll_PsW3qKwdZ_xnTrDKT_kGnpLGvtUKBQ73zvvBBM/edit?usp=sharing"
df_sb <- googlesheets4::read_sheet(link, sheet = "raw_data")

# check packaging ----
str(df_sb)
head(df_sb, n = 10)
tail(df_sb, n = 10)

# check missing data ----
visdat::vis_miss(df_sb)
visdat::vis_dat(df_sb)

# clean -----

df_tidy <- df_sb %>%
  janitor::clean_names() %>%
  mutate(spp_code = case_when(
    spp_code == "FRVE"  ~ "PONO", 
    spp_code == "ERCA"  ~ "COCA", 
    spp_code == "PACA"  ~ "PAMI", 
    spp_code == "AMSP"  ~ "VETH", 
    spp_code == "CAREX" ~ "ANGE", 
    spp_code == "SEVI"  ~ "ANGE",
    TRUE ~ spp_code)
  )
  
# plot ----

df_comm_abund <- df_sb %>%
  group_by(Section, Site, Treatment, Plot) %>%
  summarize(total_abund = sum(Abund, na.rm = TRUE))

abund_by_plot <- df_sb1 %>%
  ggplot(aes(x = Plot, y = total_abund)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Site) +
  labs(x = NULL, y = "Total Abundance") + 
  coord_flip() +
  theme_bw()

df_sb2 <- df_sb %>%
  group_by(Spp_code) %>%
  summarize(total_abund = sum(Abund, na.rm = TRUE))

(abund_by_spp <- df_sb2 %>%
  mutate(
    Spp_code = factor(Spp_code),
    Spp_code = fct_reorder(Spp_code, total_abund)) %>%
  ggplot(aes(x = Spp_code, y = total_abund)) + 
  geom_bar(stat = "identity") +
  labs(x = NULL, y = "Abundance") + 
  coord_flip() +
  theme_bw()
)

df_sr <- df_sb %>%
  group_by(Section, Site, Treatment, Plot) %>%
  summarize(species_richness = length(unique(Spp_code)))

# Useful plots for Devlin's 5 min ppt ------------------------------------------

(comm_abund_mgmt_regime <- df_tidy %>%
  group_by(section, site, treatment, plot) %>%
  summarize(total_abund = sum(abund, na.rm = TRUE)) %>%
  ggplot(aes(x = treatment, y = total_abund, color = treatment)) + 
  geom_boxplot(show.legend = FALSE) + 
  geom_point(show.legend = FALSE, alpha = 0.1) + 
  scale_x_discrete(labels=c("Maintenance Mow", "Undisturbed", "Seed Drilling")) + 
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9")) + 
  labs(
    x = "Management Regimes",
    y = "Community Abundance") + 
  theme_classic()
)

(sr_mgmt_regime <- df_tidy %>%
    group_by(section, site, treatment, plot) %>%
    summarize(sr = length(unique(spp_code))) %>%
    ggplot(aes(x = treatment, y = sr, color = treatment)) + 
    geom_boxplot(show.legend = FALSE) + 
    geom_point(show.legend = FALSE, alpha = 0.1) + 
    scale_x_discrete(labels=c("Maintenance Mow", "Undisturbed", "Seed Drilling")) + 
    scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9")) + 
    labs(
      x = "Management Regimes",
      y = "Species Richness") + 
    theme_classic()
)

(df_sb6 <- df_tidy %>%
    group_by(section, site, treatment, plot) %>%
    summarize(total_abund = sum(abund, na.rm = TRUE)) %>%
    ggplot(aes(x = site, y = total_abund, fill = treatment)) + 
    geom_boxplot() + 
    geom_point(alpha = 0.1) + 
    coord_flip() + 
    theme_bw()
)

(df_sb7 <- df_tidy %>%
    group_by(section, site, treatment, plot) %>%
    summarize(sr = length(unique(spp_code))) %>%
    ggplot(aes(x = site, y = sr, fill = treatment)) + 
    geom_boxplot() + 
    geom_point(alpha = 0.1) + 
    coord_flip() + 
    theme_bw()
)


# save to disk: plots ----------------------------------------------------------

ggsave(
  filename = here("figures", "abund_by_plot.png"),
  plot = abund_by_plot, 
  device = "png", 
  units = "in", 
  width = 5,
  height = 5
)

ggsave(
  filename = here("figures", "abund_by_spp.png"),
  plot = abund_by_spp, 
  device = "png", 
  units = "in", 
  width = 5,
  height = 5
)

ggsave(
  filename = here("figures", "comm_abund_by_mgmt_regime.png"),
  plot = comm_abund_mgmt_regime,
  device = "png", 
  units = "in", 
  width = 5,
  height = 5
)

ggsave(
  filename = here("figures", "sr_mgmt_regime.png"),
  plot = sr_mgmt_regime,
  device = "png", 
  units = "in", 
  width = 5,
  height = 5
)

# save to disk: csv's ----------------------------------------------------------

write.csv(
  x = df_comm_abund, 
  file = here("data", "final", "sb_comm_abund.csv"), 
  row.names = FALSE
)

write.csv(
  x = df_sr, 
  file = here("data", "final", "sb_sr.csv")
)
