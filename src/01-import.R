# libraries ----
library(here)
library(googlesheets4)
library(visdat)
library(ggplot2)
library(dplyr)
library(forcats)

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

# plot ----

df_sb1 <- df_sb %>%
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

# save to disk ----

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
