# libraries ----
library(here)
library(googlesheets4)
library(ggplot2)
library(janitor)
library(dplyr)
library(ggrepel)

# import ----

link <- "https://docs.google.com/spreadsheets/d/1O1Ll_PsW3qKwdZ_xnTrDKT_kGnpLGvtUKBQ73zvvBBM/edit?usp=sharing"
df_sb <- googlesheets4::read_sheet(link, sheet = "raw_data")

# data clean ----

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

df_sum <- df_tidy %>%
  group_by(site, plot, date_ymd) %>%
  summarize(total_abund = sum(abund, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(days_apart = as.numeric(as.POSIXct(Sys.Date()) - date_ymd))

df_sum %>%
  filter(site %in% c("KENN", "TIMH", "VICP")) %>%
  ggplot(aes(x = days_apart, y = total_abund)) +
  geom_point() + 
  geom_text_repel(aes(label = plot)) + 
  facet_wrap(~site) + 
  coord_flip() +
  theme_bw()


  
