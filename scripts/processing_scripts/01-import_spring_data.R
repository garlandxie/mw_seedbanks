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

# data clean -----

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
  
#
