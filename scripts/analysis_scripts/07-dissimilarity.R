# libraries -----
library(vegan)    # for doing community ecology analyses
library(ggplot2)  # for visualizing data
library(dplyr)    # for manipulating data
library(here)     # for creating relative file-paths

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
str(spr_sb)
str(fall_sb)

# clean data ----

# merge fall and spring season 
sb <- rbind(fall_sb, spr_sb)

# create spring community data matrix 
spring_comm_matrix <- sb %>%
  filter(season == "Spring") %>%
  group_by(treatment, site_code, plot, spp_code) %>%
  summarize(abundance = sum(total_abund, na.rm = TRUE)) %>%
  ungroup() %>%
  tidyr::pivot_wider(
    id_cols = c("treatment", "site_code", "plot"),  
    names_from = spp_code, 
    values_from = abundance
  ) %>%
  mutate(across(.cols = everything(), ~ tidyr::replace_na(.x, 0))) %>%
  as.data.frame()

rownames(spring_comm_matrix) <- paste(
  spring_comm_matrix$treatment, 
  c(1:length(spring_comm_matrix$treatment))
)

# calculate dissimilarities ----

# use Bray-Curtis dissimilarity to account for common and rare species
# common and rare species are determined by abundance (# of seedling emergents)
bray_spr <- spring_comm_matrix %>%
  select(-c("treatment", "site_code", "plot"))

bray_spr_dist <- vegdist(bray_spr, index = "bray")

# do principal coordinate analysis
pcoa_spr <- cmdscale(
  bray_spr_dist, 
  k = (nrow(bray_spr)-1),
  add = TRUE, # Cailliez correction
  eig = TRUE)

# visualize dissimilarities ----
pcoa_spr_df <- scores(pcoa_spr) %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "treatment") %>%
  select(treatment, Dim1, Dim2) %>%
  mutate(treatment = stringr::str_sub(treatment, start = 1L, end = 3L))

(pcoa_spr_plot <- pcoa_spr_df %>%
  ggplot(aes(x = Dim1, y = Dim2, col = treatment)) + 
  geom_point() + 
  scale_color_discrete(
    name = "Management Regime", 
    labels = c("Maintenance-Mow", "Undisturbed", "Tilled")
  ) + 
  labs(
    title = "a) Spring season",
    x = "PCoA Dim1", 
    y = "PCoA Dim2"
    ) + 
  theme_bw()
)


#beta_dist <- vegdist(t(obj$data$otu_rarefied[, sample_data$SampleID]),
#index = "bray")

#mds_data <- as.data.frame(mds$points)

#mds_data$SampleID <- rownames(mds_data)
#mds_data <- dplyr::left_join(mds_data, sample_data)

#ggplot(mds_data, aes(x = MDS1, y = MDS2, color = Type)) +
#  geom_point()