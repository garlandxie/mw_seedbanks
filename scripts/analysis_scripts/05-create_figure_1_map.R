# library ----------------------------------------------------------------------
library(ggplot2)             # for visualizing plots
library(ggsn)                # for creating scale-bars
library(sf)                  # for manipulating geospatial data
library(here)                # for creating file-paths
library(dplyr)               # for manipulating data
library(gghighlight)         # for highlighting certain ggplot features
library(opendatatoronto)     # for downloading the TO municipal boundary
library(cowplot)             # for creating inset maps
library(stringr)             # for manipulating string characters

# import data ------------------------------------------------------------------

## source functions ------------------------------------------------------------

source(here("scripts", "analysis_scripts", "functions.R"))

## plot coordinates ------------------------------------------------------------

# list all plot coordinates for every surveyed site
# there should be nine sites; three for each management regime
site_file_paths <- list.files(
  here("data", "input_data", "plot_coordinates")
)

# row bind all data frames together
site_df <- dplyr::bind_rows(lapply(site_file_paths, read_file))

## meadoway shape file ---------------------------------------------------------
mw_shp <- st_read(
  here("data", "input_data", "meadoway_shapefile", "Meadoway.shp")
)

## toronto boundary shape file -------------------------------------------------

# TO boundary 
reg_bound_ID <- "841fb820-46d0-46ac-8dcb-d20f27e57bcc"
packages <- show_package(reg_bound_ID)
resources <- list_package_resources(packages)
bound <- get_resource(resources)

# clean data -------------------------------------------------------------------

## convert coordinates to decimal degrees --------------------------------------

# fyi: I regret using decimal degree minutes for data analysis
# (but it was helpful to find the plots!)
# sorry about this convoluted code!!
site_tidy <- site_df %>%
  
  # get the center plot of each site
  # for the sake of simple data visualization
  filter(Plot == 13) %>% 
  
  # convert from decimal degree minutes to decimal degrees
  # using a custom function
  mutate(
    Latitude  = sapply(Latitude, ddm_to_dd),
    Longitude = sapply(Longitude, ddm_to_dd)*-1
    ) 

## ensure consistent coordinate system  ----------------------------------------
site_crs <- st_as_sf(
  site_tidy, 
  coords = c("Longitude", "Latitude"), 
  crs = 4326
  )

mw_shp_tidy <- mw_shp %>%
  mutate(section = factor(section)) %>%
  st_transform(crs = 4326)

bound_tidy <- bound %>%
  st_transform(crs = 4326)

# visualize data ---------------------------------------------------------------

## toronto boundary ------------------------------------------------------------
to <- ggplot() + 
  geom_sf(fill = NA, data = bound_tidy)  +
  geom_sf(fill = NA, color = "red", data = mw_shp_tidy) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(), 
    
    axis.title.x  = element_blank(),
    axis.text.x   = element_blank(),
    axis.ticks.x  = element_blank(),
    
    axis.title.y  = element_blank(), 
    axis.text.y   = element_blank(), 
    axis.ticks.y  = element_blank()
    )
  
## meadow ----------------------------------------------------------------------

mw <- ggplot() + 
  
  # meadoway sections, but just highlight the tilled and restored areas
  geom_sf(aes(fill = section), data = mw_shp_tidy) +
  gghighlight(section %in% c("2", "4")) +
  
  # map site coordinates (represented as the center plot of each site)
  # note to self: three sites are close together, so maybe use another inset map??
  geom_point(
    aes(x = Longitude, y = Latitude, shape = Treatment), 
    data = site_tidy) + 
  labs(x = "Longitude", y = "Latitude") + 
  scale_shape_manual(
    name = "Management Regime", 
    labels = c("Maintence-Mowing", "Undisturbed", "Tilled"),
    values = c(0, 1, 2)
  ) + 
  geom_text_repel(
    aes(x = Longitude, y = Latitude, label = Site),
    size = 2, 
    data = site_tidy) + 
  
  
  scalebar(
    data = mw_shp_tidy, 
    dist = 1 ,
    transform = TRUE, 
    dist_unit = "km",
    st.size = 2) +
  
  
  theme_bw() + 
  theme(
    legend.title = element_text(size = 5),
    legend.text = element_text(size = 5)
  )

## inset map -------------------------------------------------------------------
inset_map <- cowplot::ggdraw() + 
  cowplot::draw_plot(mw) + 
  cowplot::draw_plot(to, x = 0.18, y = 0.58, width = 0.2, height = 0.2)

# save to disk -----------------------------------------------------------------
ggsave(
  filename = here("output", "results", "figure-1.png"),
  plot = inset_map, 
  device = "png", 
  width = 5,
  height = 5, 
  units = "in"
)
