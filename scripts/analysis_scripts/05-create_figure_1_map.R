# library ----------------------------------------------------------------------
library(ggplot2)
library(ggsn)
library(sf)
library(here)
library(dplyr)
library(gghighlight)
library(ggmap)
library(opendatatoronto)
library(cowplot)

# import data ------------------------------------------------------------------

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

mw_shp_tidy <- mw_shp %>%
  mutate(section = factor(section)) %>%
  st_transform(crs = 4326)

bound_tidy <- bound %>%
  st_transform(crs = 4326)

# visualize data ---------------------------------------------------------------

# toronto boundary
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
  
# meadow
mw <- ggplot() + 
  geom_sf(aes(fill = section), data = mw_shp_tidy) + 
  gghighlight(section %in% c("2", "4")) +
  labs(x = "Longitude", y = "Latitude") + 
  scalebar(
    data = mw_shp_tidy, 
    dist = 1 ,
    transform = TRUE, 
    dist_unit = "km",
    st.size = 2) +
  theme_bw()

# inset map
inset_map <- cowplot::ggdraw() + 
  cowplot::draw_plot(mw) + 
  cowplot::draw_plot(to, x = 0.18, y = 0.60, width = 0.3, height = 0.3)

# save to disk -----------------------------------------------------------------
