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
library(stringr)

# import data ------------------------------------------------------------------

## source functions ------------------------------------------------------------

ddm_to_dd <- function(dms) {
  
  # split strings to seprate degrees and minutes.m
  split1 <- str_split(as.vector(dms), pattern = " ")
  
  # extract degrees
  degrees <- substr(split1[[1]][2], 1, 2)
  degrees <- as.double(degrees)
  
  # extract minutes.m
  minutes.m <- str_replace(split1[[1]][3], pattern = "'", replace = "")
  minutes.m <- as.double(minutes.m)
  
  # convert to degrees
  dot_d <- minutes.m/60
  
  # calculate decimal degrees
  decimal_degrees = degrees + dot_d
  
  
  return(decimal_degrees)
}

## plot coordinates ------------------------------------------------------------

bnsh <- read.csv(
  here("data", "input_data", "plot_coordinates", 
       "meadoway_site_BNSH_plot_coordinates.csv"), 
)

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

bsnh_tidy <- bnsh %>%
  mutate(
    Latitude  = ddm_to_dd(Latitude),
    Longitude = ddm_to_dd(Longitude)*-1
    )

## ensure consistent coordinate system  ----------------------------------------
bsnh_tidy <- st_as_sf(
  bsnh_tidy, 
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
  geom_sf(aes(fill = section), data = mw_shp_tidy) +
  gghighlight(section %in% c("2", "4")) +
  geom_sf(data = bsnh_tidy) + 
  labs(x = "Longitude", y = "Latitude") + 
  scalebar(
    data = mw_shp_tidy, 
    dist = 1 ,
    transform = TRUE, 
    dist_unit = "km",
    st.size = 2) +
  theme_bw()

## inset map -------------------------------------------------------------------
inset_map <- cowplot::ggdraw() + 
  cowplot::draw_plot(mw) + 
  cowplot::draw_plot(to, x = 0.2, y = 0.60, width = 0.3, height = 0.3)

# save to disk -----------------------------------------------------------------
