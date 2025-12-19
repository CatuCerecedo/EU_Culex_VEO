#################### Plotting Indoor and Outdoor Presences ######################
# This script generates descriptive plots of Culex presence records reported
# through the Mosquito Alert platform, stratified by indoor, outdoor and all
# observations.
################################################################################

library(tidyverse)
library(sf)
library(terra)
library(parallel)
library(data.table)
library(sf)
library(cowplot)
library(patchwork)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

sf::sf_use_s2(FALSE)
# Loading European map ---------------------------------------------------------
europe <- eurostat::get_eurostat_geospatial(resolution = 10,
                                            nuts_level = 0,
                                            year = 2021) %>%
  st_transform(4326)
plot(st_geometry(europe))

# Getting a template
my_grib = paste0(loc.output, "/ERA5_Download/ERA5_EU_monthly_")

# We need the a dataframe with pixel ID
year = "2023"
month = "01"

grib_file <- paste0(my_grib, year, "-", month, ".grib")
tmp_raster <- rast(grib_file)

template <- as.data.frame(tmp_raster[[1]], xy = TRUE)[,1:2]
colnames(template) <- c("lon", "lat")
template <- template %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(tmp_raster), remove = FALSE)
template$pixel_id <- as.factor(raster::cellFromXY(tmp_raster, st_coordinates(template)))

# Loading reports data --------------------------------------------------------
culex = read_csv(paste0(loc.data, 'culex_2020_2023.csv')) %>%
  rename("lon" = "location__point__longitude",
         "lat" = "location__point__latitude",
         "where" = "event_environment") %>%
  mutate(date = as.Date(created_at)) %>% 
  dplyr::select(date, lon, lat, where) %>%
  st_as_sf(coords = c("lon", "lat"), crs=4326, remove=FALSE)

index_intersects <- st_intersects(culex, europe) 
index_intersects <- lengths(index_intersects) > 0

culex <- culex[index_intersects, ] # Only selecting Europe data 

culex$where[is.na(culex$where)] <- "Others"
culex$where[culex$where == "vehicle"] <- "Others"
culex$where[culex$where == "outdoors"] <- "Outdoors"
culex$where[culex$where == "indoors"] <- "Indoors"

culex$where <- factor(
  culex$where,
  levels = c("Outdoors", "Indoors", "Others")
)

ggplot(culex) +
  geom_sf(data = europe, fill = "transparent", color = "grey70", linewidth = 0.2) +
  geom_point(aes(x = lon, y = lat, color = where), size = 0.5, alpha = 0.7) +
  scale_colour_manual(values = c("Outdoors" = "#653496", "Indoors" = "#DAA521", "Others" = "#74B652")) +
  coord_sf(xlim = c(-22, 43), ylim = c(30, 73), expand = FALSE) +
  facet_wrap(~ where) +
  theme_void(base_size = 16, base_family = "Helvetica") +
  theme(legend.position = "bottom",
        plot.margin = unit(c(-5, 0, -5, 0), "cm"))

ggsave(file = paste0(loc.fig, "Exploratory_analysis_reports/Culex_reports_data.pdf"), 
       width = 25, height = 25, dpi = 300, units = "cm", device = cairo_pdf)

culex %>% group_by(where) %>% summarise(n = n())
