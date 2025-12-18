###################### Preparing Mosquito Alert Data ############################
# This script links environmental variables derived from processed ERA5 data
# with Culex spp. presence records reported through the Mosquito Alert citizen
# science platform.
#
# To run this script, users must first download and preprocess ERA5 data into
# .txt format using the script 'read_and_extract_grib_to_txt.R' available in
# this repository.
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
# NOTE: we always used as template de same grid from reanalysis-era5-land, considering
# the following areas: area: [75, -20, 30, 45]. For its download, there is an jupyter
# file specific for ERA5 download in this repository
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
# Note: the data was downloaded from Mosquito Alert Data Portal
culex = readxl::read_xlsx(paste0(loc.data, 'culex_2020_2023_EnricPou.csv')) %>%
  rename("lon" = "location__point__longitude",
         "lat" = "location__point__latitude",
         "where" = "event_environment") %>%
  mutate(date = as.Date(created_at)) %>% 
  dplyr::select(date, lon, lat, where) %>%
  st_as_sf(coords = c("lon", "lat"), crs=4326, remove=FALSE)

index_intersects <- st_intersects(culex, europe) 
index_intersects <- lengths(index_intersects) > 0

culex <- culex[index_intersects, ] # Only selecting Europe data 

culex_total <- culex

for (where in c("outdoors", "indoors", "all")){
  if (where == "outdoors"){
    culex <- culex_total %>%
      filter(where == "outdoors") 
  } 
  if (where == "indoors"){
    culex <- culex_total %>%
      filter(where == "indoors") 
  } 
  if (where == "all"){
    culex <- culex_total 
  }
  
  # Plotting raw data --> checking the locations plot
  culex_raw_reports <- ggplot(culex) +
    geom_point(aes(x = lon, y = lat, color = category), size = 0.05) +
    geom_sf(data = europe, aes(fill = "transparent"), fill = "transparent") +
    labs(x = "Longitude",
         y = "Latitude",
         color = "Category") +
    xlim(c(-22, 43)) +
    ylim(c(35, 70)) +
    theme_void(base_size = 12, base_family = "Helvetica")
  
  # We have to add a new random effect: the pixel_id!
  culex$pixel_id <- as.factor(raster::cellFromXY(tmp_raster, st_coordinates(culex)))
  length(unique(culex$pixel_id)) 
  sum(is.na(culex$pixel_id)) 
  culex <- culex %>% drop_na(pixel_id)
  # Now we have incorporated the new ID cells
  
  # Grouping reports by pixel ID
  culex <- culex %>%
    st_drop_geometry() %>%
    group_by(pixel_id, date) %>% 
    summarize(n_target_reps = n(), .groups ="drop") 
  
  culex <- merge(culex, template, by = "pixel_id", all.x = TRUE)
  culex$geometry <- NULL
  summary(culex)
  
  culex_summary <- culex %>%
    group_by(pixel_id, lon, lat) %>%
    summarise(,
              n = n(),
              sum_targets = sum(n_target_reps),
              av_targets = sum_targets/n,
              presence = sum_targets > 0,
              .groups ="drop")
  
  # Plotting aggregated data
  culex_agg_plot <- ggplot(culex_summary) +
    geom_tile(aes(x = lon, y = lat, fill = sum_targets)) +
    geom_sf(data = europe, fill = "transparent") +
    scale_fill_gradientn(
      name = "N. Reports",
      colours = c("#DAA521", "#653496"),
      na.value = "white"
    ) +
    xlim(c(-22, 43)) +
    ylim(c(35, 70)) +
    theme_void(base_size = 12, base_family = "Helvetica") + 
    theme(legend.margin = margin(6, 12, 6, 6))
  
  ggpubr::ggarrange(culex_raw_reports, culex_agg_plot, labels = c("a", "b"))
  
  ggsave(file = file.path(loc.figures, sprintf("Exploratory_analysis_reports/Culex_reports_raw_and_agg_data_%s.pdf", where)), 
         width = 25, height = 12, dpi = 300, units = "cm", device = cairo_pdf)  
  
  # Adding Countries --> random effect -------------------------------------------
  culex <- st_intersection(culex %>% 
                             drop_na(lon, lat) %>%
                             st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE), 
                           europe %>% dplyr::select(NUTS_ID)) %>%
    rename("country" = "NUTS_ID") 
  
  ggplot(culex) +
    geom_sf() +
    geom_sf(data = europe, fill = "transparent") +
    xlim(c(-22, 43)) +
    ylim(c(35, 70)) +
    theme_void(base_size = 12, base_family = "Helvetica") + 
    theme(legend.margin = margin(6, 12, 6, 6))
  
  summary(culex) 
  
  # Adding weather information ---------------------------------------------------
  # To calculate weather information we need several functions:
  # One to load the weather tables to calculate accumulative data:
  # NOTE: the user should be created a folder with the ERA5 processed data
  load_weather <- function(dt){
    nd_date =  lubridate::month(dt)
    
    if (nd_date == 1){
      date_file1 <- paste0(lubridate::year(dt) - 1, "-", 12)
      date_file2 <- paste0(lubridate::year(dt), "-", "01")
    } else {
      date_file1 <- paste0(lubridate::year(dt), "-", sprintf("%02d", nd_date - 1))
      date_file2 <- paste0(lubridate::year(dt), "-", sprintf("%02d", nd_date))
    }
    
    wt1 <- fread(paste0(loc.output, "ERA5_txt_day/ERA5_", date_file1, ".txt"))[,-1]
    wt2 <- fread(paste0(loc.output, "ERA5_txt_day/ERA5_", date_file2, ".txt"))[,-1]
    wt <- rbind(wt1, wt2)
    
    return(wt)
  }
  
  calculating_weather <- function(ex_wt_row, data_row, dt){
    
    for (lag in c(0, 7, 14, 21)) {
      # Compute mean for all variables except precipitation
      mean_vars <- c("dewpoint", "mean_temperature", "u_wind", "v_wind", 
                     "max_temperature", "min_temperature", "mean_relative_humidity", "wind_speed")
      
      for (var in mean_vars) {
        
        if(lag == 0){
          data_row[[var]] <- ex_wt_row %>%
            filter(date == dt) %>%
            dplyr::select(all_of(var)) %>%
            summarise(mean_value = mean(!!sym(var), na.rm = TRUE)) %>%
            pull(mean_value)
        } else {
          new_col_name <- paste0("l", lag, "_", var)  # Generate column name dynamically
          data_row[[new_col_name]] <- ex_wt_row %>%
            filter(date >= dt - days(lag) & date < dt) %>%
            dplyr::select(all_of(var)) %>%
            summarise(mean_value = mean(!!sym(var), na.rm = TRUE)) %>%
            pull(mean_value)
        }
      }
      
      # Compute sum for precipitation
      if (lag == 0){
        data_row[["precipitation"]] <- ex_wt_row %>%
          filter(date == dt) %>%
          dplyr::select(precipitation) %>%
          summarise(sum_value = sum(precipitation, na.rm = TRUE)) %>%
          pull(sum_value)
      } else {
        precip_col_name <- paste0("l", lag, "_precipitation")
        data_row[[precip_col_name]] <- ex_wt_row %>%
          filter(date >= dt - days(lag) & date < dt) %>%
          dplyr::select(precipitation) %>%
          summarise(sum_value = sum(precipitation, na.rm = TRUE)) %>%
          pull(sum_value)
      }
    }
    return(data_row)
  }
  
  cat("------------------------- Number of rows:", nrow(culex), "\n") # 8369
  years = c("2020", "2021", "2022", "2023")
  months = c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")
  
  weather_data <- data.frame()
  for (year in years){
    for (month in months){
      culex_sample <- culex %>% filter(year(date) == as.numeric(year) & month(date) == as.numeric(month))
      
      if (nrow(culex_sample) != 0){
        # Loading weather data
        dm <- culex_sample[1,]$date
        
        cat("Loading ", month, "-", year, "\n")
        ex_wt <- load_weather(dm)
        ex_wt[, date := as.Date(gsub("^d_", "", date), format = "%Y_%m_%d")]
        
        cat("Calculating ..... \n")
        weather_sample <- mclapply(1:nrow(culex_sample), function(i){
          
          data_row <- culex_sample[i, ]
          nd_date <- as.Date(data_row$date)
          st_date <- nd_date - lubridate::days(21)
          ids <- data_row$pixel_id
          
          ex_wt_row <- ex_wt[pixel_id == ids & date >= st_date & date <= nd_date]
          
          # Calculating accumulative data
          wth <- calculating_weather(ex_wt_row, data_row, nd_date)
          
          return(wth)
        }
        , mc.cores = 10)
        
        weather_sample <- do.call(rbind, weather_sample) 
        weather_data <- rbind(weather_data, weather_sample)
      } else {
        cat("Loading ", month, "-", year, ": no data\n")
      }
    }
  }
  
  sum(duplicated(weather_data))
  saveRDS(weather_data, file = file.path(loc.output, sprintf("culex_MA_presences_%s.rds", where)))
}

