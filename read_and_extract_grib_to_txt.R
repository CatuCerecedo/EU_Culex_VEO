################# Preparing ERA5 Data: From GRIB to TXT #########################
# This script extracts meteorological variables from ERA5 GRIB files and
# converts them into plain-text (.txt) format for downstream analyses.
#
# ERA5 GRIB files must be downloaded in advance using the
# 'ERA5_daily_download.ipynb' notebook and stored in the directory:
# getwd()/ERA5_Download
#
# To sum up,
# Inputs:
#   - ERA5 GRIB files in getwd()/ERA5_Download
# Outputs:
#   - Processed ERA5 variables in .txt format
################################################################################

library(terra)
library(sf)
library(tidyverse)
library(janitor)
library(data.table)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.era5 <-paste0(getwd(), "/ERA5_Download/")

sf::sf_use_s2(FALSE)

# Elements to consider ---------------------------------------------------------
my_grib = paste0(loc.output, "/ERA5_Download/ERA5_EU_monthly_")

years = c("2020", "2021", "2022", "2023")
months = c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

# We need the a dataframe with pixel ID
year = "2023"
month = "01"

grib_file <- paste0(my_grib, year, "-", month, ".grib")
tmp_raster <- rast(grib_file)

europe <- as.data.frame(tmp_raster[[1]], xy = TRUE)[,1:2]
colnames(europe) <- c("lon", "lat")
europe <- europe %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(tmp_raster), remove = FALSE)
europe$pixel_id <- as.factor(raster::cellFromXY(tmp_raster, st_coordinates(europe)))

# Function: from grib (European level) to txt (daily-ERA5grid) -----------------
# The link og ERA5 data source: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview
# The grib have been downloaded though pyhton and the CDS API (see: https://cds.climate.copernicus.eu/how-to-api)
grib_to_text <- function(my_grib, years, months){
  for (year in years) {
    year_data <- data.frame()
    
    for (month in months) {
      date_file <- paste0(year, "-", month)
      cat(date_file, "-------------\n")
      grib_file <- paste0(my_grib, date_file, ".grib")
      
      tmp_raster <- rast(grib_file)
      vrs <- unique(names(tmp_raster))
      europe <- europe %>% st_transform(st_crs(tmp_raster))
      
      month_data <- list()
      i = 0
      for (vr in vrs){
        i = i +1
        print(vr)
        
        if (vr != "SFC (Ground or water surface); Total precipitation [m]"){
          tmp_raster_noprep <- tmp_raster[[-grep("precipitation", names(tmp_raster), ignore.case = TRUE)]]
          
          tmp_raster_vr <- subset(tmp_raster_noprep, names(tmp_raster_noprep) ==  vr)
          tmp_raster_vr <- tapp(tmp_raster_vr, index = "days", fun = mean, na.rm = NA)
          
          df <- terra::extract(tmp_raster_vr, europe) %>%  
            clean_names() %>%
            as.data.table()
          df$pixel_id <- europe$pixel_id
          df$lon <- europe$lon
          df$lat <- europe$lat
          
          # Cheeking NAs, if there is NA increase the polygon in order to capture the coast information
          # That happens only on the coast due to the NA cells from the raster
          
          df <- df %>% 
            pivot_longer(
              cols = starts_with("d_"),
              names_to = "date",
              values_to = vr
            ) %>%
            dplyr::select(-id)
          
          month_data[[i]] <- df %>% setDT()
        } else {
          # Precipitation needs another treatment
          tmp_raster_prep <- tmp_raster[[grep("precipitation", names(tmp_raster), ignore.case = TRUE)]]
          tmp_raster_vr <- tapp(tmp_raster_prep, index = "days", fun = sum, na.rm = TRUE)
          
          df <- terra::extract(tmp_raster_vr, europe) %>%  
            clean_names() %>%
            as.data.table()
          df$pixel_id <- europe$pixel_id
          df$lon <- europe$lon
          df$lat <- europe$lat
          
          df <- df %>% 
            pivot_longer(
              cols = starts_with("d_"),
              names_to = "date",
              values_to = vr
            ) %>%
            dplyr::select(-id)
          
          month_data[[i]] <- df %>% setDT()
        }
      }
      
      # We also need to calculate min and max temperature
      # Max temperature
      i = i+1
      print("Max temperature")
      
      tmp_raster_max <- tmp_raster[[grep("metre.temperature", names(tmp_raster), ignore.case = TRUE)]]
      tmp_raster_vr <- tapp(tmp_raster_max, index = "days", fun = max, na.rm = TRUE)
      
      df <- terra::extract(tmp_raster_vr, europe) %>%  
        clean_names() %>%
        as.data.table()
      df$pixel_id <- europe$pixel_id
      df$lon <- europe$lon
      df$lat <- europe$lat
      
      df <- df %>% 
        pivot_longer(
          cols = starts_with("d_"),
          names_to = "date",
          values_to = "max_temperature"
        ) %>%
        dplyr::select(-id)
      
      month_data[[i]] <- df %>% setDT()
      
      # Min temperature
      i = i+1
      print("Min temperature")
      
      tmp_raster_min <- tmp_raster[[grep("metre.temperature", names(tmp_raster), ignore.case = TRUE)]]
      tmp_raster_vr <- tapp(tmp_raster_min, index = "days", fun = min, na.rm = TRUE)
      
      df <- terra::extract(tmp_raster_vr, europe) %>%  
        clean_names() %>%
        as.data.table()
      df$pixel_id <- europe$pixel_id
      df$lon <- europe$lon
      df$lat <- europe$lat
      
      df <- df %>% 
        pivot_longer(
          cols = starts_with("d_"),
          names_to = "date",
          values_to = "min_temperature"
        ) %>%
        dplyr::select(-id)
      
      month_data[[i]] <- df %>% setDT()
      
      # Compeling all in the same data.frame
      month_data <- Reduce(function(x, y) merge(x, y, by = c("pixel_id", "lon", "lat", "date"), all = TRUE), month_data)
      
      # We need change variables names
      rename_columns <- function(col_name) {
        case_when(
          str_detect(col_name, "metre temperature") ~ "mean_temperature",
          str_detect(col_name, "dewpoint") ~ "dewpoint",
          str_detect(col_name, "u wind") ~ "u_wind",
          str_detect(col_name, "v wind") ~ "v_wind",
          str_detect(col_name, "Total precipitation") ~ "precipitation",
          TRUE ~ col_name  
        )
      }
      
      # Calculate the variables left
      month_data <-  month_data %>% 
        rename_with(rename_columns) %>%
        mutate(
          mean_temperature = mean_temperature - 273.15,
          max_temperature = max_temperature - 273.15,
          min_temperature = min_temperature - 273.15,
          dewpoint = dewpoint - 273.15
        ) %>%
        mutate(
          mean_relative_humidity = 100*10^(7.59138*((dewpoint/(dewpoint + 240.726))  - (mean_temperature/(mean_temperature + 240.726)) )),
          wind_speed = sqrt(u_wind^2 + v_wind^2)
        ) 
      
      cat("Good: saving month\n")
      write.table(month_data, paste0(loc.output, "ERA5_txt_day/ERA5_", date_file, ".txt"))
    }
  }
}

years = c("2019", "2020", "2021", "2022", "2023")
months = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")

grib_to_text(my_grib, years, months)

