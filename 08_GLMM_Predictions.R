######################### Prediction: era5, clc, model #########################
# This script generates spatial predictions from fitted models by combining
# ERA5-derived environmental variables and CORINE Land Cover (CLC) data.
#
# ERA5 data must be downloaded and preprocessed into .txt format using the
# script 'read_and_extract_grib_to_txt.R' available in this repository.
################################################################################

library(tidyverse)
library(sf)
library(parallel)
library(dplyr)
library(lubridate)
library(data.table)
library(terra)
library(glmmTMB)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")

sf::sf_use_s2(FALSE)
# Loading European map ---------------------------------------------------------
europe <- eurostat::get_eurostat_geospatial(resolution = 10,
                                            nuts_level = 0,
                                            year = 2021) %>%
  st_transform(4326)
plot(st_geometry(europe))

area <- eurostat::get_eurostat_geospatial(resolution = 10,
                                            nuts_level = 0,
                                            year = 2021) %>%
  sf::st_transform(4326) %>%
  sf::st_transform(3035) %>% 
  sf::st_as_sf() %>%
  mutate(area_km2 = as.numeric(units::set_units(sf::st_area(.), km^2))) %>%
  rename("country" = "NUTS_ID") %>%
  dplyr::select(country, area_km2) %>%
  sf::st_drop_geometry()

# Raster template
pname <- "ERA5_EU_monthly_"

# Getting a template
my_grib <- paste0("~/Documents/Mosquito_Models/ERA5_Download/EU_ERA5_Land_hourly_data/", pname)

# We need the a dataframe with pixel ID
year = "2023"
month = "01"

grib_file <- paste0(my_grib, year, "-", month, ".grib")
tmp_raster <- rast(grib_file)

template <- as.data.frame(tmp_raster[[1]], xy = TRUE)[,1:2]
colnames(template) <- c("lon", "lat")
template <- template %>% st_as_sf(coords = c("lon", "lat"), crs = st_crs(tmp_raster), remove = FALSE)
template$pixel_id <- as.factor(raster::cellFromXY(tmp_raster, st_coordinates(template)))

# Loading CLC data ------------------------------------------------------------- 
clc_surface <- readRDS(paste0(loc.output, "clc_label_0.rds")) %>%
  dplyr::select(-lon, -lat)

# Adding weather information ---------------------------------------------------
# To calculate weather information we need several functions:
# One to load the weather tables to calculate accumulative data:
load_weather <- function(dt){
  nd_date <- as.Date(dt)
  st_date <- nd_date - lubridate::days(21)
  
  months_needed <- unique(format(seq.Date(st_date, nd_date, by = "day"), "%Y-%m"))
  
  weather_list <- lapply(months_needed, function(m){
    file_path <- paste0(loc.output, "ERA5_txt_day/ERA5_", m, ".txt")
    fread(file_path)
  })
  
  wt <- rbindlist(weather_list)
  wt[, date := as.Date(date, format = "d_%Y_%m_%d")]
  
  wt <- wt[date >= st_date & date <= nd_date]
  return(wt)
}

calculate_weather_dt <- function(weather_dt, date_point){
  lags <- c(0, 21)
  mean_vars <- c("dewpoint", "mean_temperature", "u_wind", "v_wind", 
                 "max_temperature", "min_temperature", "mean_relative_humidity", "wind_speed")
  
  results <- list()
  
  for (lag in lags) {
    if (lag == 0) {
      subset_dt <- weather_dt[date == date_point]
    } else {
      subset_dt <- weather_dt[date >= date_point - days(lag) & date < date_point]
    }
    
    res <- subset_dt[, c(
      lapply(.SD[, ..mean_vars], mean, na.rm = TRUE),
      precipitation = sum(precipitation, na.rm = TRUE)
    ), by = pixel_id]
    
    # Renombrar columnas
    setnames(res, old = setdiff(names(res), "pixel_id"),
             new = if (lag == 0) {
               c(mean_vars, "precipitation")
             } else {
               c(paste0("l", lag, "_", mean_vars), paste0("l", lag, "_precipitation"))
             })
    
    results[[as.character(lag)]] <- res
  }
  
  # Fusionar todos los lags por ID
  result_final <- Reduce(function(x, y) merge(x, y, by = "pixel_id", all = TRUE), results)
  result_final <- result_final
  return(result_final)
}

# Predictions ------------------------------------------------------------------
ncores = 10

for (where in c("outdoors", "indoors", "all")){
  # Loading model
  model <- readRDS(file = file.path(loc.output, sprintf( "glmm_q3_2PA_ar1_country_truese_%s.rds", where)))
  
  mdl_name =  sprintf( "glmm_q3_2PA_ar1_country_truese_%s", where)
  
  # Preparing ERA5 folder
  for (year in c("2023")){
    
    # Months to predict
    months <- c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")
    
    for (m in months){
      
      # Empty box
      pred_points <- template %>% st_drop_geometry()
      
      print(paste("month: ", m, "year: ", year))
      
      # Calculate the first date of the month
      first_date <- lubridate::as_date(paste(year, m, "01", sep = "-"))
      
      # Calculate the last date of the month
      last_day <- ceiling_date(lubridate::as_date(paste(year, m, "01", sep = "-")), "month") - 1
      
      # Function to extract daily data -----------------------------------------------
      for (sel_date in seq.Date(first_date, last_day, 1)){
        sel_date <- as_date(sel_date)
        print(sel_date)
        
        cat("Loading and calulating weather variables -------\n")
        weather <- load_weather(sel_date)
        weather <- calculate_weather_dt(weather, sel_date)
        
        cat("Adding CLC variables -------\n")
        weather <- merge(weather %>% mutate(pixel_id = as.factor(pixel_id)), clc_surface, by = "pixel_id", all.x = TRUE)
        
        cat("------------------------------------------- Done\n")
        
        cat("Country assignation ----------------------------\n")
        data_prep_day <- merge(pred_points, weather, by = c("pixel_id"), all.x = TRUE)
        data_prep_day <- st_intersection(data_prep_day %>% 
                                           drop_na(lon, lat) %>% 
                                           st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE), 
                                         europe %>% 
                                           dplyr::select(NUTS_ID)) %>%
          rename("country" = "NUTS_ID") 
        cat("------------------------------------------- Done\n")
        
        # Function to predict daily data -----------------------------------------------
        data_prep_day <- data_prep_day %>% 
          mutate(
            SE = 0.1,
            n_traps = 1,
            trapping_effort = 7,
            year = year,
            month = as.factor(lubridate::month(sel_date))
          ) %>% 
          st_drop_geometry()
        data_prep_day <- merge(data_prep_day, area, by = "country")
        
        nrow_these_pred_points = nrow(data_prep_day)
        max_chunksize = 3000
        chunksize = min(as.integer((nrow_these_pred_points/ncores)), max_chunksize)
        
        pred <- bind_rows(mclapply(seq(1, nrow(data_prep_day[1:10]), chunksize), function(i){
          print(i)
          
          data_chunk = data_prep_day[i:min(nrow(data_prep_day), (i+(chunksize-1))), ]
          flush.console()
          
          pp <-  predict(model, newdata = data_chunk, 
                         allow.new.levels=TRUE, 
                         re.form = NA,
                         type = "response")
          
          pp <- data.frame(
            pixel_id = data_chunk$pixel_id,
            lon = data_chunk$lon,
            lat = data_chunk$lat,
            pred = pp
          )
          colnames(pp) <- c("pixel_id", "lon", "lat",  as.character(sel_date))
          
          return(pp)
          
        }, mc.cores = ncores))
        
        cat("------------------------------------- Prediction Done\n")
        
        pred_points <- merge(pred_points, pred,  by = c("pixel_id", "lon", "lat"))
        
        gc() # Clean space on R environment
        
        print(paste("Saving pred: ", m))
        saveRDS(pred_points, file = paste0(loc.output, "PREDICTIONS/GLMM/", mdl_name,
                                           "/", m, "_", year, ".rds"))
      } # Days
    } # months
  } # years
}

