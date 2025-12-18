################ Preparing Mosquito Alert Data: Pseudo-absences #################
# This script downloads and processes sampling effort (SE) information from the
# Mosquito Alert platform, which is used as a proxy to define reliable
# pseudo-absence locations in presence-only modeling.
#
# ERA5 environmental variables must be preprocessed in advance into .txt format
# using the script 'read_and_extract_grib_to_txt.R' available in this repository.
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

# Raster template
pname <- "ERA5_EU_monthly_"

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

# # Checking point cells
# ggplot(template) + geom_point(aes(x = lon, y = lat), size = 0.001)

# Download Sampling effort data -----------------------------------------------
trs_daily <- read_csv("https://github.com/Mosquito-Alert/sampling_effort_data/raw/main/sampling_effort_daily_cellres_025.csv.gz") %>%
  as_tibble() %>% 
  mutate(
    date = as_date(date),
    year = year(as_date(date)),
    month = month(as_date(date)),
    n_total_reports = n_reports_albopictus + n_reports_bite + n_reports_culex + n_reports_japonicus + n_reports_aegypti + n_reports_koreicus,
    n_total_reporters = n_reporters_albopictus + n_reporters_bite + n_reporters_culex + n_reporters_japonicus + n_reporters_aegypti + n_reporters_koreicus
  ) %>%
  filter(year %in% c("2020", "2021", "2022", "2023", )) %>%
  st_as_sf(coords = c("masked_lon", "masked_lat"), remove = FALSE, crs = 4326)

index_intersects <- st_intersects(trs_daily, europe) 
index_intersects <- lengths(index_intersects) > 0

trs_daily <- trs_daily[index_intersects, ] # Only selecting Europe data

# We have to add a new random effect: the pixel_id!
trs_daily$pixel_id <- as.factor(raster::cellFromXY(tmp_raster, st_coordinates(trs_daily)))

# Now aggregate the data by pixel_id, lon and lat
se <- trs_daily %>%
  st_drop_geometry() %>%
  group_by(pixel_id, date) %>%
  summarise(
    SE = 1-prod(1-SE),
    SE_expected = sum(SE_expected),
    n_total_reports = sum(n_total_reports),
    .groups ="drop"
  ) 
summary(se)
cor(se$SE, se$SE_expected)
cor(se$SE, se$n_total_reports)
cor(se$SE_expected, se$n_total_reports)

se <- merge(se, template, by = "pixel_id", all.x = TRUE)
se$geometry <- NULL
summary(se)

# Adding Countries --> random effect -------------------------------------------
se <- st_intersection(se %>% 
                        drop_na(lon, lat) %>% 
                        st_as_sf(coords = c("lon", "lat"), crs = 4326, remove = FALSE), 
                     europe %>% 
                       dplyr::select(NUTS_ID)) %>%
  rename("country" = "NUTS_ID") 

ggplot(se) +
  geom_point(aes(x = lon, y = lat)) +
  geom_sf(data = europe, fill = "transparent") +
  xlim(c(-22, 43)) +
  ylim(c(35, 70)) +
  theme_void(base_size = 12, base_family = "Helvetica") + 
  theme(legend.margin = margin(6, 12, 6, 6))

summary(se) 

# Plotting historical SE on Europe
trs_summary <- setDT(se)
trs_summary <- trs_summary[, .(SE = 1-prod(1-SE),
                               SE_expected = sum(SE_expected),
                               n_total_reports = sum(n_total_reports)),
                           by = .(pixel_id, lon, lat)]
summary(trs_summary)

aSE <- ggplot(trs_summary) +
  geom_tile(aes(x = lon, y = lat, fill = SE)) +
  geom_sf(data = europe, aes(fill = "transparent"), fill = "transparent") +
  scale_fill_gradientn(
    name = "Sampling Effort\n(probability)",
    colours = c("#DAA521", "#653496"),
    na.value = "white"
  ) +
  xlim(c(-22, 43)) +
  ylim(c(35, 70)) +
  theme_void(base_size = 10, base_family = "Helvetica") +
  theme(plot.margin = margin(l = -1000, r = -1000))

dens_SE <- ggplot(trs_summary) + 
  geom_density(aes(x = SE), fill = "#653496") +
  labs(y = "Sampling Effort\n(probability)") +
  theme_classic(base_size = 10, base_family = "Helvetica") 

ggpubr::ggarrange(
  aSE,
  dens_SE, 
  ncol = 1,
  nrow = 2,
  heights = c(1.5, 0.6),
  widths = c(1, 0.2),
  labels = c("a", "b"),
  label.x = c(0, 0),  
  label.y = c(1, 1)
)

ggsave(file = paste0(loc.fig, "Exploratory_analysis_reports/Culex_se_raw_data_1.pdf"), 
       width = 25, height = 25, dpi = 300, units = "cm", device = cairo_pdf) 


# Adding weather information ---------------------------------------------------
# To calculate weather information we need several functions:
# One to load the weather tables to calculate accumalative data:
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

cat("------------------------- Number of rows:", nrow(se), "\n") # 173765
years = c("2020", "2021", "2022", "2023")
months = c("01", "02", "03", "04", "05", "06" ,"07", "08", "09", "10", "11", "12")

for (year in years){
  weather_data <- data.frame()
  
  for (month in months){
    se_sample <- se %>% filter(year(date) == as.numeric(year) & month(date) == as.numeric(month))
    
    if (nrow(se_sample) != 0){
      # Loading weather data
      dm <- se_sample[1,]$date
      
      cat("Loading ", month, "-", year, "\n")
      ex_wt <- load_weather(dm)
      ex_wt[, date := as.Date(gsub("^d_", "", date), format = "%Y_%m_%d")]
      
      cat("Calculating ..... \n")
      weather_sample <- mclapply(1:nrow(se_sample), function(i){
        
        data_row <- se_sample[i, ]
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
    
    print("Saving...")
    saveRDS(weather_data, file = paste0(loc.output, "culex_MA_pseudoabsences_", year, ".rds"))
  }
}

# Joining tables:
culex_MA_pseudoabsences <- data.frame()
for (year in years) {
  print(year)
  weather_data <- readRDS(file = paste0(loc.output, "culex_MA_pseudoabsences_", year, ".rds"))
  
  culex_MA_pseudoabsences <- rbind(culex_MA_pseudoabsences, weather_data)
}

culex_MA_pseudoabsences <- st_drop_geometry(culex_MA_pseudoabsences)
sum(duplicated(culex_MA_pseudoabsences))
saveRDS(culex_MA_pseudoabsences, file = paste0(loc.output, "culex_MA_pseudoabsences.rds"))
