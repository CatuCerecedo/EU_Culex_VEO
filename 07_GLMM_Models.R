############################ Background Model (GLMM) ############################
# This script implements the informed-pseudoabsences modeling framework for Mosquito
# Alert data using Generalized Linear Mixed Models (GLMMs).
#
# Models include temporal autocorrelation through an AR(1) structure to account
# for serial dependence in mosquito occurrence.
################################################################################

library(lme4)
library(terra)
library(sf)
library(MuMIn)
library(ggplot2)
library(glmmTMB)
library(caret)
library(parallel)
library(DHARMa)
library(tidyverse)
library(pROC)

rm(list = ls())
sf::sf_use_s2(FALSE)
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

# Functions --------------------------------------------------------------------
# Marginal effects
marginaleffects <- function(model, train, var = "mean_temperature"){
  temp_seq <- seq(min(train[[var]], na.rm = TRUE),
                  max(train[[var]], na.rm = TRUE),
                  length.out = 100)
  
  newdata <- data.frame(
    mean_relative_humidity = rep(mean(train$mean_relative_humidity, na.rm = TRUE), 100),
    l21_precipitation = rep(mean(train$l21_precipitation, na.rm = TRUE), 100),
    agricultural = rep(mean(train$agricultural, na.rm = TRUE), 100),
    discont_urban_fabric = rep(mean(train$discont_urban_fabric, na.rm = TRUE), 100),
    forests_scrub = rep(mean(train$forests_scrub, na.rm = TRUE), 100),
    cont_urban_fabric = rep(mean(train$cont_urban_fabric, na.rm = TRUE), 100),
    SE = 1,
    area_km2 = rep(mean(train$area_km2, na.rm = TRUE), 100),
    country = NA,  # Para efectos aleatorios
    pixel_id = NA,
    wt = 1,
    year = NA,
    month = NA
  )
  
  newdata[[var]] <- temp_seq
  
  newdata$pred <- predict(model, newdata = newdata, type = "response", allow.new.levels = TRUE)
  
  return(newdata)
}

# We need the a dataframe with pixel ID
# Getting a template
my_grib = paste0(loc.output, "/ERA5_Download/ERA5_EU_monthly_")

year = "2023"
month = "01"

grib_file <- paste0(my_grib, year, "-", month, ".grib")
tmp_raster <- rast(grib_file)

# GLMM modeling ----------------------------------------------------------------
# Adding country_surface
europe <- eurostat::get_eurostat_geospatial(resolution = 10,
                                            nuts_level = 0,
                                            year = 2021) %>%
  sf::st_transform(4326) %>%
  sf::st_transform(3035) %>% 
  sf::st_as_sf() %>%
  mutate(area_km2 = as.numeric(units::set_units(sf::st_area(.), km^2))) %>%
  rename("country" = "NUTS_ID") %>%
  dplyr::select(country, area_km2) %>%
  sf::st_drop_geometry()

## Informative pseudo-absences --------------------------------------------------
# MA with clean data (obtained by Mosquito Alert map) 

for (where in c("outdoors", "indoors", "all")){
  ma <- readRDS(file.path(loc.output, sprintf("culex_MA_clean_%s_rm.rds", where)))  %>%
    mutate(year = as.factor(lubridate::year(date)),
           month = 12L * (lubridate::year(date) - min(lubridate::year(date), na.rm = TRUE)) + lubridate::month(date)
    ) %>%
    drop_na(l21_mean_temperature)
  ma <- merge(ma, europe, by = "country") 
  
  ## Taking RM in consideration to avoid less "reliable" zeros -------------------
  ma2 <- rbind(
    ma %>% filter(presence == FALSE),
    ma %>% filter(presence == TRUE & R0_cul_temp >= 1)
  )
  
  ma2 <- ma2 %>%
    arrange(pixel_id, month) %>%
    mutate(presence = as.logical(presence),
           month = as.factor(month))
  
  # Using Pseudo-absences > Q3
  q3 <- ma2 %>% filter(presence == FALSE)
  
  # Balancing number of zeros
  q3 <- q3 %>% filter(SE > quantile(q3$SE, 0.75))
  
  q3 <- q3 %>% sample_n(size = 2*(ma2 %>% filter(presence == TRUE)) %>% nrow())
  ma_q3 <- rbind(ma2 %>% filter(presence == TRUE),
                 q3)
  
  summary(ma_q3$presence)
  
  set.seed(126831623)
  sample <- sample(c(TRUE, FALSE), nrow(ma_q3), replace=TRUE, prob=c(0.8, 0.2))
  train  <- ma_q3[sample, ]
  test   <- ma_q3[!sample, ]
  
  glmm_q3_se <- glmmTMB(presence ~ poly(l21_mean_temperature, 2) +
                          mean_relative_humidity +
                          discont_urban_fabric + forests_scrub +
                          offset(log(SE)) + 
                          ar1(month + 0 | country)  + # WARNING: ranked month by group
                          (1 | country),
                        data = train,
                        family = binomial("logit"))
  summary(glmm_q3_se)
  
  newdata <- marginaleffects(glmm_q3_se, train, var = "l21_mean_temperature")
  ggplot(newdata, aes(x = l21_mean_temperature, y = pred)) +
    geom_line(color = "#653496", linewidth = 1) +
    labs(y = "Predicted ocurrance probability",
         x = "Mean Temperature (lag 21)") +
    theme_classic(base_size = 14, base_family = "Helvetica")
  
  
  # Model assessment
  simulationOutput <- simulateResiduals(fittedModel = glmm_q3_se)
  
  testUniformity(simulationOutput)
  testDispersion(simulationOutput)
  testOutliers(simulationOutput)
  
  site_id <- interaction(train$lon, train$lat, drop = TRUE)
  
  resid_agg <- recalculateResiduals(simulationOutput, group = site_id)
  lon_site <- tapply(train$lon, site_id, FUN = function(z) z[1])
  lat_site <- tapply(train$lat, site_id, FUN = function(z) z[1])
  testSpatialAutocorrelation(resid_agg,
                             x = as.numeric(lon_site),
                             y = as.numeric(lat_site))
  
  time_id <- as.Date(train$date)
  resid_agg <- recalculateResiduals(simulationOutput, group = time_id)
  time_unique <- tapply(time_id, interaction(time_id, drop = TRUE), `[`, 1)
  time_unique <- as.Date(time_unique) 
  testTemporalAutocorrelation(resid_agg, time = time_unique)
  
  # Cross-validation
  pred_probs <- predict(glmm_q3_se, newdata = test, type = "response", allow.new.levels=TRUE)
  roc_obj <- roc(test$presence, pred_probs)
  plot(roc_obj)
  auc(roc_obj)
  
  saveRDS(glmm_q3_se, file = file.path(loc.output, sprintf( "glmm_q3_2PA_ar1_country_truese_%s.rds", where)))
}

