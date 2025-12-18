################# Preparing Mosquito Alert Data: joining #######################
library(tidyverse)
library(parallel)
library(data.table)
library(cowplot)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

sf::sf_use_s2(FALSE)
# Loading data -----------------------------------------------------------------

for (where in c("outdoors", "indoors", "all")){
  # Presences
  presence <- readRDS(file = paste0(loc.output, "culex_MA_presences_", where,".rds")) %>%
    sf::st_drop_geometry()
  
  # Pseudo-absences
  absence <- readRDS( file = paste0(loc.output, "culex_MA_pseudoabsences.rds")) %>%
    sf::st_drop_geometry() %>% dplyr::select(-geometry.x, -geometry.y)
  
  ## checking the overlapping ----------------------------------------------------
  pres <- presence %>% 
    dplyr::select(pixel_id, date, n_target_reps, lon, lat, country) %>%
    mutate(
      presence = n_target_reps > 0
    ) %>%
    dplyr::select(-n_target_reps)
  abs <- absence %>% 
    mutate(
      presence = FALSE
    ) %>%
    dplyr::select(pixel_id, date, lon, lat, country, SE, n_total_reports) 
  
  
  # How many presences have SE?
  pres_se <- merge(pres, abs, by = c("pixel_id", "date", "lon", "lat", "country")) %>%
    mutate(presence = ifelse(is.na(presence) == TRUE, FALSE, presence))
  
  pres <- presence %>% 
    mutate(
      presence = n_target_reps > 0
    ) %>%
    dplyr::select(-n_target_reps) 
  abs <- absence %>% 
    dplyr::select(-SE_expected) 
  
  # We started with that
  culex = abs %>% 
    full_join(pres, by = names(abs)[c(-6, -43)], multiple = "all") %>% # merging without SE
    replace_na(list(presence = FALSE)) %>%
    mutate(SE = if_else(presence == TRUE, 1, SE)) %>%
    filter(SE>0) 
  
  culex <- culex %>% group_by(date, pixel_id) %>%
    filter(if (any(presence)) presence else TRUE) %>%
    ungroup()
  sum(duplicated(culex))
  summary(culex$presence)
  
  # Saving with the calculated SE
  culex_setrue <- merge(abs, pres_se, by = c("pixel_id", "date", "lon", "lat", "country", "SE", "n_total_reports"), all = TRUE)
  culex_setrue <- culex_setrue %>%
    replace_na(list(presence = FALSE)) %>%
    filter(SE>0) %>%
    drop_na(SE)
  sum(duplicated(culex_setrue))
  summary(culex_setrue$presence)
  
  ## Adding CORINE land COVER ----------------------------------------------------
  corine <- readRDS(file = paste0(loc.output, "clc_label_0.rds"))
  
  culex_setrue <- merge(culex_setrue, 
                        corine %>% dplyr::select(-lon, -lat),
                        by = "pixel_id",
                        all.x = TRUE) 
  summary(culex_setrue)
  
  summary(culex_setrue$presence)
  
  culex_setrue <- culex_setrue %>% drop_na()
  
  saveRDS(culex_setrue, file = file.path(loc.output, sprintf("culex_MA_clean_%s.rds", where)))
}

