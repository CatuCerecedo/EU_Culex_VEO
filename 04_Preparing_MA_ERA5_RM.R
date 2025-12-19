################################# RM filter ####################################
# This script merges the cleaned Mosquito Alert dataset with the RM filter.
#
# The RM filter is an epidemiological screening criterion used to reduce the
# influence of likely misidentifications in Mosquito Alert reports.
# 
# RM is calculated monthly using the methodology described by Pardo-Araujo et al.
# 2024. Proceedings B.291: 20241960
################################################################################

library(terra)
library(sf)
library(ggplot2)
library(parallel)
library(tidyverse)

rm(list = ls())
sf::sf_use_s2(FALSE)
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.fig <- paste0(getwd(), "/FIGURES/")

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
  ma <- readRDS(file.path(loc.output, sprintf("culex_MA_random_clean_%s.rds", where)))
  ma <- merge(ma, europe, by = "country") 
  
  ## Taking RM in consideration to avoid less "reliable" zeros -------------------
  rm_adding <- function(ma){
    ma_rm <- data.frame()
    for (year in c("2020", "2021", "2022", "2023")) {
      for (month in c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12")){
        cat(year, "-", month, "--------\n")
        rm <- readRDS(paste0(loc.data, "RM/RM_culex_temp_", year, "-", month, ".Rds")) %>%
          st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE)
        rm$pixel_id <- as.factor(raster::cellFromXY(tmp_raster, st_coordinates(rm)))
        rm <- rm[c("date", "R0_cul_temp", "pixel_id")] %>% st_drop_geometry()
        rm$date <- as.Date(rm$date, format = c("%Y.%m.%d"))
        
        ma_sample <- merge(ma, rm, by = c("date", "pixel_id"), all.x = TRUE)
        ma_sample <- ma_sample %>% drop_na(R0_cul_temp)
        
        if (nrow(ma_sample) > 0){
          ma_rm <- rbind(ma_rm, ma_sample)
        } else {
          cat("No data in this month\n")
        }
      }
    }
    return(ma_rm)
  }
  
  ma_rm_full <- rm_adding(ma)
  
  saveRDS(ma_rm_full, file.path(loc.output, sprintf("culex_MA_random_clean_%s_rm.rds", where))) 
}

# Plots ------------------------------------------------------------------------
# Outdoors
ma_out <- readRDS(paste0(loc.output, "culex_MA_clean_outdoors_rm2.rds")) %>%
  mutate(source = "Outdoors",
         filter_status = "Unfiltered")

ma_out_filt <- ma_out %>%
  filter((presence == FALSE & R0_cul_temp < 1) |
           (presence == TRUE  & R0_cul_temp >= 1)) %>%
  mutate(filter_status = "Filtered")

# Indoors
ma_in      <- readRDS(paste0(loc.output, "culex_MA_clean_indoors_rm.rds")) %>%
  mutate(source = "Indoors",
         filter_status = "Unfiltered")

ma_in_filt <- ma_in %>%
  filter((presence == FALSE & R0_cul_temp < 1) |
           (presence == TRUE  & R0_cul_temp >= 1)) %>%
  mutate(filter_status = "Filtered")

# All
ma_all      <- readRDS(paste0(loc.output, "culex_MA_clean_rm.rds")) %>%
  mutate(source = "All",
         filter_status = "Unfiltered")

ma_all_filt <- ma_all %>%
  filter((presence == FALSE & R0_cul_temp < 1) |
           (presence == TRUE  & R0_cul_temp >= 1)) %>%
  mutate(filter_status = "Filtered")

df_plot <- bind_rows(ma_out, ma_out_filt,
                     ma_in, ma_in_filt,
                     ma_all, ma_all_filt
) %>%
  mutate(
    filter_status = factor(filter_status, levels = c("Unfiltered", "Filtered")),
    source = factor(source, levels = c("Outdoors", "Indoors", "All")),
  )

Temp_plot <- ggplot(df_plot %>% filter(presence == TRUE), 
                    aes(x = filter_status, y = mean_temperature, fill = filter_status)) +
  annotate("rect",
           xmin = -Inf, xmax = Inf,
           ymin = 22, ymax = 23,
           alpha = 0.8,
           fill  = "#EAE2D6") +
  
  geom_violin(alpha = 0.4, color = NA, trim = FALSE) +
  geom_boxplot(width = 0.2, alpha = 0.9) +
  geom_hline(yintercept = 8.8,  linetype = 2, linewidth = 0.3) +
  geom_hline(yintercept = 34,  linetype = 2, linewidth = 0.3) +
  scale_fill_manual(values = c("Unfiltered" = "#653496",
                               "Filtered"  = "#DAA521"),
                    labels = c("Unfiltered", "Filtered"),
                    name   = "Culex pipiens\npresences") +
  scale_x_discrete(labels = c("Unfiltered" = "Unfiltered",
                              "Filtered"  = "Filtered")) +
  
  labs(x = NULL,
       y = "Temperature (ºC)") +
  
  facet_grid(~ source) +  # FILAS = Outdoors/Indoors/All; COLUMNAS = Raw/Filtered
  theme_classic(base_size = 16, base_family = "Helvetica") 
Temp_plot

combo <- readRDS(paste0(loc.output, "marginal_effects_plot_all.rds"))
combo

ggpubr::ggarrange(Temp_plot, combo, nrow = 2, labels = c("a", ""))

ggsave(file = paste0(loc.fig, "RM_filter/filter_presence_combinations_temperature.pdf"), 
       width = 25, height = 25, dpi = 300, units = "cm", device = cairo_pdf)
