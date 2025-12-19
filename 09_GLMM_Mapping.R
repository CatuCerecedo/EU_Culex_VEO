########################### Mapping the predictions ############################
# This script generates spatial maps of model-based predictions.
################################################################################

library(dplyr)
library(sf)
library(ggplot2)
library(janitor)
library(scales)
library(purrr)
library(tidyr)
library(terra)

rm(list = ls())
# Directories ------------------------------------------------------------------
loc.output <- paste0(getwd(), "/OUTPUT/")
loc.data <- paste0(getwd(), "/DATA/")
loc.figures <- paste0(getwd(), "/FIGURES/")

sf::sf_use_s2(FALSE)
# Europe map -------------------------------------------------------------------
europe <- eurostat::get_eurostat_geospatial(resolution = 10,
                                            nuts_level = 0,
                                            year = 2021) %>%
  st_transform(4326)

# Load data --------------------------------------------------------------------
read_preds <- function(model_dir, label) {
  f <- readRDS(file.path(loc.output, "PREDICTIONS", "GLMM", model_dir, paste0(m, "_", year, ".rds")))
  f <- clean_names(f)
  # compute 'total' as mean across all numeric columns except lon/lat
  num_cols <- setdiff(names(f)[sapply(f, is.numeric)], c("pixel_id", "lon", "lat"))
  f$total  <- rowMeans(f[, num_cols], na.rm = TRUE)
  f$total  <- f$total/max(f$total, na.rm = TRUE)
  f$source <- label
  f <- f[c("pixel_id", "lon", "lat", "total", "source")]
  st_as_sf(f, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

# Plot by subset ---------------------------------------------------------------
year = "2023"

for (m in c("02", "04", "07", "11")){
  
  pred_outdoors <- read_preds("glmm_q3_2PA_ar1_country_truese_outdoors", "Outdoors")
  pred_indoors  <- read_preds("glmm_q3_2PA_ar1_country_truese_indoors",  "Indoors")
  pred_all      <- read_preds("glmm_q3_2PA_ar1_country_truese_all",      "All")
  
  preds <- dplyr::bind_rows(pred_outdoors, pred_indoors, pred_all) %>% 
    dplyr::mutate(
      source = factor(source, levels = c("Outdoors", "Indoors", "All")) 
    )
  
  p <- ggplot() +
    geom_tile(
      data = preds,
      aes(x = lon, y = lat, color = total), alpha = 0.8,
      na.rm = TRUE 
    ) +
    scale_color_gradientn(
      colours = c("#2E0C47",  
                  "#653496",  
                  "#2EC4B6",
                  "#B18910",  
                  "#DAA521"), 
      values  = scales::rescale(c(0, 0.2, 0.8, 1)),
      limits = c(0, 1),
      na.value = NA,                   
      name = "Predicted suitability",
      labels = scales::label_percent(accuracy = 1)
    ) +
    geom_sf(data = europe, fill = "transparent", color = "white", linewidth = 0.1) +
    coord_sf(xlim = c(-13, 45), ylim = c(34, 75), expand = FALSE) +
    facet_wrap(~ source, nrow = 1) +
    theme_void(base_size = 12, base_family = "Helvetica") +
    theme(
      strip.text      = element_text(face = "bold", size = 12),
      strip.background= element_rect(fill = "white", color = NA),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold"),
      plot.margin     = margin(6, 6, 6, 6)
    ) +
    guides(
      color = guide_colorbar(
        title.position = "top",
        barwidth = unit(10, "lines"),
        barheight = unit(0.6, "lines")
      )
    )
  
  p
  
  ggsave(
    file = file.path(loc.figures, sprintf("Predicted_maps/GLMM_predicted_suitability_%s_%s.pdf", m, year)),
    plot = p, width = 22, height = 16, units = "cm", dpi = 300, device = cairo_pdf
  )  
}

## By country ------------------------------------------------------------------
m = "07"

pred_outdoors <- read_preds("glmm_q3_2PA_ar1_country_truese_outdoors", "Outdoors")
pred_indoors  <- read_preds("glmm_q3_2PA_ar1_country_truese_indoors",  "Indoors")
pred_all      <- read_preds("glmm_q3_2PA_ar1_country_truese_all",      "All")

preds <- dplyr::bind_rows(pred_outdoors, pred_indoors, pred_all) %>% 
  dplyr::mutate(
    source = factor(source, levels = c("Outdoors", "All", "Indoors")) 
  )

preds <- st_intersection(preds, 
                         europe %>% dplyr::select(NUTS_ID)) %>%
  rename("country" = "NUTS_ID") %>%
  st_drop_geometry()

preds <- preds %>% filter(country %in% c("ES", "IT", "NL", "EL"))
preds <- preds %>% 
  mutate(country = case_when(
    country == "ES" ~ "Spain",
    country == "IT" ~ "Italy",
    country == "NL" ~ "Netherlands",
    country == "EL" ~ "Greece",
    .default = country))

p <- ggplot(preds, 
       aes(x = source, y = total, fill = source)) +
  geom_boxplot(width = 0.2, alpha = 0.9, outliers = FALSE) +
  scale_fill_manual(values = c( "Outdoors" = "#653496",
                                "Indoors"  = "#1F78B4",
                                "All"      = "#DAA521"),
                    name   = "Culex pipiens\ndata models") +
  labs(x = NULL,
       y = "Avrg predicted suitability") +
  facet_grid(~ country) + 
  theme_classic(base_size = 16, base_family = "Helvetica")  +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

p

ggsave(
  file = file.path(loc.figures, sprintf("Predicted_maps/GLMM_predicted_suitability_%s_%s_by_country.pdf", m, year)),
  plot = p, width = 22, height = 16, units = "cm", dpi = 300, device = cairo_pdf
)  

# Plot by month ----------------------------------------------------------------
read_preds_month <- function(m, year, mdl) {
  f <- readRDS(paste0(loc.output, "PREDICTIONS/GLMM/", mdl, "/", paste0(m, "_", year, ".rds")))
  f <- clean_names(f)
  num_cols <- setdiff(names(f)[sapply(f, is.numeric)], c("pixel_id", "lon", "lat"))
  f$total  <- rowMeans(f[, num_cols], na.rm = TRUE)
  f <- f[,c("pixel_id", "lon", "lat", "total")]
  # f$total  <- rank_order_conversion(f$total)
  f$total  <- f$total/max(f$total, na.rm = TRUE)
  f$month  <- m
  st_as_sf(f, coords = c("lon", "lat"), crs = 4326, remove = FALSE)
}

year <- "2023"
months <- c("02", "04", "07", "11")
mdl <- "glmm_q3_2PA_ar1_country_truese_outdoors"
mdl <- "glmm_q3_2PA_ar1_country_truese_indoors"
mdl <- "glmm_q3_2PA_ar1_country_truese_all"
mdl <- "glmm_ar1_country_truese_traps"

preds_out <- purrr::map_df(months, ~ read_preds_month(.x, year, mdl)) %>%
  mutate(month = factor(month, levels = months))
# preds_out$total  <- preds_out$total/max( preds_out$total, na.rm = TRUE)

p_out <- ggplot() +
  geom_tile(
    data = preds_out,
    aes(x = lon, y = lat, color = total), alpha = 0.2,
    na.rm = TRUE 
  ) +
  scale_color_gradientn(
    colours = c("#2E0C47",  
                "#653496",  
                "#2EC4B6",
                "#B18910",  
                "#DAA521"), 
    values  = scales::rescale(c(0, 0.2, 0.8, 1)),
    limits = c(0, 1),
    na.value = NA,                   
    name = "Predicted suitability",
    labels = scales::label_percent(accuracy = 1)
  ) +
  geom_sf(data = europe, fill = "transparent", color = "white", linewidth = 0.1) +
  coord_sf(xlim = c(-13, 45), ylim = c(34, 75), expand = FALSE) +
  facet_wrap(~ month, nrow = 1) +
  theme_void(base_size = 12, base_family = "Helvetica") +
  theme(
    strip.text       = element_text(face = "bold", size = 12),
    strip.background = element_rect(fill = "white", color = NA),
    legend.position  = "bottom",
    legend.title     = element_text(face = "bold"),
    plot.margin      = margin(6, 6, 6, 6)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      barwidth = unit(12, "lines"),
      barheight = unit(0.6, "lines")
    )
  )

p_out

ggsave(
  file = file.path(loc.figures, sprintf("Predicted_maps/GLMM_predicted_suitability_Traps_%s-%s.pdf",
                                    paste(months, collapse = "_"), year)),
  plot = p_out, width = 22, height = 16, units = "cm", dpi = 300, device = cairo_pdf
)

