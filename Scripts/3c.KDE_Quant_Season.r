#
# Seasonal KDE Analysis of Baleen Whale PAM Detections with Quantile Binning
# Updated 2025
#
# This script performs kernel density estimation (KDE) on passive acoustic monitoring (PAM)
# data for baleen whale species, stratified by season (Winter, Spring, Summer, Fall).
#
# Key features:
# - Aggregates daily PAM detections by site, species, and season
# - Calculates proportion of days with detections per season at each station
# - Generates weighted KDE surfaces using detection proportions as weights
# - Bins continuous KDE values into quantiles (50%, 80%, 90%, 95%)
# - Converts quantile-binned rasters to vector polygon shapefiles
# - Creates seasonal comparison plots (2x2 grid) for each species
# - Normalizes density values within each species for visual comparison
# - Uses common spatial window across all species/seasons for consistency
# - Outputs: continuous rasters (.tif), quantile shapefiles (.shp), and plots (.png)
#
# Species analyzed: Sei (Bb), Blue (Bm), Fin (Bp), Humpback (Mn), Minke (Ba), 
#                   North Atlantic Right Whale (Eg)
#
# Seasons defined: Winter (Dec-Feb), Spring (Mar-May), Summer (Jun-Aug), Fall (Sep-Nov)
#

options(scipen = 999)

pacman::p_load(sf, tidyverse, terra, spatstat, lubridate, patchwork, viridis, classInt, stars)

#========================
# INPUT / OUTPUT SETUP
#========================

filepath      <- "input/2025/baleen_presence_laura_2025.csv"
raster_dir    <- "output/tif/seasonal/"
shapefile_dir <- "output/shapes/seasonal/"
plot_dir      <- "output/FIGS/seasonal/"

if (!dir.exists(raster_dir)) dir.create(raster_dir, recursive = TRUE)
if (!dir.exists(shapefile_dir)) dir.create(shapefile_dir, recursive = TRUE)
if (!dir.exists(plot_dir))   dir.create(plot_dir, recursive = TRUE)

UTM20 <- "EPSG:32620"

species_list_pam <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")
species_names <- list(
  "Bb" = "Sei Whale",
  "Bm" = "Blue Whale",
  "Bp" = "Fin Whale",
  "Mn" = "Humpback Whale",
  "Ba" = "Minke Whale",
  "Eg" = "North Atlantic Right Whale"
)

season_levels <- c("Winter", "Spring", "Summer", "Fall")

#========================
# BUILD SEASONAL DATA
#========================

baleen_raw <- read.csv(filepath)

baleen_season <- baleen_raw %>%
  mutate(
    rec_date = as.Date(rec_date),
    month    = month(rec_date),
    Season   = case_when(
      month %in% c(12, 1, 2)  ~ "Winter",
      month %in% c(3, 4, 5)   ~ "Spring",
      month %in% c(6, 7, 8)   ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall",
      TRUE ~ NA_character_
    ),
    Season = factor(Season, levels = season_levels)
  ) %>%
  filter(!is.na(Season)) %>%
  group_by(site, latitude, longitude, species, Season) %>%
  summarise(
    total_days     = n_distinct(rec_date),
    detection_days = sum(presence > 0, na.rm = TRUE),
    proportion_det = detection_days / total_days,
    .groups = "drop"
  ) %>%
  filter(total_days > 0)

# sf version in UTM
baleen_season_sf <- baleen_season %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

shapefile_crs <- st_crs(baleen_season_sf)$wkt

#========================
# COMMON WINDOW
#========================

pam_coords <- st_coordinates(ess_study_final)
x_range <- range(pam_coords[, "X"])
y_range <- range(pam_coords[, "Y"])

buffer_percent <- 0.2
x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
y_buffer <- (y_range[2] - y_range[1]) * buffer_percent

expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

common_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)

#========================
# QUANTILE BINNING FUNCTION
#========================

convertRasterToQuantiles <- function(kde_raster, 
                                     species_code, 
                                     season_name,
                                     output_dir = shapefile_dir,
                                     probs = c(0.50, 0.80, 0.90, 0.95, 1)) {
  
  # Extract CRS
  raster_crs <- crs(kde_raster)
  
  # Set values <= 0 to NA
  kde_raster[kde_raster <= 0] <- NA
  
  # Calculate quantile breaks
  kde_values <- values(kde_raster, mat = FALSE)
  quantiles <- quantile(kde_values, probs = probs, type = 6, na.rm = TRUE)
  
  # Create labels
  quant_labels <- paste0(probs[-length(probs)] * 100, "%")
  
  # Bin values
  kde_binned_values <- cut(kde_values, 
                           breaks = quantiles, 
                           labels = quant_labels,
                           include.lowest = TRUE)
  
  # Create binned raster
  kde_quant_rast <- kde_raster
  values(kde_quant_rast) <- as.numeric(kde_binned_values)
  
  # Convert to polygons
  kde_polys <- as.polygons(kde_quant_rast, dissolve = TRUE)
  kde_sf <- st_as_sf(kde_polys)
  kde_sf <- st_make_valid(kde_sf)
  kde_sf <- st_set_crs(kde_sf, raster_crs)
  
  # Add quantile labels
  names(kde_sf)[1] <- "quant_level"
  kde_sf <- kde_sf %>%
    mutate(Quantile = quant_labels[quant_level],
           species = species_code,
           season = season_name)
  
  # Save shapefile
  sanitize_name <- function(x) gsub("[/\\:*?\"<>|]", "_", gsub(" ", "", trimws(x)))
  output_path <- file.path(output_dir, 
                           paste0("KDE_pam_", sanitize_name(species_code), 
                                  "_season_", sanitize_name(season_name), ".shp"))
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  cat("    Created quantile shapefile for", season_name, "\n")
  
  return(invisible(kde_sf))
}

#========================
# SEASONAL KDE FUNCTION WITH QUANTILES
#========================

performSeasonalKDE <- function(data_sf,
                               species_col,
                               season_col,
                               species_list,
                               species_names = NULL,
                               weight_col = NULL,
                               sigma_val = 10000,
                               window = NULL,
                               output_prefix = "pam",
                               output_dir = raster_dir,
                               shapefile_output_dir = shapefile_dir,
                               plot_output_dir = plot_dir,
                               create_quantiles = TRUE) {
  
  # helper for filenames
  sanitize_filename <- function(name) {
    name <- trimws(name)
    name <- gsub("[/\\:*?\"<>|]", "_", name)
    name <- gsub(" ", "", name)
    name
  }
  
  seasons <- sort(na.omit(unique(data_sf[[season_col]])))
  
  for (sp in species_list) {
    
    sp_name <- if (!is.null(species_names[[sp]])) species_names[[sp]] else sp
    cat("\n=== Species:", sp_name, "(", sp, ") ===\n")
    
    species_data <- data_sf[data_sf[[species_col]] == sp, ]
    if (nrow(species_data) == 0) {
      cat("  No data for", sp, "\n")
      next
    }
    
    seasonal_rasters <- list()
    seasonal_plots   <- list()
    global_min <- Inf
    global_max <- -Inf
    
    #---------------------
    # FIRST PASS: KDE BY SEASON
    #---------------------
    
    for (ss in seasons) {
      season_data <- species_data[species_data[[season_col]] == ss, ]
      if (nrow(season_data) == 0) {
        cat("  Season", ss, ": no data\n")
        next
      }
      
      cat("  Season", ss, ":", nrow(season_data), "stations\n")
      
      coords <- st_coordinates(season_data)
      
      # window
      win <- if (is.null(window)) {
        xr <- range(coords[, "X"])
        yr <- range(coords[, "Y"])
        xb <- (xr[2] - xr[1]) * 0.2
        yb <- (yr[2] - yr[1]) * 0.2
        owin(xrange = c(xr[1] - xb, xr[2] + xb),
             yrange = c(yr[1] - yb, yr[2] + yb))
      } else {
        window
      }
      
      ppp_obj <- ppp(x = coords[, "X"], y = coords[, "Y"], window = win)
      ppp_obj <- rjitter(ppp_obj, retry = TRUE, nsim = 1, drop = TRUE)
      
      if (!is.null(weight_col)) {
        wts <- season_data[[weight_col]] * 10000
        marks(ppp_obj) <- data.frame(weights = wts)
        w_arg <- marks(ppp_obj)
      } else {
        w_arg <- NULL
      }
      
      kd <- density.ppp(
        ppp_obj,
        sigma   = sigma_val,
        positive = TRUE,
        kernel  = "gaussian",
        weights = w_arg,
        dimyx   = 500,
        diggle  = TRUE
      )
      
      r_kde <- rast(kd)
      crs(r_kde) <- shapefile_crs
      
      max_val <- max(values(r_kde), na.rm = TRUE)
      if (!is.na(max_val) && max_val > 0) {
        r_kde <- r_kde / max_val
      }
      
      local_min <- min(values(r_kde), na.rm = TRUE)
      local_max <- max(values(r_kde), na.rm = TRUE)
      global_min <- min(global_min, local_min)
      global_max <- max(global_max, local_max)
      
      # Save continuous raster
      out_rast <- file.path(
        output_dir,
        paste0("KDE_", output_prefix, "_", sanitize_filename(sp),
               "_season_", sanitize_filename(ss), ".tif")
      )
      writeRaster(r_kde, out_rast, overwrite = TRUE)
      
      # Create quantile shapefile if requested
      if (create_quantiles) {
        tryCatch({
          convertRasterToQuantiles(
            kde_raster = r_kde,
            species_code = sp,
            season_name = as.character(ss),
            output_dir = shapefile_output_dir
          )
        }, error = function(e) {
          cat("    Error creating quantiles for", ss, ":", e$message, "\n")
        })
      }
      
      seasonal_rasters[[as.character(ss)]] <- r_kde
    }
    
    # nothing for this species
    if (!length(seasonal_rasters)) {
      cat("  No seasons produced KDE for", sp, "\n")
      next
    }
    
    #---------------------
    # SECOND PASS: PLOTS
    #---------------------
    
    for (ss in names(seasonal_rasters)) {
      r_kde <- seasonal_rasters[[ss]]
      kde_df <- as.data.frame(r_kde, xy = TRUE) %>%
        filter(!is.na(lyr.1))
      
      p <- ggplot(kde_df, aes(x = x, y = y, fill = lyr.1)) +
        geom_tile() +
        scale_fill_viridis_c(limits = c(global_min, global_max),
                             option = "viridis") +
        coord_equal() +
        theme_minimal() +
        theme(
          axis.title = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank()
        ) +
        labs(
          title = paste0(ss),
          fill  = "Density"
        )
      
      seasonal_plots[[ss]] <- p
    }
    
    combined_plot <- wrap_plots(seasonal_plots, ncol = 2) +
      plot_annotation(
        title = paste("Seasonal KDE -", sp_name),
        theme = theme(plot.title = element_text(size = 16, face = "bold"))
      )
    
    out_png <- file.path(
      plot_dir,
      paste0("combined_seasonal_KDE_", output_prefix, "_", sanitize_filename(sp), ".png")
    )
    ggsave(out_png, combined_plot, width = 12, height = 10, dpi = 300, bg = "white")
    cat("  Saved:", out_png, "\n")
  }
  
  cat("\n=== Seasonal KDE complete for all species ===\n")
  if (create_quantiles) {
    cat("Continuous rasters saved to:", output_dir, "\n")
    cat("Quantile shapefiles saved to:", shapefile_output_dir, "\n")
  }
}

#========================
# RUN SEASONAL KDE WITH QUANTILES
#========================

performSeasonalKDE(
  data_sf              = baleen_season_sf,
  species_col          = "species",
  season_col           = "Season",
  species_list         = species_list_pam,
  species_names        = species_names,
  weight_col           = "proportion_det",
  sigma_val            = 10000,
  window               = common_window,
  output_prefix        = "pam",
  create_quantiles     = TRUE  # Set to FALSE to skip quantile creation
)

cat("\n=== ALL OUTPUTS COMPLETE ===\n")
cat("Continuous rasters:", raster_dir, "\n")
cat("Quantile shapefiles:", shapefile_dir, "\n")
cat("Plots:", plot_dir, "\n")