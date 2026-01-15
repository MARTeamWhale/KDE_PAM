# Complete Monthly KDE analysis 
# This version creates cleaner, more comparable monthly maps
# Updated to clip data to study area

options(scipen = 999)

pacman::p_load(sf, tidyverse, readxl, here, leaflet, 
               scales, terra, ggrepel, viridis, ggspatial, spatstat, 
               dplyr, patchwork, lubridate)

# Directory setup
raster_dir <- "output/tif/monthly/"
plot_dir <- "output/FIGS/monthly/"

# Ensure directories exist
if (!dir.exists(raster_dir)) dir.create(raster_dir, recursive = TRUE)
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

# Projections
UTM20 <- crs("+init=epsg:32620")

#............................
# MONTHLY KDE FUNCTION (NO LOG TRANSFORM)------
#............................

performMonthlyKDE <- function(data_sf, 
                              species_col,
                              month_col,
                              species_list, 
                              species_names = NULL,
                              weight_col = NULL,
                              buffer_percent = 0.2, 
                              sigma_val = 10000,
                              window = NULL,
                              output_prefix,
                              output_dir = "output/tif/monthly/",
                              plot_output_dir = "output/FIGS/monthly/") {
  
  # Ensure output directories exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
  
  shapefile_crs <- st_crs(data_sf)$wkt
  
  # Sanitize filename helper
  sanitize_filename <- function(name) {
    name <- trimws(name)
    name <- gsub("[/\\:*?\"<>|]", "_", name)
    name <- gsub(" ", "", name)
    return(name)
  }
  
  # Month names for labeling
  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  
  #.............................
  # PROCESS EACH SPECIES------
  #.............................
  
  for(species in species_list) {
    # Get common name if provided
    species_display <- species
    if(!is.null(species_names) && species %in% names(species_names)) {
      species_display <- species_names[[species]]
    }
    
    cat("\n=== Processing species:", species_display, "(", species, ") ===\n")
    
    # Filter data for current species
    species_data <- data_sf[data_sf[[species_col]] == species, ]
    
    if(nrow(species_data) == 0) {
      cat("No data found for", species, "\n")
      next
    }
    
    # Initialize storage for this species
    monthly_rasters <- list()
    monthly_plots <- list()
    global_min_species <- Inf
    global_max_species <- -Inf
    
    #.............................
    # FIRST PASS: Create KDE for each month------
    #.............................
    
    for(month in 1:12) {
      month_data <- species_data[species_data[[month_col]] == month, ]
      
      if(nrow(month_data) == 0) {
        cat("  Month", month, "(", month_names[month], "): No data\n")
        next
      }
      
      cat("  Month", month, "(", month_names[month], "):", nrow(month_data), "stations\n")
      
      tryCatch({
        month_coords <- st_coordinates(month_data)
        
        # Create or use provided window
        if(is.null(window)) {
          x_range <- range(month_coords[, "X"])
          y_range <- range(month_coords[, "Y"])
          x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
          y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
          expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
          expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
          month_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
        } else {
          month_window <- window
        }
        
        # Create ppp object
        points_ppp <- ppp(x = month_coords[,"X"], 
                          y = month_coords[,"Y"], 
                          window = month_window)
        
        # Handle duplicates with jitter
        points_ppp <- rjitter(points_ppp, retry = TRUE, nsim = 1, drop = TRUE)
        
        # Apply weights if provided
        if(!is.null(weight_col)) {
          weights <- month_data[[weight_col]] * 10000
          marks(points_ppp) <- data.frame(weights)
          weight_arg <- marks(points_ppp)
        } else {
          weight_arg <- NULL
        }
        
        # Kernel Density Estimation
        kd_result <- density.ppp(points_ppp, 
                                 sigma = sigma_val, 
                                 positive = TRUE,
                                 kernel = "gaussian", 
                                 weights = weight_arg,
                                 dimyx = 500, 
                                 diggle = TRUE)
        
        # Convert to SpatRaster
        raster_kd <- rast(kd_result)
        crs(raster_kd) <- shapefile_crs
        
        # Normalize by month max (NO LOG TRANSFORM)
        max_val <- max(values(raster_kd), na.rm = TRUE)
        if (!is.na(max_val) && max_val > 0) {
          raster_kd <- raster_kd / max_val
        }
        
        # Update global min/max for this species
        local_min <- min(values(raster_kd), na.rm = TRUE)
        local_max <- max(values(raster_kd), na.rm = TRUE)
        global_min_species <- min(global_min_species, local_min)
        global_max_species <- max(global_max_species, local_max)
        
        # Save raster
        raster_filename <- paste0(output_dir, "KDE_", output_prefix, "_", 
                                  sanitize_filename(species), "_month", 
                                  sprintf("%02d", month), ".tif")
        writeRaster(raster_kd, filename = raster_filename, overwrite = TRUE)
        
        # Store raster
        monthly_rasters[[month]] <- raster_kd
        
      }, error = function(e) {
        cat("  Error processing month", month, ":", e$message, "\n")
      })
    }
    
    #.............................
    # SECOND PASS: Create plots with consistent scale------
    #.............................
    
    for(month in 1:12) {
      if(!is.null(monthly_rasters[[month]])) {
        raster_kd <- monthly_rasters[[month]]
        
        # Convert raster to data frame for ggplot
        kde_df <- as.data.frame(raster_kd, xy = TRUE) %>%
          filter(!is.na(lyr.1))
        
        # Create ggplot using species-wide global range for colour scale
        p <- ggplot(kde_df, aes(x = x, y = y, fill = lyr.1)) +
          geom_tile() +
          scale_fill_viridis_c(limits = c(global_min_species, global_max_species),
                               option = "viridis") +
          theme_minimal() +
          theme(axis.title = element_blank(),
                axis.ticks = element_blank(),
                axis.text = element_blank(), 
                axis.line = element_blank()) +
          labs(fill = "Density", title = month_names[month]) +
          annotate("text", x = Inf, y = Inf, 
                   label = paste("Bandwidth:", round(sigma_val, 0)),
                   hjust = 1.5, vjust = 2, size = 3, color = "orange") +
          coord_equal()
        
        monthly_plots[[month]] <- p
      }
    }
    
    # Combine all monthly plots for this species
    if(length(monthly_plots) > 0) {
      plot_combined <- wrap_plots(monthly_plots, ncol = 4, guides = "collect") +
        plot_annotation(title = paste("Monthly KDE -", species_display),
                        theme = theme(plot.title = element_text(size = 16, face = "bold")))
      
      # Save combined plot
      combined_filename <- paste0(plot_output_dir, "combined_monthly_KDE_", 
                                  output_prefix, "_", sanitize_filename(species), ".png")
      ggsave(filename = combined_filename, 
             plot = plot_combined, 
             width = 20, height = 15, dpi = 300)
      
      cat("Saved combined plot:", combined_filename, "\n")
    }
  }
  
  cat("\n=== All species processing complete ===\n")
}


#.............................
# LOAD STUDY AREA AND CREATE COMMON WINDOW-------
#.............................

cat("\n=== LOADING STUDY AREA AND CREATING COMMON WINDOW ===\n\n")

# Load study area shapefile
study_area_path <- "shapefiles/studyArea/ESS_study_area_simple.shp"
study_area <- st_read(study_area_path, quiet = TRUE)

cat("Study area CRS:", st_crs(study_area)$input, "\n")

# Transform to UTM20 if needed
if(st_crs(study_area)$epsg != 32620) {
  study_area <- st_transform(study_area, UTM20)
}

# Create window from study area extents
study_bbox <- st_bbox(study_area)

buffer_percent <- 0.05  # Small buffer for edge cases
x_range <- c(study_bbox["xmin"], study_bbox["xmax"])
y_range <- c(study_bbox["ymin"], study_bbox["ymax"])

x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
y_buffer <- (y_range[2] - y_range[1]) * buffer_percent

expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)

common_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)

cat("Common window created from study area\n")
cat("X range:", expanded_x_range, "\n")
cat("Y range:", expanded_y_range, "\n\n")

#.............................
# PROCESS MONTHLY PAM DATA---------
#.............................

cat("\n=== PROCESSING MONTHLY PAM DATA ===\n\n")

# Read baleen whale PAM data
filepath <- "input/2025/baleen_presence_laura_2025.csv"
baleen_PA_raw <- read.csv(filepath)

# Transform data: calculate detection rate per station-month-species
baleen_PA <- baleen_PA_raw %>%
  mutate(rec_date = as.Date(rec_date),
         month = month(rec_date),
         year = year(rec_date)) %>%
  group_by(site, latitude, longitude, species, month) %>%
  summarise(
    total_days = n(),
    detection_days = sum(presence, na.rm = TRUE),
    proportion_det = detection_days / total_days,
    .groups = 'drop'
  ) %>%
  filter(total_days > 0)

cat("Data summary (before clipping):\n")
cat("  Stations:", length(unique(baleen_PA$site)), "\n")
cat("  Species:", paste(unique(baleen_PA$species), collapse = ", "), "\n")
cat("  Station-month records:", nrow(baleen_PA), "\n\n")

# Create spatial object
baleen_sf_full <- baleen_PA %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = UTM20)

# Clip PAM data to study area
cat("Clipping PAM data to study area...\n")
baleen_sf <- st_intersection(baleen_sf_full, study_area)

cat("  Before clipping:", nrow(baleen_sf_full), "station-month records\n")
cat("  After clipping:", nrow(baleen_sf), "station-month records\n")
cat("  Removed:", nrow(baleen_sf_full) - nrow(baleen_sf), "records outside study area\n\n")

# Summary by species after clipping
species_summary <- baleen_sf %>%
  st_drop_geometry() %>%
  group_by(species) %>%
  summarise(
    n_records = n(),
    n_months = n_distinct(month),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

cat("Records by species (after clipping):\n")
print(species_summary)
cat("\n")

species_list_pam <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")

# Create species name mapping
species_names <- list(
  "Bb" = "Sei Whale",
  "Bm" = "Blue Whale",
  "Bp" = "Fin Whale",
  "Mn" = "Humpback Whale",
  "Ba" = "Minke Whale",
  "Eg" = "North Atlantic Right Whale"
)

cat("Processing", length(species_list_pam), "species from PAM data\n")
cat("Creating monthly KDEs for each species\n\n")

# Run monthly KDE for PAM using common window
performMonthlyKDE(
  data_sf = baleen_sf, 
  species_col = "species",
  month_col = "month",
  species_list = species_list_pam,
  species_names = species_names,
  weight_col = "proportion_det",
  window = common_window,
  output_prefix = "pam_clipped",
  sigma_val = 10000
)

cat("\n=== MONTHLY ANALYSIS COMPLETE ===\n")
cat("✓ Monthly PAM data processed for all species\n")
cat("\nOutputs saved to:\n")
cat("  Rasters:", raster_dir, "\n")
cat("  Plots:", plot_dir, "\n")