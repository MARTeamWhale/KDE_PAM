# HEADER --------------------------------------------
#
# Author: Laura Joan Feyrer
# Email:  ljfeyrer@dal.ca
# Date Created: 2026-01-12
# Most Recent Date Updated: 2026-01-15
#
# Script Name: KDES_CE_Collab.R
#
# Description:
## Comprehensive KDE analysis script for collaborative deliverables
## Generates all KDEs (PAM monthly, PAM overall, species, guilds) with:
##   - Organized folder structure
##   - 0.90 contour shapefiles
##   - Complete metadata index CSV
##   - Quality check maps with coastline, OWA, and PAM stations
#
# Outputs:
##   - /collab_CE/shapefiles/study_area/
##   - /collab_CE/kde_surfaces/PAM/monthly/
##   - /collab_CE/kde_surfaces/PAM/overall/
##   - /collab_CE/kde_surfaces/species/
##   - /collab_CE/kde_surfaces/guilds/
##   - /collab_CE/KDE_maps/ (QC visualizations)
##   - /collab_CE/kde_index.csv (master metadata file)
#

# SET OPTIONS ---------------------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(scipen = 999)
options(encoding = "UTF-8")

# Load necessary libraries
pacman::p_load(sf, tidyverse, readxl, here, leaflet, 
               scales, terra, ggrepel, viridis, ggspatial, spatstat, 
               dplyr, patchwork, lubridate)

# PROJECTIONS
UTM20 <- crs("+init=epsg:32620")
WGS84 <- crs("+init=epsg:4326")

# ................................
# SETUP OUTPUT DIRECTORY STRUCTURE------
# ................................

cat("\n=== SETTING UP OUTPUT DIRECTORY STRUCTURE ===\n\n")

base_dir <- "collab_CE"

# Clean existing output directory if it exists
if (dir.exists(base_dir)) {
  cat("Removing existing output directory...\n")
  unlink(base_dir, recursive = TRUE)
  cat("Cleaned:", base_dir, "\n")
}

# Define directory structure with separate rasters and shapes
dirs <- list(
  base = base_dir,
  study_area = file.path(base_dir, "shapefiles", "study_area"),
  
  # Species - separate rasters and shapes
  species_rasters = file.path(base_dir, "kde_surfaces", "species", "rasters"),
  species_shapes = file.path(base_dir, "kde_surfaces", "species", "shapes"),
  
  # Guilds - separate rasters and shapes
  guilds_rasters = file.path(base_dir, "kde_surfaces", "guilds", "rasters"),
  guilds_shapes = file.path(base_dir, "kde_surfaces", "guilds", "shapes"),
  
  # PAM overall - separate rasters and shapes
  pam_overall_rasters = file.path(base_dir, "kde_surfaces", "PAM", "overall", "rasters"),
  pam_overall_shapes = file.path(base_dir, "kde_surfaces", "PAM", "overall", "shapes"),
  
  # PAM monthly - separate rasters and shapes
  pam_monthly_rasters = file.path(base_dir, "kde_surfaces", "PAM", "monthly", "rasters"),
  pam_monthly_shapes = file.path(base_dir, "kde_surfaces", "PAM", "monthly", "shapes"),
  
  # Maps
  maps = file.path(base_dir, "KDE_maps")
)

# Create all directories
for (dir_name in names(dirs)) {
  if (!dir.exists(dirs[[dir_name]])) {
    dir.create(dirs[[dir_name]], recursive = TRUE)
    cat("Created:", dirs[[dir_name]], "\n")
  }
}

cat("\n")

# Initialize metadata index
metadata_index <- tibble(
  dataset = character(),
  species_group = character(),
  month = character(),
  point_count = integer(),
  kde_raster_filename = character(),
  contour_90_filename = character(),
  crs = character(),
  resolution_m = numeric(),
  bandwidth_m = numeric(),
  notes = character()
)

# ................................
# CONFIGURATION PARAMETERS------
# ................................

# Bandwidth settings
SIGHTINGS_SIGMA <- 10000  # 10km for sightings
PAM_SIGMA <- 10000         # 10km for PAM

# Proximity weighting for individual species (not guilds)
USE_PROXIMITY_WEIGHTING <- TRUE
PROXIMITY_K_NEIGHBORS <- 8
PROXIMITY_MAX_DISTANCE <- 5000
PROXIMITY_MIN_WEIGHT <- 0.05

# ................................
# HELPER FUNCTIONS------
# ................................

# Sanitize filename
sanitize_filename <- function(name) {
  name <- trimws(name)
  name <- gsub("[/\\:*?\"<>|]", "_", name)
  name <- gsub(" ", "_", name)
  return(name)
}

# Calculate proximity weights
calculate_proximity_weights <- function(data_sf, 
                                        k_neighbors = PROXIMITY_K_NEIGHBORS, 
                                        max_distance = PROXIMITY_MAX_DISTANCE,
                                        min_weight = PROXIMITY_MIN_WEIGHT) {
  
  coords <- st_coordinates(data_sf)
  dist_matrix <- as.matrix(dist(coords))
  
  proximity_scores <- apply(dist_matrix, 1, function(row) {
    sorted_dists <- sort(row)
    nearest_k <- sorted_dists[2:min(k_neighbors + 1, length(sorted_dists))]
    mean(nearest_k)
  })
  
  weights <- 1 / (1 + (proximity_scores / max_distance))
  weights <- (weights - min(weights)) / (max(weights) - min(weights))
  weights <- pmax(weights, min_weight)
  
  return(weights)
}

# Convert KDE raster to 0.90 contour shapefile
create_contour_90 <- function(kde_raster, output_path) {
  
  # Set values <= 0 to NA
  kde_raster[kde_raster <= 0] <- NA
  
  # Calculate 0.90 quantile
  kde_values <- values(kde_raster, mat = FALSE, na.rm = TRUE)
  threshold_90 <- quantile(kde_values, probs = 0.90, type = 6, na.rm = TRUE)
  
  # Create binary raster (1 = above threshold, NA = below)
  kde_binary <- kde_raster
  kde_binary[kde_binary < threshold_90] <- NA
  kde_binary[!is.na(kde_binary)] <- 1
  
  # Convert to polygon
  kde_poly <- as.polygons(kde_binary, dissolve = TRUE)
  kde_sf <- st_as_sf(kde_poly)
  kde_sf <- st_make_valid(kde_sf)
  
  # Add attribute
  kde_sf <- kde_sf %>%
    mutate(contour = "0.90") %>%
    select(contour, geometry)
  
  # Save shapefile
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  return(invisible(kde_sf))
}

# ................................
# UNIFIED KDE FUNCTION------
# ................................

performKDE_complete <- function(data_sf, 
                                species_col,
                                species_list, 
                                dataset_name,
                                weight_col = NULL,
                                use_proximity = FALSE,
                                sigma_val = 10000,
                                window = NULL,
                                raster_dir,
                                shape_dir,
                                month_val = NULL,
                                resolution = 500) {  # Add resolution parameter (500m default)
  
  shapefile_crs <- st_crs(data_sf)$wkt
  results_list <- list()
  
  for(species in species_list) {
    current_sf <- data_sf[data_sf[[species_col]] == species, ]
    
    if(nrow(current_sf) == 0) {
      cat("  No data for", species, "\n")
      next
    }
    
    cat("  Processing:", species, "(n =", nrow(current_sf), ")\n")
    
    tryCatch({
      species_coords <- st_coordinates(current_sf)
      
      # Create window
      if(is.null(window)) {
        x_range <- range(species_coords[, "X"])
        y_range <- range(species_coords[, "Y"])
        x_buffer <- (x_range[2] - x_range[1]) * 0.2
        y_buffer <- (y_range[2] - y_range[1]) * 0.2
        expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
        expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
        species_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
      } else {
        species_window <- window
      }
      
      # Apply proximity weighting if requested
      if(use_proximity && is.null(weight_col)) {
        current_sf$proximity_weight <- calculate_proximity_weights(current_sf)
        weight_col_use <- "proximity_weight"
      } else {
        weight_col_use <- weight_col
      }
      
      # Create ppp object
      if(!is.null(weight_col_use)) {
        weights <- current_sf[[weight_col_use]] * 10000
        points_ppp <- ppp(x = species_coords[,"X"], 
                          y = species_coords[,"Y"], 
                          window = species_window,
                          marks = data.frame(weights = weights))
      } else {
        points_ppp <- ppp(x = species_coords[,"X"], 
                          y = species_coords[,"Y"], 
                          window = species_window)
      }
      
      # Handle duplicates
      points_ppp <- rjitter(points_ppp, retry = TRUE, nsim = 1, drop = TRUE)
      
      # Extract weights after jittering
      if(!is.null(weight_col_use)) {
        weight_arg <- marks(points_ppp)
      } else {
        weight_arg <- NULL
      }
      
      # Calculate dimensions for specified resolution
      x_extent <- diff(species_window$xrange)
      y_extent <- diff(species_window$yrange)
      dimx <- round(x_extent / resolution)
      dimy <- round(y_extent / resolution)
      
      # Kernel Density Estimation with specified resolution
      kd_result <- density.ppp(points_ppp, 
                               sigma = sigma_val, 
                               positive = TRUE,
                               kernel = "gaussian", 
                               weights = weight_arg,
                               dimyx = c(dimy, dimx),  # Set explicit dimensions for 500m resolution
                               diggle = TRUE)
      
      # Convert to SpatRaster
      raster_kd <- rast(kd_result)
      crs(raster_kd) <- shapefile_crs
      
      # Normalize by max
      max_val <- max(values(raster_kd), na.rm = TRUE)
      if (!is.na(max_val) && max_val > 0) {
        raster_kd <- raster_kd / max_val
      }
      
      # Generate filenames
      species_clean <- sanitize_filename(species)
      
      if(!is.null(month_val)) {
        month_str <- sprintf("_month%02d", month_val)
        raster_filename <- paste0(species_clean, month_str, ".tif")
        contour_filename <- paste0(species_clean, month_str, "_contour90.shp")
      } else {
        raster_filename <- paste0(species_clean, ".tif")
        contour_filename <- paste0(species_clean, "_contour90.shp")
      }
      
      raster_path <- file.path(raster_dir, raster_filename)
      contour_path <- file.path(shape_dir, contour_filename)
      
      # Save raster
      writeRaster(raster_kd, filename = raster_path, overwrite = TRUE)
      
      # Create and save 0.90 contour
      create_contour_90(raster_kd, contour_path)
      
      # Get actual resolution
      actual_resolution <- mean(res(raster_kd))
      
      # Add to metadata index
      month_str_index <- if(!is.null(month_val)) as.character(month_val) else "all"
      
      new_row <- tibble(
        dataset = dataset_name,
        species_group = species,
        month = month_str_index,
        point_count = nrow(current_sf),
        kde_raster_filename = raster_filename,
        contour_90_filename = contour_filename,
        crs = "EPSG:32620",
        resolution_m = round(actual_resolution, 2),
        bandwidth_m = sigma_val,
        notes = if(use_proximity) "proximity weighted" else if(!is.null(weight_col)) "detection weighted" else "unweighted"
      )
      
      metadata_index <<- bind_rows(metadata_index, new_row)
      
      # Store result
      results_list[[species]] <- list(
        raster = raster_kd,
        contour = st_read(contour_path, quiet = TRUE),
        n_points = nrow(current_sf)
      )
      
    }, error = function(e) {
      cat("  ERROR processing", species, ":", e$message, "\n")
    })
  }
  
  return(invisible(results_list))
}

# ................................
# LOAD BASE SPATIAL DATA------
# ................................

cat("\n=== LOADING BASE SPATIAL DATA ===\n\n")

# Study area
study_area_path <- "shapefiles/studyArea/ESS_study_area_simple.shp"
study_area <- st_read(study_area_path, quiet = TRUE) %>%
  st_transform(UTM20)

# Copy study area to output
write_sf(study_area, file.path(dirs$study_area, "ESS_study_area.shp"), delete_dsn = TRUE)
cat("Copied study area shapefile\n")

# Coastline (for visualization)
coastline_path <- "shapefiles/Canada/NE_10mLand.shp"
if(file.exists(coastline_path)) {
  coastline <- st_read(coastline_path, quiet = TRUE) %>%
    st_transform(UTM20) %>%
    st_crop(st_bbox(study_area))
  cat("Loaded coastline\n")
} else {
  coastline <- NULL
  cat("Coastline file not found - maps will not include coastline\n")
}

# OWA polygons (for visualization)
owa_path <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"
if(file.exists(owa_path)) {
  owa <- st_read(owa_path, quiet = TRUE) %>%
    st_transform(UTM20) %>%
    st_crop(st_bbox(study_area))
  cat("Loaded offshore wind areas\n")
} else {
  owa <- NULL
  cat("OWA file not found - maps will not include OWA\n")
}
if(file.exists(owa_path)) {
  owa <- st_read(owa_path, quiet = TRUE) %>%
    st_transform(UTM20) %>%
    st_crop(st_bbox(study_area))
  cat("Loaded offshore wind areas\n")
} else {
  owa <- NULL
  cat("OWA file not found - maps will not include OWA\n")
}

# Create common window
study_bbox <- st_bbox(study_area)
buffer_percent <- 0.05
x_range <- c(study_bbox["xmin"], study_bbox["xmax"])
y_range <- c(study_bbox["ymin"], study_bbox["ymax"])
x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
common_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)

cat("Common window created\n\n")

# ................................
# PROCESS SIGHTINGS DATA------
# ................................

cat("\n=== PROCESSING SIGHTINGS DATA ===\n\n")

# Load combined sightings
combined_data <- read_csv("input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv",
                          show_col_types = FALSE) %>%
  filter(as.Date(date_utc) >= as.Date("2015-01-01"))

cat("Total sightings records (2015+):", nrow(combined_data), "\n")

# Create sf object
combined_sf <- st_as_sf(combined_data, 
                        coords = c("lon", "lat"), 
                        crs = 4326) %>%
  st_transform(UTM20) %>%
  st_make_valid()

# Filter to cetaceans only
cetacean_sf <- combined_sf %>%
  filter(cetacean_family %in% c("odontocete", "mysticete")) %>%
  mutate(common_name = case_when(
    str_detect(common_name, "SEI") ~ "WHALE-SEI",
    TRUE ~ common_name
  ))

cat("Cetacean records:", nrow(cetacean_sf), "\n")

# Get species counts
species_counts <- cetacean_sf %>%
  st_drop_geometry() %>%
  group_by(common_name, cetacean_family) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  arrange(desc(n_records))

cat("\nSpecies counts:\n")
print(species_counts, n = Inf)

# Filter species with >= 25 records
species_25plus <- species_counts %>%
  filter(n_records >= 25) %>%
  pull(common_name)

cat("\n", length(species_25plus), "species with ≥25 records\n")

#.............................
# INDIVIDUAL SPECIES KDEs------
#.............................

cat("\n=== CREATING INDIVIDUAL SPECIES KDEs ===\n\n")

species_sf <- cetacean_sf %>%
  filter(common_name %in% species_25plus)

performKDE_complete(
  data_sf = species_sf,
  species_col = "common_name",
  species_list = species_25plus,
  dataset_name = "sightings_individual",
  weight_col = NULL,
  use_proximity = USE_PROXIMITY_WEIGHTING,
  sigma_val = SIGHTINGS_SIGMA,
  window = common_window,
  raster_dir = dirs$species_rasters,
  shape_dir = dirs$species_shapes,
  resolution = 500
)

cat("✓ Individual species KDEs complete\n")

#.............................
# GUILD-LEVEL KDEs------
#.............................

cat("\n=== CREATING GUILD-LEVEL KDEs ===\n\n")

# Odontocetes
odontocete_sf <- cetacean_sf %>%
  filter(cetacean_family == "odontocete") %>%
  mutate(guild = "Odontocetes")

cat("Processing Odontocetes (n =", nrow(odontocete_sf), ")\n")

performKDE_complete(
  data_sf = odontocete_sf,
  species_col = "guild",
  species_list = "Odontocetes",
  dataset_name = "sightings_guild",
  weight_col = NULL,
  use_proximity = FALSE,  # No proximity weighting for guilds
  sigma_val = SIGHTINGS_SIGMA,
  window = common_window,
  raster_dir = dirs$guilds_rasters,
  shape_dir = dirs$guilds_shapes,
  resolution = 500
)

# Mysticetes
mysticete_sf <- cetacean_sf %>%
  filter(cetacean_family == "mysticete") %>%
  mutate(guild = "Mysticetes")

cat("Processing Mysticetes (n =", nrow(mysticete_sf), ")\n")

performKDE_complete(
  data_sf = mysticete_sf,
  species_col = "guild",
  species_list = "Mysticetes",
  dataset_name = "sightings_guild",
  weight_col = NULL,
  use_proximity = FALSE,  # No proximity weighting for guilds
  sigma_val = SIGHTINGS_SIGMA,
  window = common_window,
  raster_dir = dirs$guilds_rasters,
  shape_dir = dirs$guilds_shapes,
  resolution = 500
)

cat("✓ Guild-level KDEs complete\n")

# ................................
# PROCESS PAM DATA------
# ................................

cat("\n=== PROCESSING PAM DATA ===\n\n")

# Load PAM data
filepath <- 'input/2025/baleen_presence_laura_2025.csv'
baleen_PA_raw <- read.csv(filepath)

# Calculate station-month summaries
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

# Create spatial object
baleen_sf_full <- baleen_PA %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = UTM20)

# Clip to study area
baleen_sf <- st_intersection(baleen_sf_full, study_area)

cat("PAM records (after clipping):", nrow(baleen_sf), "\n")

# Get unique stations for visualization
# baleen_sf is already spatial, so just extract one point per site
pam_stations <- baleen_sf %>%
  group_by(site) %>%
  slice(1) %>%
  ungroup() %>%
  select(site, geometry)

cat("Unique PAM stations:", nrow(pam_stations), "\n")

# Species list
species_list_pam <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")
species_names_pam <- c(
  "Bb" = "WHALE-SEI",
  "Bm" = "WHALE-BLUE",
  "Bp" = "WHALE-FIN",
  "Mn" = "WHALE-HUMPBACK",
  "Ba" = "WHALE-MINKE",
  "Eg" = "WHALE-NORTH_ATLANTIC_RIGHT"
)

#.............................
# PAM MONTHLY KDEs------
#.............................

cat("\n=== CREATING PAM MONTHLY KDEs ===\n\n")

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

for(species in species_list_pam) {
  species_common <- species_names_pam[species]
  cat("\nProcessing", species_common, "(", species, ")\n")
  
  species_data <- baleen_sf %>%
    filter(species == !!species) %>%
    mutate(species_common = species_common)  # Add common name column
  
  for(month in 1:12) {
    month_data <- species_data %>%
      filter(month == !!month)
    
    if(nrow(month_data) > 0) {
      cat("  ", month_names[month], ": n =", nrow(month_data), "\n")
      
      performKDE_complete(
        data_sf = month_data,
        species_col = "species_common",  # Use common name column
        species_list = species_common,    # Use common name
        dataset_name = "PAM_monthly",
        weight_col = "proportion_det",
        use_proximity = FALSE,
        sigma_val = PAM_SIGMA,
        window = common_window,
        raster_dir = dirs$pam_monthly_rasters,
        shape_dir = dirs$pam_monthly_shapes,
        month_val = month,
        resolution = 500
      )
    }
  }
}

cat("\n✓ PAM monthly KDEs complete\n")

#.............................
# PAM OVERALL KDEs (all months combined)------
#.............................

cat("\n=== CREATING PAM OVERALL KDEs ===\n\n")

# For overall PAM, aggregate across all months
baleen_overall <- baleen_sf %>%
  group_by(site, species) %>%
  summarise(
    total_days = sum(total_days),
    detection_days = sum(detection_days),
    proportion_det = detection_days / total_days,
    .groups = 'drop'
  ) %>%
  st_as_sf()

for(species in species_list_pam) {
  species_common <- species_names_pam[species]
  cat("Processing", species_common, "(", species, ")\n")
  
  species_data <- baleen_overall %>%
    filter(species == !!species) %>%
    mutate(species_common = species_common)  # Add common name column
  
  if(nrow(species_data) > 0) {
    performKDE_complete(
      data_sf = species_data,
      species_col = "species_common",  # Use common name column
      species_list = species_common,    # Use common name
      dataset_name = "PAM_overall",
      weight_col = "proportion_det",
      use_proximity = FALSE,
      sigma_val = PAM_SIGMA,
      window = common_window,
      raster_dir = dirs$pam_overall_rasters,
      shape_dir = dirs$pam_overall_shapes,
      month_val = NULL,
      resolution = 500
    )
  }
}

cat("✓ PAM overall KDEs complete\n")

# ................................
# SAVE METADATA INDEX------
# ................................

cat("\n=== SAVING METADATA INDEX ===\n\n")

index_path <- file.path(base_dir, "kde_index.csv")
write_csv(metadata_index, index_path)

cat("Metadata index saved:", index_path, "\n")
cat("Total entries:", nrow(metadata_index), "\n\n")

# ................................
# CREATE QUALITY CHECK MAPS------
# ................................

cat("\n=== CREATING QUALITY CHECK MAPS ===\n\n")

# Helper function to create a map with improved scaling
create_kde_map <- function(raster_path, title, output_filename, 
                           show_stations = FALSE,
                           show_contour = FALSE,
                           contour_path = NULL,
                           show_points = FALSE,
                           points_sf = NULL) {
  
  # Load raster
  r <- rast(raster_path)
  
  # Convert to data frame, apply square root transform for better visibility
  # Square root is less extreme than log and still helps reveal low values
  r_df <- as.data.frame(r, xy = TRUE) %>%
    filter(!is.na(lyr.1)) %>%
    mutate(lyr.1_sqrt = sqrt(lyr.1))  # Square root transform
  
  # Create base plot
    # Add KDE tiles
  p <- ggplot() + geom_tile(data = r_df, aes(x = x, y = y, fill = lyr.1_sqrt)) +
    scale_fill_viridis_c(option = "viridis", 
                         name = "sqrt(Density)",
                         limits = c(0, max(r_df$lyr.1_sqrt)),
                         begin = 0.15,  # Trim dark end of palette (0.15 instead of 0)
                         end = 1.0) +
    
  # Add coastline on top of raster
   geom_sf(data = coastline, fill = "grey90", color = "grey70", size = 0.3)+
    theme_minimal() +
    labs(title = title) +
    coord_sf(crs = st_crs(UTM20))
  
  # Add OWA
  if(!is.null(owa)) {
    p <- p + geom_sf(data = owa, fill = NA, color = "red", size = 0.5, linetype = "dashed")
  }
  
  # Add 0.90 contour if provided
  if(show_contour && !is.null(contour_path) && file.exists(contour_path)) {
    contour_sf <- st_read(contour_path, quiet = TRUE)
    p <- p + geom_sf(data = contour_sf, fill = NA, color = "white", size = 1, linetype = "solid")
  }
  
  # Add sighting points if provided
  if(show_points && !is.null(points_sf)) {
    p <- p + geom_sf(data = points_sf, color = "white", size = 0.8, 
                     shape = 21, fill = "black", stroke = 0.3, alpha = 0.6)
  }
  
  # Add PAM stations if requested
  if(show_stations && exists("pam_stations")) {
    p <- p + geom_sf(data = pam_stations, color = "white", size = 1.5, shape = 21, 
                     fill = "black", stroke = 0.5, alpha = 0.8)
  }
  
  # Add study area outline (top layer)
  p <- p + geom_sf(data = study_area, fill = NA, color = "black", size = 0.8)
  
  # Theme adjustments
  p <- p + theme(
    axis.title = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 12, face = "bold"),
    panel.grid = element_line(color = "grey95")
  )
  
  # Save
  ggsave(output_filename, p, width = 10, height = 8, dpi = 150)
  
  return(invisible(p))
}

# Create maps for each dataset type

# 1. Individual species sightings (with contours and points)
cat("Creating maps for individual species...\n")
species_rasters <- list.files(dirs$species_rasters, pattern = "\\.tif$", full.names = TRUE)

for(raster_file in species_rasters) {
  species_name <- tools::file_path_sans_ext(basename(raster_file))
  
  # Get corresponding contour shapefile
  contour_file <- file.path(dirs$species_shapes, paste0(species_name, "_contour90.shp"))
  
  # Get sighting points for this species
  species_points <- species_sf %>%
    filter(sanitize_filename(common_name) == species_name)
  
  map_filename <- file.path(dirs$maps, paste0("map_species_", species_name, ".png"))
  create_kde_map(raster_file, 
                 title = paste("Sightings KDE:", gsub("_", " ", species_name)),
                 output_filename = map_filename,
                 show_stations = FALSE,
                 show_contour = TRUE,
                 contour_path = contour_file,
                 show_points = TRUE,
                 points_sf = species_points)
}

# 2. Guilds (with contours)
cat("Creating maps for guilds...\n")
guild_rasters <- list.files(dirs$guilds_rasters, pattern = "\\.tif$", full.names = TRUE)
for(raster_file in guild_rasters) {
  guild_name <- tools::file_path_sans_ext(basename(raster_file))
  
  # Get corresponding contour shapefile
  contour_file <- file.path(dirs$guilds_shapes, paste0(guild_name, "_contour90.shp"))
  
  map_filename <- file.path(dirs$maps, paste0("map_guild_", guild_name, ".png"))
  create_kde_map(raster_file, 
                 title = paste("Guild KDE:", gsub("_", " ", guild_name)),
                 output_filename = map_filename,
                 show_stations = FALSE,
                 show_contour = TRUE,
                 contour_path = contour_file,
                 show_points = FALSE)
}

# 3. PAM overall (with contours and stations)
cat("Creating maps for PAM overall...\n")
pam_overall_rasters <- list.files(dirs$pam_overall_rasters, pattern = "\\.tif$", full.names = TRUE)
for(raster_file in pam_overall_rasters) {
  species_name <- tools::file_path_sans_ext(basename(raster_file))
  
  # Get corresponding contour shapefile
  contour_file <- file.path(dirs$pam_overall_shapes, paste0(species_name, "_contour90.shp"))
  
  map_filename <- file.path(dirs$maps, paste0("map_PAM_overall_", species_name, ".png"))
  create_kde_map(raster_file, 
                 title = paste("PAM Overall:", gsub("_", " ", species_name)),
                 output_filename = map_filename,
                 show_stations = TRUE,
                 show_contour = TRUE,
                 contour_path = contour_file,
                 show_points = FALSE)
}

# 4. PAM monthly (create faceted 12-month plots per species)
cat("Creating faceted monthly maps for PAM (one plot per species, 12 months)...\n")

# Get all PAM monthly rasters
pam_monthly_rasters <- list.files(dirs$pam_monthly_rasters, pattern = "\\.tif$", full.names = TRUE)

# Group by species
for(species_code in species_list_pam) {
  species_common <- species_names_pam[species_code]
  species_clean <- sanitize_filename(species_common)
  
  cat("Creating monthly facet plot for", species_common, "...\n")
  
  # Get all months for this species
  species_rasters <- pam_monthly_rasters[grepl(paste0(species_clean, "_month"), pam_monthly_rasters)]
  
  if(length(species_rasters) == 0) {
    cat("  No rasters found for", species_common, "\n")
    next
  }
  
  # Prepare data for all months
  monthly_data_list <- list()
  
  for(raster_file in species_rasters) {
    filename_base <- tools::file_path_sans_ext(basename(raster_file))
    month_num <- as.numeric(gsub(".*month(\\d+).*", "\\1", filename_base))
    
    # Load raster and convert to dataframe
    r <- rast(raster_file)
    r_df <- as.data.frame(r, xy = TRUE) %>%
      filter(!is.na(lyr.1)) %>%
      mutate(
        lyr.1_sqrt = sqrt(lyr.1),
        month = month_num,
        month_name = factor(month_names[month_num], levels = month_names)
      )
    
    monthly_data_list[[month_num]] <- r_df
  }
  
  # Combine all months
  all_months_df <- bind_rows(monthly_data_list)
  
  # Get global scale limits for this species across all months
  global_min <- min(all_months_df$lyr.1_sqrt, na.rm = TRUE)
  global_max <- max(all_months_df$lyr.1_sqrt, na.rm = TRUE)
  
  # Create faceted plot
  p <- ggplot()
  
  # Add coastline first (bottom layer)
  if(!is.null(coastline)) {
    p <- p + geom_sf(data = coastline, fill = "grey90", color = "grey70", size = 0.2)
  }
  
  # Add KDE tiles
  p <- p + geom_tile(data = all_months_df, aes(x = x, y = y, fill = lyr.1_sqrt)) +
    scale_fill_viridis_c(option = "viridis", 
                         name = "sqrt(Density)",
                         limits = c(global_min, global_max),
                         begin = 0.15,
                         end = 1.0)
  
  # Add OWA
  if(!is.null(owa)) {
    p <- p + geom_sf(data = owa, fill = NA, color = "red", size = 0.3, linetype = "dashed")
  }
  
  # Add contours for each month
  for(raster_file in species_rasters) {
    filename_base <- tools::file_path_sans_ext(basename(raster_file))
    month_num <- as.numeric(gsub(".*month(\\d+).*", "\\1", filename_base))
    contour_file <- file.path(dirs$pam_monthly_shapes, paste0(filename_base, "_contour90.shp"))
    
    if(file.exists(contour_file)) {
      contour_sf <- st_read(contour_file, quiet = TRUE) %>%
        mutate(month_name = factor(month_names[month_num], levels = month_names))
      
      p <- p + geom_sf(data = contour_sf, aes(geometry = geometry), 
                       fill = NA, color = "white", size = 0.6, linetype = "solid")
    }
  }
  
  # Add PAM stations
  if(exists("pam_stations")) {
    p <- p + geom_sf(data = pam_stations, color = "white", size = 0.8, shape = 21, 
                     fill = "black", stroke = 0.3, alpha = 0.7)
  }
  
  # Add study area outline
  p <- p + geom_sf(data = study_area, fill = NA, color = "black", size = 0.5)
  
  # Facet by month
  p <- p + facet_wrap(~ month_name, ncol = 4) +
    coord_sf(crs = st_crs(UTM20)) +
    theme_minimal() +
    labs(title = paste("PAM Monthly:", gsub("_", " ", species_common))) +
    theme(
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.grid = element_blank(),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey95", color = "grey70")
    )
  
  # Save
  map_filename <- file.path(dirs$maps, paste0("map_PAM_monthly_faceted_", species_clean, ".png"))
  ggsave(map_filename, p, width = 16, height = 12, dpi = 150)
  
  cat("  Saved:", map_filename, "\n")
}

cat("Created faceted monthly maps for", length(species_list_pam), "PAM species\n")

cat("✓ Quality check maps created\n")

# ................................
# CREATE README FILES------
# ................................

cat("\n=== CREATING README FILES ===\n\n")

# Main README
main_readme <- "# KDE Analysis Outputs for Collaborative Project

This directory contains Kernel Density Estimation (KDE) surfaces and associated shapefiles 
for cetacean distribution analysis in the Eastern Scotian Shelf study area.

## Directory Structure

```
collab_CE/
├── shapefiles/
│   └── study_area/           # Study area boundary shapefile
├── kde_surfaces/
│   ├── species/
│   │   ├── rasters/          # Individual species KDE rasters (GeoTIFF)
│   │   └── shapes/           # Individual species 0.90 contour shapefiles
│   ├── guilds/
│   │   ├── rasters/          # Guild-level KDE rasters (GeoTIFF)
│   │   └── shapes/           # Guild-level 0.90 contour shapefiles
│   └── PAM/
│       ├── overall/
│       │   ├── rasters/      # Overall PAM KDE rasters (all months)
│       │   └── shapes/       # Overall PAM 0.90 contour shapefiles
│       └── monthly/
│           ├── rasters/      # Monthly PAM KDE rasters
│           └── shapes/       # Monthly PAM 0.90 contour shapefiles
├── KDE_maps/                 # Quality check visualizations
└── kde_index.csv             # Master metadata file
```

## File Naming Convention

### Rasters (GeoTIFF)
- Species sightings: `{SPECIES_NAME}.tif`
- Guilds: `{GUILD_NAME}.tif`
- PAM overall: `{SPECIES_CODE}.tif`
- PAM monthly: `{SPECIES_CODE}_month{MM}.tif`

### Shapefiles (0.90 quantile contours)
- Corresponding shapefile: `{RASTER_NAME}_contour90.shp`
- Located in parallel `/shapes/` directories

## KDE Parameters

### Sightings Data
- **Bandwidth**: 10,000 m (10 km)
- **Kernel**: Gaussian
- **Resolution**: ~500 m
- **Weighting**: Proximity-weighted for individual species (isolated observations downweighted)
- **Normalization**: Normalized to [0, 1] by species maximum
- **CRS**: EPSG:32620 (WGS84 UTM Zone 20N)

### PAM Data
- **Bandwidth**: 10,000 m (10 km)
- **Kernel**: Gaussian
- **Resolution**: ~500 m
- **Weighting**: Detection rate (proportion of days with detections)
- **Normalization**: Normalized to [0, 1] by species maximum
- **CRS**: EPSG:32620 (WGS84 UTM Zone 20N)

## Species Codes (PAM)

- **Bb** - Sei Whale (*Balaenoptera borealis*)
- **Bm** - Blue Whale (*Balaenoptera musculus*)
- **Bp** - Fin Whale (*Balaenoptera physalus*)
- **Mn** - Humpback Whale (*Megaptera novaeangliae*)
- **Ba** - Minke Whale (*Balaenoptera acutorostrata*)
- **Eg** - North Atlantic Right Whale (*Eubalaena glacialis*)

## Metadata Index

The file `kde_index.csv` contains complete metadata for all KDE outputs including:
- Dataset type (sightings_individual, sightings_guild, PAM_monthly, PAM_overall)
- Species/group name
- Month (for monthly analyses)
- Number of input points/stations
- Raster filename
- Contour shapefile filename
- CRS
- Spatial resolution
- Bandwidth
- Weighting method

## Quality Check Maps

The `KDE_maps/` directory contains PNG visualizations of KDE surfaces with:
- Sqrt-transformed density values for improved visibility
- Coastline overlay
- Offshore Wind Area (OWA) boundaries
- PAM station locations (for PAM analyses)
- Study area boundary

## Data Sources

- **Sightings**: Combined deduplicated dataset from WSDB and aerial surveys (2015-present)
- **PAM**: Autonomous acoustic recorder deployments (baleen/ beaked whale presence)

## Contact

For questions about these outputs, contact:
Laura Joan Feyrer (ljfeyrer@dal.ca)

## Date Generated

"

main_readme <- paste0(main_readme, format(Sys.Date(), "%Y-%m-%d"), "\n")

writeLines(main_readme, file.path(base_dir, "README.md"))
cat("Main README created\n")

# Dataset-specific READMEs
species_readme <- "# Individual Species KDEs (Sightings Data)

This directory contains KDE surfaces for individual cetacean species with ≥25 sighting records.

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Processing Details

- **Data**: Combined sightings from WSDB and aerial surveys (2015-present)
- **Weighting**: Proximity-weighted (isolated observations receive lower weight)
- **Proximity parameters**:
  - K nearest neighbors: 8
  - Reference distance: 5000 m
  - Minimum weight: 0.05
- **Bandwidth**: 10,000 m (fixed)

## Files

Each species has:
- `rasters/{SPECIES}.tif` - KDE raster surface
- `shapes/{SPECIES}_contour90.shp` - 0.90 quantile contour polygon

See ../../../kde_index.csv for complete metadata.
"

writeLines(species_readme, file.path(dirs$species_rasters, "..", "README.md"))

guilds_readme <- "# Guild-Level KDEs (Sightings Data)

This directory contains KDE surfaces for cetacean guilds (functional groups).

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Guilds

- **Odontocetes**: All toothed whales
- **Mysticetes**: All baleen whales

## Processing Details

- **Data**: Combined sightings from WSDB and aerial surveys (2015-present)
- **Weighting**: None (guilds have naturally broad distributions)
- **Bandwidth**: 10,000 m (fixed)

## Files

Each guild has:
- `rasters/{GUILD}.tif` - KDE raster surface
- `shapes/{GUILD}_contour90.shp` - 0.90 quantile contour polygon

See ../../../kde_index.csv for complete metadata.
"

writeLines(guilds_readme, file.path(dirs$guilds_rasters, "..", "README.md"))

pam_overall_readme <- "# PAM Overall KDEs (All Months Combined)

This directory contains KDE surfaces for passive acoustic monitoring (PAM) detections 
aggregated across all months.

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Processing Details

- **Data**: Autonomous acoustic recorder deployments
- **Weighting**: Detection rate (proportion of deployment days with species detection)
- **Bandwidth**: 10,000 m (fixed)
- **Temporal coverage**: All months combined

## Species Codes

See main README for species code definitions.

## Files

Each species has:
- `rasters/{SPECIES_CODE}.tif` - KDE raster surface
- `shapes/{SPECIES_CODE}_contour90.shp` - 0.90 quantile contour polygon

See ../../../../kde_index.csv for complete metadata.
"

writeLines(pam_overall_readme, file.path(dirs$pam_overall_rasters, "..", "README.md"))

pam_monthly_readme <- "# PAM Monthly KDEs

This directory contains monthly KDE surfaces for passive acoustic monitoring (PAM) detections.

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Processing Details

- **Data**: Autonomous acoustic recorder deployments
- **Weighting**: Detection rate (proportion of deployment days with species detection)
- **Bandwidth**: 10,000 m (fixed)
- **Temporal resolution**: Monthly

## File Naming

Format: `{SPECIES_CODE}_month{MM}.tif`
- MM = 01 (Jan) through 12 (Dec)

## Species Codes

See main README for species code definitions.

## Files

Each species-month combination has:
- `rasters/{SPECIES_CODE}_month{MM}.tif` - KDE raster surface
- `shapes/{SPECIES_CODE}_month{MM}_contour90.shp` - 0.90 quantile contour polygon

See ../../../../kde_index.csv for complete metadata.
"

writeLines(pam_monthly_readme, file.path(dirs$pam_monthly_rasters, "..", "README.md"))

cat("✓ README files created\n")

# ................................
# FINAL SUMMARY------
# ................................

cat("\n")
cat("==============================================================================\n")
cat("ALL PROCESSING COMPLETE\n")
cat("==============================================================================\n\n")

cat("Summary of outputs:\n")
cat("  Individual species KDEs:", length(list.files(dirs$species_rasters, pattern = "\\.tif$")), "\n")
cat("  Guild KDEs:", length(list.files(dirs$guilds_rasters, pattern = "\\.tif$")), "\n")
cat("  PAM overall KDEs:", length(list.files(dirs$pam_overall_rasters, pattern = "\\.tif$")), "\n")
cat("  PAM monthly KDEs:", length(list.files(dirs$pam_monthly_rasters, pattern = "\\.tif$")), "\n")
cat("  Total KDE surfaces:", nrow(metadata_index), "\n")
cat("  Quality check maps:", length(list.files(dirs$maps, pattern = "\\.png$")), "\n\n")

cat("Main outputs:\n")
cat("  Metadata index:", index_path, "\n")
cat("  Main README:", file.path(base_dir, "README.md"), "\n")
cat("  Output directory:", base_dir, "\n\n")

cat("All files ready for collaborative delivery!\n")
cat("==============================================================================\n")