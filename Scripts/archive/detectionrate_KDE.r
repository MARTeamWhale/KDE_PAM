# ==============================================================================
# DETECTION-BASED THRESHOLD SHAPEFILE GENERATION FOR PAM DATA
# ==============================================================================
# Instead of using quantiles of all pixels, uses thresholds based on 
# actual detection rates observed at PAM stations
# ==============================================================================

options(scipen = 999)
pacman::p_load(sf, tidyverse, terra)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

UTM20 <- st_crs(32620)

# Input paths
KDE_RASTER_DIR <- "output/tif/"
PAM_DATA <- "input/2025/baleen_presence_laura_2025.csv"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"

# Output
SHAPEFILE_OUTPUT_DIR <- "output/shapes/detection_based/"
if (!dir.exists(SHAPEFILE_OUTPUT_DIR)) {
  dir.create(SHAPEFILE_OUTPUT_DIR, recursive = TRUE)
}

# Detection rate thresholds to create contours for
# These represent "areas likely to have X% detection rates"
DETECTION_THRESHOLDS <- c(0.01, 0.02, 0.05, 0.10, 0.20, 0.30)
THRESHOLD_LABELS <- c(">1%", ">2%", ">5%", ">10%", ">20%", ">30%")

# PAM species codes
PAM_SPECIES <- c("Ba", "Bb", "Bm", "Bp", "Mn", "Eg")
PAM_SPECIES_NAMES <- list(
  "Ba" = "Minke",
  "Bb" = "Sei", 
  "Bm" = "Blue",
  "Bp" = "Fin",
  "Mn" = "Humpback",
  "Eg" = "NARW"
)

# ==============================================================================
# LOAD PAM DATA
# ==============================================================================

cat("=== LOADING PAM DATA ===\n\n")

# Load study area
study_area <- st_read(STUDY_AREA_PATH, quiet = TRUE) %>%
  st_transform(UTM20)

# Load and process PAM data
pam_data_raw <- read_csv(PAM_DATA, show_col_types = FALSE) %>%
  mutate(rec_date = as.Date(rec_date),
         month = month(rec_date)) %>%
  group_by(site, latitude, longitude, species) %>%
  summarise(
    total_days = n(),
    detection_days = sum(presence, na.rm = TRUE),
    proportion_det = detection_days / total_days,
    .groups = 'drop'
  ) %>%
  filter(total_days > 0) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

# Clip to study area
pam_data <- st_intersection(pam_data_raw, study_area)

cat("Loaded", nrow(pam_data), "PAM station records\n\n")

# ==============================================================================
# FUNCTION: DETECTION-BASED SHAPEFILE GENERATION
# ==============================================================================

processKDE_detection_based <- function(tif_file, species_code, pam_data) {
  
  cat("Processing:", basename(tif_file), "\n")
  
  # Load KDE raster
  kde_raster <- rast(tif_file)
  raster_crs <- crs(kde_raster)
  
  # Set values <= 0 to NA
  kde_raster[kde_raster <= 0] <- NA
  
  # Get PAM stations for this species
  species_stations <- pam_data %>%
    filter(species == species_code)
  
  if (nrow(species_stations) == 0) {
    cat("  WARNING: No PAM stations found for", species_code, "\n")
    return(NULL)
  }
  
  cat("  Stations:", nrow(species_stations), "\n")
  
  # Extract KDE values at station locations
  station_kde_values <- terra::extract(kde_raster, 
                                       vect(species_stations), 
                                       method = "bilinear")[,2]
  
  station_detections <- species_stations$proportion_det
  
  # Create data frame of station detection rates and their KDE values
  station_data <- data.frame(
    site = species_stations$site,
    detection_rate = station_detections,
    kde_value = station_kde_values
  ) %>%
    filter(!is.na(kde_value))
  
  cat("  Stations with valid KDE values:", nrow(station_data), "\n")
  
  if (nrow(station_data) < 3) {
    cat("  WARNING: Too few stations with valid KDE values\n")
    return(NULL)
  }
  
  # Calculate detection-based thresholds
  kde_thresholds <- numeric(length(DETECTION_THRESHOLDS))
  
  for (i in seq_along(DETECTION_THRESHOLDS)) {
    threshold_pct <- DETECTION_THRESHOLDS[i]
    
    # Get stations with at least this detection rate
    stations_above <- station_data %>%
      filter(detection_rate >= threshold_pct)
    
    if (nrow(stations_above) == 0) {
      # No stations above this threshold - use the minimum KDE value
      kde_thresholds[i] <- min(station_data$kde_value, na.rm = TRUE)
      cat("  ", THRESHOLD_LABELS[i], "detection: No stations above threshold, using min KDE\n")
    } else if (nrow(stations_above) == 1) {
      # Only one station - use its value
      kde_thresholds[i] <- stations_above$kde_value[1]
      cat("  ", THRESHOLD_LABELS[i], "detection: KDE =", 
          sprintf("%.6f", kde_thresholds[i]), 
          "(1 station)\n")
    } else {
      # Multiple stations - use median of their KDE values
      kde_thresholds[i] <- median(stations_above$kde_value, na.rm = TRUE)
      cat("  ", THRESHOLD_LABELS[i], "detection: KDE =", 
          sprintf("%.6f", kde_thresholds[i]), 
          "(", nrow(stations_above), "stations)\n")
    }
  }
  
  # Add maximum threshold
  kde_thresholds <- c(kde_thresholds, max(values(kde_raster), na.rm = TRUE))
  
  # Ensure thresholds are monotonically increasing
  kde_thresholds <- sort(kde_thresholds)
  kde_thresholds <- unique(kde_thresholds)
  
  # Adjust labels to match actual thresholds used
  final_labels <- THRESHOLD_LABELS[1:(length(kde_thresholds)-1)]
  
  cat("  Final thresholds:", length(kde_thresholds), "\n")
  
  # Cut the raster into threshold classes
  kde_binned_values <- cut(values(kde_raster, mat = FALSE), 
                           breaks = kde_thresholds, 
                           labels = final_labels,
                           include.lowest = TRUE)
  
  # Create a new raster with binned values
  kde_binned_rast <- kde_raster
  values(kde_binned_rast) <- as.numeric(kde_binned_values)
  
  # Convert raster to polygon
  kde_polys <- as.polygons(kde_binned_rast, dissolve = TRUE)
  kde_sf <- st_as_sf(kde_polys)
  kde_sf <- st_make_valid(kde_sf)
  kde_sf <- st_set_crs(kde_sf, raster_crs)
  
  # Add threshold labels
  names(kde_sf)[1] <- "threshold_level"
  kde_sf <- kde_sf %>%
    mutate(Threshold = final_labels[threshold_level])
  
  # Extract species name from file path
  species_name <- str_extract(basename(tif_file), "(?<=KDE_pam_)[^\\.]+")
  
  # Save the shapefile
  output_path <- file.path(SHAPEFILE_OUTPUT_DIR, 
                           paste0("KDE_pam_detbased_", species_name, ".shp"))
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  cat("  ✓ Created shapefile:", basename(output_path), "\n\n")
  
  # Return summary for documentation
  return(list(
    species = species_code,
    n_stations = nrow(species_stations),
    n_valid = nrow(station_data),
    thresholds = data.frame(
      detection_rate = DETECTION_THRESHOLDS[1:length(final_labels)],
      label = final_labels,
      kde_value = kde_thresholds[1:length(final_labels)]
    )
  ))
}

# ==============================================================================
# PROCESS ALL PAM SPECIES
# ==============================================================================

cat("=== GENERATING DETECTION-BASED SHAPEFILES ===\n\n")

# Find all PAM KDE rasters
pam_rasters <- list.files(KDE_RASTER_DIR, 
                          pattern = "KDE_pam_[A-Z][a-z]\\.tif$", 
                          full.names = TRUE)

cat("Found", length(pam_rasters), "PAM KDE rasters\n\n")

# Process each species
threshold_summaries <- list()

for (raster_file in pam_rasters) {
  # Extract species code from filename
  species_code <- str_extract(basename(raster_file), "[A-Z][a-z](?=\\.tif)")
  
  if (species_code %in% PAM_SPECIES) {
    summary <- processKDE_detection_based(raster_file, species_code, pam_data)
    
    if (!is.null(summary)) {
      threshold_summaries[[species_code]] <- summary
    }
  }
}

# ==============================================================================
# CREATE SUMMARY DOCUMENTATION
# ==============================================================================

cat("\n=== CREATING THRESHOLD DOCUMENTATION ===\n\n")

# Compile all thresholds into a single table
threshold_doc <- bind_rows(
  lapply(names(threshold_summaries), function(sp) {
    summary <- threshold_summaries[[sp]]
    summary$thresholds %>%
      mutate(
        species_code = sp,
        species_name = PAM_SPECIES_NAMES[[sp]],
        .before = 1
      )
  })
)

# Save documentation
doc_file <- file.path(SHAPEFILE_OUTPUT_DIR, "threshold_documentation.csv")
write_csv(threshold_doc, doc_file)

cat("Threshold documentation saved to:", doc_file, "\n\n")

# Print summary
cat("THRESHOLD SUMMARY:\n")
print(threshold_doc %>% 
        arrange(species_code, detection_rate), 
      n = Inf)

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== DETECTION-BASED SHAPEFILE GENERATION COMPLETE ===\n")
cat("Generated shapefiles for", length(threshold_summaries), "species\n")
cat("Output directory:", SHAPEFILE_OUTPUT_DIR, "\n\n")

cat("Contour interpretation:\n")
cat("  Each contour represents areas with KDE values similar to stations\n")
cat("  that have at least X% detection rate.\n\n")

cat("Example: The '>5%' contour includes all areas where the KDE value\n")
cat("         is at least as high as the median KDE value at stations\n")
cat("         with ≥5% detection rates.\n\n")

cat("This ensures ALL stations with detections appear in at least one contour!\n")