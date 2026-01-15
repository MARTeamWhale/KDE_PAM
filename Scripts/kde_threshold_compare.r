# ==============================================================================
# COMPLETE KDE THRESHOLD ANALYSIS WORKFLOW
# ==============================================================================
# Part 1: Process KDE rasters with absolute density threshold
# Part 2: Create side-by-side comparison maps (standard vs threshold)
#
# This script does everything in one go
# ==============================================================================

options(scipen = 999)

pacman::p_load(sf, stars, dplyr, stringr, terra, classInt, tidyverse, 
               ggplot2, viridis, ggspatial, patchwork, lubridate)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input/Output paths for threshold processing
INPUT_RASTERS <- "output/tif/"  # Overall KDE rasters
OUTPUT_SHAPEFILES <- "output/shapes/percentile_threshold/"

# Threshold settings - separate for PAM vs Sightings due to different data characteristics
PAM_PERCENTILE_THRESHOLD <- 0.25  # PAM: Remove bottom 25% (less aggressive, PAM already spatially constrained)
SIGHTINGS_PERCENTILE_THRESHOLD <- 0.50  # Sightings: Remove bottom 50% (more aggressive to reduce edge effects)
QUANTILE_PROBS <- c(0.85, 0.9, 0.95, 1)  # Skip 50%, focus on hotspots

# Paths for comparison
SHAPEFILE_DIR_STANDARD <- "output/shapes/"
SIGHTINGS_DATA <- "input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv"
PAM_DATA <- "input/2025/baleen_presence_laura_2025.csv"
LAND_PATH <- "/Users/chirp/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"
WEA_PATH <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"

QA_OUTPUT_DIR <- "output/FIGS/QA_threshold_comparison/"
if (!dir.exists(QA_OUTPUT_DIR)) dir.create(QA_OUTPUT_DIR, recursive = TRUE)

# Visualization settings
QUANTILE_PALETTE <- "mako"
POINT_COLOR <- "orange"
STATION_COLOR <- "orange"
WEA_COLOR <- "red"
UTM20 <- st_crs(32620)

month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

species_names_pam <- list(
  "Bb" = "Sei Whale", "Bm" = "Blue Whale", "Bp" = "Fin Whale",
  "Mn" = "Humpback Whale", "Ba" = "Minke Whale", "Eg" = "North Atlantic Right Whale"
)

# ==============================================================================
# PART 1: PROCESS KDE RASTERS WITH PERCENTILE THRESHOLD
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("PART 1: PROCESSING KDE RASTERS WITH PERCENTILE THRESHOLD\n")
cat("==============================================================================\n")
cat("Input rasters:", INPUT_RASTERS, "\n")
cat("Output shapefiles:", OUTPUT_SHAPEFILES, "\n")
cat("PAM percentile threshold:", PAM_PERCENTILE_THRESHOLD * 100, "th percentile\n")
cat("Sightings percentile threshold:", SIGHTINGS_PERCENTILE_THRESHOLD * 100, "th percentile\n")
cat("Quantile breaks:", paste(QUANTILE_PROBS * 100, collapse = "%, "), "%\n")
cat("==============================================================================\n\n")

processKDE_PercentileThreshold <- function(tif_file, 
                                           output_dir = OUTPUT_SHAPEFILES,
                                           quantile_probs = QUANTILE_PROBS) {
  
  cat("Processing:", basename(tif_file), "\n")
  
  # Determine if this is PAM or sightings based on filename
  is_pam <- grepl("pam", basename(tif_file), ignore.case = TRUE)
  percentile_threshold <- ifelse(is_pam, PAM_PERCENTILE_THRESHOLD, SIGHTINGS_PERCENTILE_THRESHOLD)
  
  cat("  Data type:", ifelse(is_pam, "PAM", "Sightings"), "\n")
  cat("  Using percentile threshold:", percentile_threshold * 100, "%\n")
  
  kde_raster <- rast(tif_file)
  raster_crs <- crs(kde_raster)  
  kde_values <- values(kde_raster, mat = FALSE)
  
  original_non_na <- sum(!is.na(kde_values))
  
  # Calculate percentile threshold
  percentile_cutoff <- quantile(kde_values, probs = percentile_threshold, na.rm = TRUE)
  cat("  Percentile cutoff value:", round(percentile_cutoff, 6), "\n")
  
  # Apply percentile threshold - remove bottom X% of density values
  kde_raster[kde_raster < percentile_cutoff] <- NA
  
  kde_values_filtered <- values(kde_raster, mat = FALSE)
  cells_remaining <- sum(!is.na(kde_values_filtered))
  
  cat("  Kept", cells_remaining, "cells (", 
      round(100 * cells_remaining / original_non_na, 1), "%)\n")
  
  if (cells_remaining == 0) {
    cat("  WARNING: No cells remain after filtering. Skipping...\n\n")
    return(invisible(NULL))
  }
  
  # Calculate quantiles on filtered data
  quantiles <- quantile(kde_values_filtered, probs = quantile_probs, 
                        type = 6, na.rm = TRUE)
  
  quant_labels <- paste0(quantile_probs[-length(quantile_probs)] * 100, "%")
  
  kde_binned_values <- cut(kde_values_filtered, breaks = quantiles, 
                           labels = quant_labels, include.lowest = TRUE)
  
  kde_quant_rast <- kde_raster
  values(kde_quant_rast) <- as.numeric(kde_binned_values)
  
  # Convert to polygon
  kde_polys <- as.polygons(kde_quant_rast, dissolve = TRUE)
  kde_sf <- st_as_sf(kde_polys)
  kde_sf <- st_make_valid(kde_sf)
  kde_sf <- st_set_crs(kde_sf, raster_crs)
  
  names(kde_sf)[1] <- "quant_level"
  kde_sf <- kde_sf %>%
    mutate(Quantile = quant_labels[quant_level])
  
  # Save shapefile - extract species name from filename
  # Try sightings pattern first: KDE_combined_species_SPECIESNAME.tif
  species_name <- str_extract(basename(tif_file), "(?<=KDE_combined_species_)[^\\.]+")
  
  # If that failed, try PAM pattern: KDE_pam_SPECIESCODE.tif
  if (is.na(species_name)) {
    species_name <- str_extract(basename(tif_file), "(?<=KDE_pam_)[^\\.]+")
  }
  
  # If still NA, try combined families: KDE_combined_odontocetes or KDE_combined_mysticetes
  if (is.na(species_name)) {
    species_name <- str_extract(basename(tif_file), "(?<=KDE_combined_odontocetes_)[^\\.]+")
  }
  if (is.na(species_name)) {
    species_name <- str_extract(basename(tif_file), "(?<=KDE_combined_mysticetes_)[^\\.]+")
  }
  
  # If still NA, try generic pattern: KDE_ANYTHING.tif
  if (is.na(species_name)) {
    species_name <- str_extract(basename(tif_file), "(?<=KDE_)[^\\.]+")
  }
  
  # Last resort: use full filename without extension
  if (is.na(species_name)) {
    species_name <- tools::file_path_sans_ext(basename(tif_file))
    cat("  WARNING: Could not extract species name, using full filename\n")
  }
  
  cat("  Extracted species name:", species_name, "\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_path <- file.path(output_dir, 
                           paste0(species_name, 
                                  "_pct", round(percentile_threshold * 100), ".shp"))
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  cat("  ✓ Saved:", basename(output_path), "\n\n")
  
  return(invisible(kde_sf))
}

# Process all rasters
tif_files <- list.files(INPUT_RASTERS, pattern = "\\.tif$", full.names = TRUE)
cat("Found", length(tif_files), "raster files to process\n\n")

success_count <- 0
for (tif_file in tif_files) {
  tryCatch({
    result <- processKDE_PercentileThreshold(tif_file)
    if (!is.null(result)) success_count <- success_count + 1
  }, error = function(e) {
    cat("ERROR:", basename(tif_file), "-", conditionMessage(e), "\n\n")
  })
}

cat("Successfully processed:", success_count, "files\n")

# ==============================================================================
# PART 2: LOAD BASE LAYERS FOR COMPARISON MAPS
# ==============================================================================

cat("\n")
cat("==============================================================================\n")
cat("PART 2: CREATING COMPARISON MAPS\n")
cat("==============================================================================\n\n")

cat("Loading base layers...\n")

land <- st_read(LAND_PATH, quiet = TRUE) %>%
  filter(CONTINENT == "North America") %>%
  st_transform(UTM20)

study_area <- st_read(STUDY_AREA_PATH, quiet = TRUE) %>%
  st_transform(UTM20)

wea <- st_read(WEA_PATH, quiet = TRUE) %>%
  st_transform(UTM20)

bbox <- st_bbox(study_area)
xlims <- c(bbox["xmin"], bbox["xmax"])
ylims <- c(bbox["ymin"], bbox["ymax"])

cat("✓ Base layers loaded\n\n")

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

find_quantile_column <- function(shp) {
  possible_names <- c("Quantile", "Quantil", "quant_level", "quantile", "QUANTILE")
  for (col_name in possible_names) {
    if (col_name %in% names(shp)) return(col_name)
  }
  return(NULL)
}

create_single_kde_map <- function(kde_shp, point_data, map_title, custom_palette = NULL) {
  quantile_col <- find_quantile_column(kde_shp)
  if (is.null(quantile_col)) return(NULL)
  
  kde_shp <- kde_shp %>% rename(Quantile = !!quantile_col)
  
  unique_quantiles <- unique(kde_shp$Quantile)
  unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
  
  if (is.numeric(unique_quantiles)) {
    unique_quantiles <- sort(unique_quantiles)
    kde_shp <- kde_shp %>% mutate(Quantile = paste0(Quantile * 100, "%"))
    unique_quantiles <- paste0(unique_quantiles * 100, "%")
  } else {
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    unique_quantiles <- unique_quantiles[order(quantile_nums)]
  }
  
  kde_shp$Quantile <- factor(kde_shp$Quantile, levels = unique_quantiles, ordered = TRUE)
  
  # Use custom palette if provided, otherwise generate one
  if (is.null(custom_palette)) {
    n_quantiles <- length(unique_quantiles)
    pal <- rev(get(QUANTILE_PALETTE)(n_quantiles + 1))
  } else {
    pal <- custom_palette
  }
  
  p <- ggplot() +
    geom_sf(data = kde_shp, aes(fill = Quantile), color = NA, alpha = 0.7) +
    geom_sf(data = point_data, color = POINT_COLOR, size = 0.8, alpha = 0.6) +
    geom_sf(data = wea, fill = NA, color = WEA_COLOR, linewidth = 1, linetype = "solid") +
    geom_sf(data = land, fill = "grey20", color = NA) +
    geom_sf(data = study_area, fill = NA, color = "black", 
            linewidth = 0.8, linetype = "dashed") +
    scale_fill_manual(values = pal, name = "KDE Quantile", drop = FALSE) +
    coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
    labs(title = map_title) +
    guides(fill = guide_legend(order = 1)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "right") +
    annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.5)
  
  return(p)
}

# ==============================================================================
# COMPARISON FUNCTIONS
# ==============================================================================

compare_sightings_kde <- function(species_cd, species_name, sightings_data) {
  cat("Comparing:", species_name, "\n")
  
  species_sightings <- sightings_data %>% filter(species_cd == !!species_cd)
  if (nrow(species_sightings) == 0) {
    cat("  No sightings\n\n")
    return(NULL)
  }
  
  # Remove spaces and underscores but keep hyphens and apostrophes
  species_name_clean <- toupper(gsub("[ _]", "", species_name))
  cat("  Searching for pattern:", species_name_clean, "\n")
  
  # List all files
  all_standard <- list.files(SHAPEFILE_DIR_STANDARD, pattern = "\\.shp$", full.names = TRUE)
  all_threshold <- list.files(OUTPUT_SHAPEFILES, pattern = "\\.shp$", full.names = TRUE)
  
  cat("  Standard directory has", length(all_standard), "files\n")
  cat("  Threshold directory has", length(all_threshold), "files\n")
  
  # Debug: show what we're looking for
  cat("  Looking for files containing:", species_name_clean, "\n")
  
  # Try to find matching files - search for species name anywhere in filename
  shp_standard <- NULL
  shp_threshold <- NULL
  
  for (f in all_standard) {
    fname_upper <- toupper(gsub("[ _]", "", basename(f)))  # Remove spaces AND underscores
    if (grepl(species_name_clean, fname_upper, fixed = TRUE)) {
      shp_standard <- f
      cat("  Found standard:", basename(f), "\n")
      break
    }
  }
  
  for (f in all_threshold) {
    fname_upper <- toupper(gsub("[ _]", "", basename(f)))  # Remove spaces AND underscores
    if (grepl(species_name_clean, fname_upper, fixed = TRUE)) {
      shp_threshold <- f
      cat("  Found threshold:", basename(f), "\n")
      break
    }
  }
  
  if (is.null(shp_standard)) {
    cat("  Standard shapefile not found\n")
    cat("  First 5 standard files:\n")
    for (i in 1:min(5, length(all_standard))) {
      cat("    -", basename(all_standard[i]), "\n")
    }
    cat("\n")
    return(NULL)
  }
  
  if (is.null(shp_threshold)) {
    cat("  Threshold shapefile not found\n")
    cat("  First 5 threshold files:\n")
    for (i in 1:min(5, length(all_threshold))) {
      cat("    -", basename(all_threshold[i]), "\n")
    }
    cat("\n")
    return(NULL)
  }
  
  kde_standard <- st_read(shp_standard, quiet = TRUE) %>% st_transform(UTM20)
  kde_threshold <- st_read(shp_threshold, quiet = TRUE) %>% st_transform(UTM20)
  
  # Get quantile columns
  qcol_std <- find_quantile_column(kde_standard)
  qcol_thr <- find_quantile_column(kde_threshold)
  if (is.null(qcol_std) || is.null(qcol_thr)) {
    cat("  Error: Could not find quantile columns\n\n")
    return(NULL)
  }
  
  # Debug: show what quantiles exist before filtering
  cat("  Standard quantiles before filter:", paste(unique(kde_standard[[qcol_std]]), collapse = ", "), "\n")
  
  # Filter out 50% quantile from standard for better visual comparison
  quants_std_raw <- kde_standard[[qcol_std]]
  if (is.numeric(quants_std_raw)) {
    # Numeric format (0.5, 0.85, 0.9, 0.95)
    kde_standard <- kde_standard %>% filter(!!sym(qcol_std) > 0.5)
  } else {
    # String format - check both with and without space
    kde_standard <- kde_standard %>% 
      filter(!(!!sym(qcol_std) %in% c("50%", "50 %", "0.5")))
  }
  
  # Debug: show what remains
  cat("  Standard quantiles after filter:", paste(unique(kde_standard[[qcol_std]]), collapse = ", "), "\n")
  
  if (nrow(kde_standard) == 0) {
    cat("  ERROR: All data filtered out! Skipping...\n\n")
    return(NULL)
  }
  
  # Get all unique quantiles from both datasets
  quants_std <- unique(kde_standard[[qcol_std]])
  quants_thr <- unique(kde_threshold[[qcol_thr]])
  
  # Convert to percentages if numeric
  if (is.numeric(quants_std)) {
    quants_std <- paste0(quants_std * 100, "%")
  }
  if (is.numeric(quants_thr)) {
    quants_thr <- paste0(quants_thr * 100, "%")
  }
  
  # Combine and sort all unique quantiles
  all_quantiles <- unique(c(quants_std, quants_thr))
  all_quantiles <- all_quantiles[!is.na(all_quantiles)]
  quantile_nums <- as.numeric(gsub("%", "", all_quantiles))
  all_quantiles <- all_quantiles[order(quantile_nums)]
  
  # Create a unified color palette for all quantiles
  n_colors <- length(all_quantiles)
  unified_palette <- rev(get(QUANTILE_PALETTE)(n_colors + 1))
  names(unified_palette) <- all_quantiles
  
  cat("  Using unified palette for quantiles:", paste(all_quantiles, collapse = ", "), "\n")
  
  p1 <- create_single_kde_map(kde_standard, species_sightings, "Standard Quantiles", unified_palette)
  p2 <- create_single_kde_map(kde_threshold, species_sightings, "Percentile Threshold", unified_palette)
  
  if (is.null(p1) || is.null(p2)) {
    cat("  Error creating maps\n\n")
    return(NULL)
  }
  
  p_combined <- p1 + p2 + 
    plot_annotation(
      title = paste("QA Comparison:", species_name),
      subtitle = paste0(nrow(species_sightings), " sighting locations"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 12, hjust = 0.5))
    )
  
  output_file <- file.path(QA_OUTPUT_DIR, paste0("compare_sightings_", species_cd, ".png"))
  ggsave(output_file, p_combined, width = 16, height = 8, dpi = 300, bg = "white")
  
  cat("  ✓ Saved:", basename(output_file), "\n\n")
  return(p_combined)
}

compare_pam_kde <- function(species_code, species_name, pam_data) {
  cat("Comparing PAM:", species_name, "\n")
  
  species_pam <- pam_data %>% filter(species == species_code)
  if (nrow(species_pam) == 0) {
    cat("  No PAM data\n\n")
    return(NULL)
  }
  
  # Split into zero and non-zero detections
  species_pam_zero <- species_pam %>% filter(proportion_det == 0)
  species_pam_detected <- species_pam %>% filter(proportion_det > 0)
  
  cat("  Stations with detections:", nrow(species_pam_detected), "\n")
  cat("  Stations with zero detections:", nrow(species_pam_zero), "\n")
  
  # For PAM, we need to search by species code (e.g., "Ba", "Bb", "Bp")
  cat("  Looking for PAM species code:", species_code, "\n")
  
  shp_standard <- file.path(SHAPEFILE_DIR_STANDARD, paste0("KDE_pam_", species_code, ".shp"))
  all_threshold <- list.files(OUTPUT_SHAPEFILES, pattern = "\\.shp$", full.names = TRUE)
  
  # Look for threshold files that START with the species code
  shp_threshold <- NULL
  for (f in all_threshold) {
    fname <- basename(f)
    # Match files like "Ba_pct25.shp" or "KDE_pam_Ba_pct25.shp"
    if (grepl(paste0("^", species_code, "_pct"), fname) || 
        grepl(paste0("pam_", species_code, "_pct"), fname)) {
      shp_threshold <- f
      cat("  Found threshold:", basename(f), "\n")
      break
    }
  }
  
  cat("  Looking for standard:", basename(shp_standard), "\n")
  cat("  Standard exists:", file.exists(shp_standard), "\n")
  
  if (!file.exists(shp_standard) || is.null(shp_threshold)) {
    cat("  Shapefiles not found\n")
    if (length(all_threshold) > 0) {
      cat("  Example threshold files:", paste(head(basename(all_threshold), 3), collapse = ", "), "\n")
    }
    cat("\n")
    return(NULL)
  }
  
  kde_standard <- st_read(shp_standard, quiet = TRUE) %>% st_transform(UTM20)
  kde_threshold <- st_read(shp_threshold, quiet = TRUE) %>% st_transform(UTM20)
  
  # Process quantiles
  qcol_std <- find_quantile_column(kde_standard)
  qcol_thr <- find_quantile_column(kde_threshold)
  if (is.null(qcol_std) || is.null(qcol_thr)) return(NULL)
  
  # Filter out 50% quantile from standard for better visual comparison
  cat("  Standard quantiles before filter:", paste(unique(kde_standard[[qcol_std]]), collapse = ", "), "\n")
  
  quants_std_raw <- kde_standard[[qcol_std]]
  if (is.numeric(quants_std_raw)) {
    kde_standard <- kde_standard %>% filter(!!sym(qcol_std) > 0.5)
  } else {
    kde_standard <- kde_standard %>% 
      filter(!(!!sym(qcol_std) %in% c("50%", "50 %", "0.5")))
  }
  
  cat("  Standard quantiles after filter:", paste(unique(kde_standard[[qcol_std]]), collapse = ", "), "\n")
  
  if (nrow(kde_standard) == 0) {
    cat("  ERROR: All data filtered out! Skipping...\n\n")
    return(NULL)
  }
  
  kde_standard <- kde_standard %>% rename(Quantile = !!qcol_std)
  kde_threshold <- kde_threshold %>% rename(Quantile = !!qcol_thr)
  
  uq_std <- sort(unique(kde_standard$Quantile[!is.na(kde_standard$Quantile)]))
  uq_thr <- sort(unique(kde_threshold$Quantile[!is.na(kde_threshold$Quantile)]))
  
  kde_standard$Quantile <- factor(kde_standard$Quantile, levels = uq_std, ordered = TRUE)
  kde_threshold$Quantile <- factor(kde_threshold$Quantile, levels = uq_thr, ordered = TRUE)
  
  pal_std <- rev(get(QUANTILE_PALETTE)(length(uq_std) + 1))
  pal_thr <- rev(get(QUANTILE_PALETTE)(length(uq_thr) + 1))
  
  p1 <- ggplot() +
    geom_sf(data = kde_standard, aes(fill = Quantile), color = NA, alpha = 0.7) +
    geom_sf(data = species_pam_zero, size = 1.5, color = "gray50", alpha = 0.7) +
    geom_sf(data = species_pam_detected, aes(size = proportion_det), color = STATION_COLOR, alpha = 0.7) +
    geom_sf(data = wea, fill = NA, color = WEA_COLOR, linewidth = 1, linetype = "solid") +
    geom_sf(data = land, fill = "grey20", color = NA) +
    geom_sf(data = study_area, fill = NA, color = "black", linewidth = 0.8, linetype = "dashed") +
    scale_fill_manual(values = pal_std, name = "KDE Quantile", drop = FALSE) +
    scale_size_continuous(name = "Detection Days/\nRecording Days",
                          range = c(1, 4), limits = c(0, 1), breaks = c(0.2, 0.5, 0.9)) +
    coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
    labs(title = "Standard Quantiles") +
    guides(fill = guide_legend(order = 1), size = guide_legend(order = 2)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "right") +
    annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.5)
  
  p2 <- ggplot() +
    geom_sf(data = kde_threshold, aes(fill = Quantile), color = NA, alpha = 0.7) +
    geom_sf(data = species_pam_zero, size = 1.5, color = "gray50", alpha = 0.7) +
    geom_sf(data = species_pam_detected, aes(size = proportion_det), color = STATION_COLOR, alpha = 0.7) +
    geom_sf(data = wea, fill = NA, color = WEA_COLOR, linewidth = 1, linetype = "solid") +
    geom_sf(data = land, fill = "grey20", color = NA) +
    geom_sf(data = study_area, fill = NA, color = "black", linewidth = 0.8, linetype = "dashed") +
    scale_fill_manual(values = pal_thr, name = "KDE Quantile", drop = FALSE) +
    scale_size_continuous(name = "Detection Days/\nRecording Days",
                          range = c(1, 4), limits = c(0, 1), breaks = c(0.2, 0.5, 0.9)) +
    coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
    labs(title = "Percentile Threshold") +
    guides(fill = guide_legend(order = 1), size = guide_legend(order = 2)) +
    theme_bw() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), axis.title = element_blank(),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          legend.position = "right") +
    annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.5)
  
  p_combined <- p1 + p2 + 
    plot_annotation(
      title = paste("PAM QA Comparison:", species_name),
      subtitle = paste0(nrow(species_pam), " station-month records"),
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                    plot.subtitle = element_text(size = 12, hjust = 0.5))
    )
  
  output_file <- file.path(QA_OUTPUT_DIR, paste0("compare_pam_", species_code, ".png"))
  ggsave(output_file, p_combined, width = 16, height = 8, dpi = 300, bg = "white")
  
  cat("  ✓ Saved:", basename(output_file), "\n\n")
  return(p_combined)
}

# ==============================================================================
# LOAD SOURCE DATA
# ==============================================================================

cat("Loading source data...\n")

sightings_data_raw <- read_csv(SIGHTINGS_DATA, show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat))

cat("Total sightings loaded:", nrow(sightings_data_raw), "\n")

# Filter to 2015 onwards (to match KDE generation)
sightings_data <- sightings_data_raw %>%
  filter(as.Date(date_utc) >= as.Date("2015-01-01")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(UTM20)

cat("✓ Loaded", nrow(sightings_data), "sightings (2015+)\n")

pam_data_raw <- read_csv(PAM_DATA, show_col_types = FALSE) %>%
  mutate(rec_date = as.Date(rec_date), month = month(rec_date)) %>%
  group_by(site, latitude, longitude, species, month) %>%
  summarise(total_days = n(), detection_days = sum(presence, na.rm = TRUE),
            proportion_det = detection_days / total_days, .groups = 'drop') %>%
  filter(total_days > 0) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

pam_data <- st_intersection(pam_data_raw, study_area)
cat("✓ Loaded", nrow(pam_data), "PAM station-month records\n\n")

# ==============================================================================
# RUN COMPARISONS
# ==============================================================================

cat("=== COMPARING SIGHTINGS ===\n\n")

top_sightings <- sightings_data %>%
  st_drop_geometry() %>%
  filter(cetacean_family %in% c("odontocete", "mysticete")) %>%
  count(species_cd, common_name) %>%
  filter(n >= 25) %>%  # Changed from 50 to 25
  arrange(desc(n))  # Removed slice_head(n = 10) to include all species

cat("Comparing", nrow(top_sightings), "species with ≥25 records\n\n")

for (i in 1:nrow(top_sightings)) {
  compare_sightings_kde(top_sightings$species_cd[i], 
                        top_sightings$common_name[i], 
                        sightings_data)
}

cat("=== COMPARING PAM ===\n\n")

for (sp_code in c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")) {
  compare_pam_kde(sp_code, species_names_pam[[sp_code]], pam_data)
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n==============================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================================================================\n")
cat("Threshold shapefiles:", OUTPUT_SHAPEFILES, "\n")
cat("Comparison maps:", QA_OUTPUT_DIR, "\n")
cat("\nReview side-by-side maps to see if percentile threshold approach reduces\n")
cat("isolated contours while preserving core hotspots.\n")
cat("==============================================================================\n")
