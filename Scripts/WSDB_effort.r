################################################################################
# Script: Opportunistic Cetacean Distribution Mapping
# Author: Laura Joan Feyrer
# Date Updated: 2025-01-30
# 
# Purpose: Generate gridded maps of cetacean sightings from opportunistic data
#          Creates 3 types of maps:
#          1) Coverage map (all cetacean records)
#          2) Target species maps (raw counts + effort-normalized)
#          3) Data confidence/consistency assessment
#
# Changes from previous version:
# - Streamlined code structure with better organization
# - Added comprehensive annotations throughout
# - Reduced redundancy in spatial overlay code
# - Simplified confidence calculation logic
# - Improved error handling and validation
# - Consolidated mapping functions
################################################################################

# ==============================================================================
# 1. SETUP: PACKAGES & CONFIGURATION
# ==============================================================================

#' Load Required Packages
#' Installs missing packages automatically
setup_packages <- function() {
  # Core packages needed
  pkgs <- c("dplyr", "readr", "lubridate", "ggplot2", "stringr", 
            "tidyr", "sf", "ggspatial", "ggpattern", "writexl", "png")
  
  # Optional performance packages
  if (exists("use_datatable") && use_datatable) pkgs <- c(pkgs, "data.table")
  if (exists("use_parallel") && use_parallel) pkgs <- c(pkgs, "future", "furrr")
  
  # Install missing packages
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) {
    message("Installing: ", paste(to_install, collapse = ", "))
    install.packages(to_install)
  }
  
  # Load all packages
  invisible(lapply(pkgs, library, character.only = TRUE))
}

# ---------------------------
# USER CONFIGURATION
# ---------------------------

# Performance options (use for large datasets)
use_datatable <- TRUE    # Faster aggregation with data.table
use_parallel <- TRUE     # Parallel processing for multiple targets

# Input/output paths
infile <- "input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv"
out_dir <- "output/Effort/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Grid and filtering parameters
grid_km <- 25                    # Grid cell size in kilometers
min_records_all <- 5             # Minimum records to show normalized maps
year_min <- 2010                 # Start year for data filter
year_max <- 2025                 # End year for data filter
bbox_lonlat <- NULL              # Optional: c(xmin, xmax, ymin, ymax)
log_transform_all <- TRUE        # Log-transform coverage map to reduce bias
weight_col <- NULL               # Column for record weighting (NULL = count each row as 1)

# Spatial reference files
owa_polygons_path <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"
survey_polygons_path <- "shapefiles/WEA_5km_buffer_combined/WEA_5km_buffer_combined.shp"
study_area_path <- "shapefiles/studyArea/ESS_study_area_simple.shp"
land_path <- "shapefiles/Canada/NE_10mLand.shp"

# Map styling
owa_color <- "#FF6B35"           # Color for offshore wind area boundaries
owa_alpha <- 0.3                 # Transparency for OWA fill
owa_linewidth <- 1               # Line width for OWA borders
map_width <- 10                  # Output width in inches
map_height <- 8                  # Output height in inches
map_dpi <- 300                   # Resolution

# Seasonal definition
spring_months <- c(3, 4, 5)      # March, April, May

# Target species/groups (regex patterns for matching)
targets <- list(
  NARW = "Eubalaena\\s+glacialis|RIGHT\\s+WHALE|N\\s*ATLANTIC\\s+RIGHT",
  Blue_whale = "Balaenoptera\\s+musculus|BLUE\\s+WHALE",
  Shelf_baleen = "Megaptera|HUMPBACK|acutorostrata|MINKE|physalus|FIN\\s+WHALE|borealis|SEI|musculus|BLUE|glacialis|RIGHT|RORQUAL|BALEEN",
  Harbour_porpoise = "Phocoena\\s+phocoena|HARBOU?R\\s+PORPOISE",
  Small_odontocetes = "DOLPHIN|Delphin|Stenella|Tursiops|Lagenorhynchus|Grampus|WHITE[-\\s]?SIDED|WHITE[-\\s]?BEAKED|PILOT",
  Deep_divers = "Physeter|SPERM\\s+WHALE|Mesoplodon|Hyperoodon|Ziphius|BEAKED\\s+WHALE|SOWERBY|CUVIER|BOTTLENOSE"
)

# Human-readable titles for maps (must match targets list order)
target_titles <- list(
  NARW = "North Atlantic Right Whale",
  Blue_whale = "Blue Whale",
  Shelf_baleen = "Baleen Whales",
  Harbour_porpoise = "Harbour Porpoise",
  Small_odontocetes = "Small Odontocetes",
  Deep_divers = "Deep Divers"
)

# Color palettes for each target (helps distinguish different analyses)
target_palettes <- list(
  NARW = "inferno",
  Blue_whale = "cividis",
  Shelf_baleen = "magma",
  Harbour_porpoise = "turbo",
  Small_odontocetes = "plasma",
  Deep_divers = "rocket"
)

# Species-level outputs
species_field <- "common_name"              # Field to use for species summaries
appendix_xlsx <- file.path(out_dir, "APPENDIX_species_summary.xlsx")  # Changed to XLSX
make_species_maps <- TRUE                   # Generate individual species maps?
top_n_species_maps <- 30                    # How many species to map
species_map_dir <- file.path(out_dir, "species_maps")

# Species map styling
species_map_bins <- c(0, 5, 25, 50, 100, Inf)  # Bin breaks for categorical color scale
species_map_labels <- c("1-5", "6-25", "26-50", "51-100", "100+")  # Bin labels

# Species map filters
exclude_zero_owa <- TRUE                    # Exclude species with 0 OWA sightings?
exclude_beaked_whales <- TRUE               # Exclude beaked whale species?

dir.create(species_map_dir, showWarnings = FALSE, recursive = TRUE)

# Initialize packages
setup_packages()

# ==============================================================================
# 2. HELPER FUNCTIONS
# ==============================================================================

#' Create filename-safe string
#' Removes special characters and replaces with underscores
nm_safe <- function(x) {
  stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
}

#' Calculate longest consecutive year streak
#' Used for temporal consistency metrics
longest_streak <- function(years_vec) {
  y <- sort(unique(years_vec[!is.na(years_vec)]))
  if (length(y) == 0) return(0L)
  grp <- cumsum(c(TRUE, diff(y) != 1))
  max(as.integer(table(grp)))
}

#' Standardized map theme
#' Returns consistent ggplot theme for all maps
theme_map <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10),
      plot.caption = element_text(size = 8, hjust = 0)
    )
}

#' Add spatial overlays to ggplot map
#' Adds land, OWA boundaries, and study area in correct order
#' @param use_survey_not_owa If TRUE, shows survey area instead of OWA polygons (for species maps)
add_spatial_layers <- function(p, land, osw_wind, study_area, survey_area = NULL, use_survey_not_owa = FALSE) {
  # Order matters: land first, then wind areas/survey area, then boundaries
  if (!is.null(land)) {
    p <- p + geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
  }
  
  # Show either OWA or survey area (not both)
  if (use_survey_not_owa) {
    # For species maps: show survey area boundary
    if (!is.null(survey_area)) {
      p <- p + geom_sf(data = survey_area, fill = NA, color = "blue", 
                       linetype = "solid", alpha = 0.7, 
                       linewidth = 0.6, inherit.aes = FALSE)
    }
  } else {
    # For other maps: show OWA polygons
    if (!is.null(osw_wind)) {
      p <- p + geom_sf(data = osw_wind, fill = NA, color = owa_color, 
                       linetype = "dashed", alpha = owa_alpha, 
                       linewidth = owa_linewidth, inherit.aes = FALSE)
    }
  }
  
  if (!is.null(study_area)) {
    p <- p + geom_sf(data = study_area, fill = NA, color = "black", 
                     linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
  }
  # Add scale bar
  p <- p + ggspatial::annotation_scale(location = "br", width_hint = 0.25)
  return(p)
}

# ==============================================================================
# 3. DATA LOADING & CLEANING
# ==============================================================================

message("\n=== LOADING DATA ===")
message("Reading: ", infile)

# Define required columns (selective reading for efficiency)
required_cols <- c("lon", "lat", "x_m", "y_m", "date_utc", 
                   "scientific_name", "common_name")
optional_cols <- c("is_duplicate", if (!is.null(weight_col)) weight_col else NULL)

# Read data
df <- readr::read_csv(
  infile, 
  col_select = all_of(c(required_cols, optional_cols)),
  show_col_types = FALSE
)

message("Initial rows: ", nrow(df))

# Clean and filter data in one pipeline
df <- df %>%
  # Convert data types
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    x_m = as.numeric(x_m),
    y_m = as.numeric(y_m),
    date_utc = lubridate::as_date(date_utc),
    year = lubridate::year(date_utc),
    month = lubridate::month(date_utc),
    common_name = as.character(common_name),
    scientific_name = as.character(scientific_name),
    is_duplicate = if ("is_duplicate" %in% names(.)) as.logical(is_duplicate) else FALSE
  ) %>%
  # Apply filters
  filter(
    !is.na(x_m), !is.na(y_m), !is.na(date_utc),
    year >= year_min, year <= year_max
  )

# Optional spatial bbox filter
if (!is.null(bbox_lonlat)) {
  df <- df %>%
    filter(
      !is.na(lon), !is.na(lat),
      between(lon, bbox_lonlat["xmin"], bbox_lonlat["xmax"]),
      between(lat, bbox_lonlat["ymin"], bbox_lonlat["ymax"])
    )
}

message("Rows after filtering: ", nrow(df))

# Add weighting column (default = 1 per record)
df <- df %>%
  mutate(
    w = if (!is.null(weight_col) && weight_col %in% names(.)) {
      as.numeric(.data[[weight_col]])
    } else {
      1
    },
    w = ifelse(is.na(w) | w <= 0, 1, w)
  )

# Create combined name field for target matching
df <- df %>%
  mutate(
    name_blob = paste0(scientific_name, " | ", common_name),
    is_spring = month %in% spring_months
  )

# Aggregate unidentified species for cleaner mapping
# 1. All seal species → "All Seals"
# 2. Whale (ns) AND Cetacean (ns) → "Unidentified Whale" (COMBINED)
# 3. Fin + Sei whales → "Fin/Sei Whale"

df <- df %>%
  mutate(
    common_name = case_when(
      # Aggregate all true seals (exclude sea lions/fur seals)
      grepl("SEAL", common_name, ignore.case = TRUE) & 
        !grepl("SEA LION|FUR SEAL", common_name, ignore.case = TRUE) ~ "All Seals",
      
      # Aggregate ONLY "Whale-SEI" and "Whale-FIN/SEI" (but NOT "Whale-FIN")
      grepl("WHALE[-\\s]*FIN[-/\\s]*SEI|WHALE[-\\s]*SEI(?!.*FIN)", common_name, ignore.case = TRUE, perl = TRUE) ~ "Fin/Sei Whale",
      
      # Aggregate ALL unidentified whales AND cetaceans into ONE category
      # Combines: Whale (ns), Cetacean (ns), Unidentified Whale, Unidentified Cetacean
      grepl("WHALE.*\\(NS\\)|\\(NS\\).*WHALE|CETACEAN.*\\(NS\\)|\\(NS\\).*CETACEAN|UNIDENTIFIED.*WHALE|WHALE.*UNIDENTIFIED|UNIDENTIFIED.*CETACEAN|CETACEAN.*UNIDENTIFIED|WHALE.*NOT\\s+SPECIFIED|CETACEAN.*NOT\\s+SPECIFIED", 
            common_name, ignore.case = TRUE) ~ "Unidentified Whale",
      
      # Keep everything else as is
      TRUE ~ common_name
    ),
    
    # Update scientific names for aggregated groups
    scientific_name = case_when(
      common_name == "All Seals" ~ "Phocidae spp.",
      common_name == "Fin/Sei Whale" ~ "Balaenoptera physalus/borealis",
      common_name == "Unidentified Whale" ~ "Cetacea spp.",
      TRUE ~ scientific_name
    ),
    
    # Update name_blob after aggregation
    name_blob = paste0(scientific_name, " | ", common_name)
  )

message("  Aggregated species groups:")
message("    - Seals → 'All Seals'")
message("    - Fin + Sei whales → 'Fin/Sei Whale'")
message("    - Whale (ns) + Cetacean (ns) → 'Unidentified Whale'")

# ==============================================================================
# 4. GRID CALCULATION
# ==============================================================================

message("\n=== CALCULATING GRID CELLS ===")
cell_m <- grid_km * 1000  # Convert km to meters

df <- df %>%
  mutate(
    gx = floor(x_m / cell_m),
    gy = floor(y_m / cell_m),
    cell_id = paste0(gx, "_", gy),
    x0 = gx * cell_m,
    y0 = gy * cell_m,
    xc = x0 + cell_m/2,
    yc = y0 + cell_m/2
  )

# ==============================================================================
# 5. LOAD SPATIAL LAYERS
# ==============================================================================

message("\n=== LOADING SPATIAL LAYERS ===")

#' Load spatial file with error handling
load_spatial <- function(path, name) {
  if (!is.null(path) && file.exists(path)) {
    message("  Loading ", name, "...")
    shp <- sf::st_read(path, quiet = TRUE)
    return(sf::st_transform(shp, crs = 32620))  # Transform to UTM Zone 20N
  }
  return(NULL)
}

osw_wind <- load_spatial(owa_polygons_path, "offshore wind areas")
survey_area <- load_spatial(survey_polygons_path, "survey area")
study_area <- load_spatial(study_area_path, "study area")
land <- load_spatial(land_path, "land polygons")

# ==============================================================================
# 6. SPATIAL TAGGING
# ==============================================================================

message("\n=== TAGGING RECORDS BY LOCATION ===")

# Convert to sf points (using projected coordinates)
pts_sf <- sf::st_as_sf(df, coords = c("x_m", "y_m"), crs = 32620, remove = FALSE)

# Tag points within survey area (5km buffer around OWA)
# NOTE: Data is already clipped to study area in the input file
# This tagging identifies which records fall within the survey area subset
if (!is.null(survey_area)) {
  survey_union <- sf::st_union(survey_area)
  pts_sf$in_survey <- lengths(sf::st_intersects(pts_sf, survey_union)) > 0
  message("  Tagged ", sum(pts_sf$in_survey), " records in survey area (5km buffer)")
} else if (!is.null(osw_wind)) {
  # Fallback to OWA polygons if survey area not available
  osw_union <- sf::st_union(osw_wind)
  pts_sf$in_survey <- lengths(sf::st_intersects(pts_sf, osw_union)) > 0
  message("  Tagged ", sum(pts_sf$in_survey), " records in OWA polygons (fallback)")
} else {
  pts_sf$in_survey <- FALSE
  message("  No survey area or OWA polygons - all records tagged as outside survey area")
}

# Return to regular dataframe
df <- pts_sf %>% sf::st_drop_geometry()

# Optional: Convert to data.table for faster aggregation
if (use_datatable) {
  message("\n=== Converting to data.table ===")
  dt <- data.table::as.data.table(df)
  data.table::setkey(dt, cell_id)
} else {
  dt <- df
}

# ==============================================================================
# 7. SET MAP EXTENT
# ==============================================================================

message("\n=== Calculating map extent ===")

# Use data extent (not study area) to avoid cropping
xlims <- range(dt$xc, na.rm = TRUE) + c(-cell_m/2, cell_m/2)
ylims <- range(dt$yc, na.rm = TRUE) + c(-cell_m/2, cell_m/2)
message(sprintf("  X: [%.0f, %.0f] | Y: [%.0f, %.0f]", 
                xlims[1], xlims[2], ylims[1], ylims[2]))

# Crop land to map extent for efficiency
if (!is.null(land)) {
  bbox_crop <- sf::st_bbox(
    c(xmin = xlims[1], xmax = xlims[2], ymin = ylims[1], ymax = ylims[2]), 
    crs = sf::st_crs(32620)
  )
  land <- sf::st_crop(land, bbox_crop)
}

# ==============================================================================
# 8. AGGREGATE DATA: ALL CETACEAN RECORDS
# ==============================================================================

message("\n=== AGGREGATING ALL CETACEAN RECORDS ===")

if (use_datatable && data.table::is.data.table(dt)) {
  # Fast aggregation with data.table
  cell_all <- dt[, .(
    n_all = sum(w, na.rm = TRUE),
    n_years = data.table::uniqueN(year, na.rm = TRUE)
  ), by = .(cell_id, gx, gy, xc, yc)]
  cell_all <- as_tibble(cell_all)
} else {
  # Standard dplyr aggregation
  cell_all <- dt %>%
    group_by(cell_id, gx, gy, xc, yc) %>%
    summarise(
      n_all = sum(w, na.rm = TRUE),
      n_years = n_distinct(year, na.rm = TRUE),
      .groups = "drop"
    )
}

# ==============================================================================
# 9. CALCULATE CONFIDENCE INDEX
# ==============================================================================

message("\n=== CALCULATING DATA CONFIDENCE ===")

# Calculate quantile thresholds (using cells with data only)
cells_with_data <- cell_all %>% filter(n_all > 0)

# Coverage thresholds (33rd, 67th percentiles)
n_all_q33 <- quantile(cells_with_data$n_all, 0.33, na.rm = TRUE)
n_all_q67 <- quantile(cells_with_data$n_all, 0.67, na.rm = TRUE)

# Temporal thresholds (33rd, 67th percentiles)
n_years_q33 <- quantile(cells_with_data$n_years, 0.33, na.rm = TRUE)
n_years_q67 <- quantile(cells_with_data$n_years, 0.67, na.rm = TRUE)

message(sprintf("  Coverage: Low < %.0f | Med %.0f-%.0f | High > %.0f", 
                n_all_q33, n_all_q33, n_all_q67, n_all_q67))
message(sprintf("  Temporal: Low < %.0f yrs | Med %.0f-%.0f yrs | High > %.0f yrs", 
                n_years_q33, n_years_q33, n_years_q67, n_years_q67))

# Assign confidence categories
cell_all <- cell_all %>%
  mutate(
    # Categorize coverage
    coverage_cat = case_when(
      n_all == 0 ~ "None",
      n_all < n_all_q33 ~ "Low",
      n_all < n_all_q67 ~ "Medium",
      TRUE ~ "High"
    ),
    # Categorize temporal consistency
    temporal_cat = case_when(
      n_years == 0 ~ "None",
      n_years < n_years_q33 ~ "Low",
      n_years < n_years_q67 ~ "Medium",
      TRUE ~ "High"
    ),
    # Combined confidence score
    confidence = case_when(
      n_all == 0 | n_years == 0 ~ "No Data",
      coverage_cat == "High" & temporal_cat == "High" ~ "High",
      coverage_cat == "High" | temporal_cat == "High" ~ "Medium-High",
      coverage_cat == "Medium" & temporal_cat == "Medium" ~ "Medium",
      coverage_cat == "Medium" | temporal_cat == "Medium" ~ "Low-Medium",
      TRUE ~ "Low"
    ),
    confidence = factor(confidence, levels = c("High", "Medium-High", "Medium", 
                                               "Low-Medium", "Low"))
  )

# ==============================================================================
# 10. CREATE COVERAGE MAP (ALL RECORDS)
# ==============================================================================

message("\n=== CREATING COVERAGE MAP ===")

# Apply log transform if requested (reduces bias from high-density areas)
if (log_transform_all) {
  cell_all <- cell_all %>%
    mutate(n_all_display = log10(n_all + 1))
  
  # Create intuitive breaks at powers of 10
  log_breaks <- seq(0, ceiling(max(cell_all$n_all_display, na.rm = TRUE)), by = 1)
  actual_values <- 10^log_breaks
  actual_values[1] <- 0
  
  scale_label <- "log scale"
} else {
  cell_all <- cell_all %>%
    mutate(n_all_display = n_all)
  log_breaks <- NULL
  actual_values <- NULL
  scale_label <- ""
}

# Create plot
p_all <- ggplot(cell_all) +
  geom_tile(aes(x = xc, y = yc, fill = n_all_display)) +
  {if (log_transform_all) {
    scale_fill_viridis_c(
      option = "viridis", na.value = "grey90", name = "Count",
      breaks = log_breaks, labels = actual_values,
      begin = 0.15, end = 0.95
    )
  } else {
    scale_fill_viridis_c(
      option = "viridis", na.value = "grey90", name = "Count",
      begin = 0.15, end = 0.95
    )
  }} +
  labs(title = "Density of all cetacean sightings") +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid",
                          if(log_transform_all) " | log scale" else ""),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  theme_map()

# Add spatial layers
p_all <- add_spatial_layers(p_all, land, osw_wind, study_area, survey_area)

# Save
ggsave(
  file.path(out_dir, "00_all_cetacean_records.png"),
  p_all, width = map_width, height = map_height, dpi = map_dpi
)

message("  ✓ Saved: 00_all_cetacean_records.png")

# ==============================================================================
# 11. CREATE CONFIDENCE MAP
# ==============================================================================

message("\n=== CREATING CONFIDENCE MAP ===")

# Define color scheme (grey to green = low to high confidence)
confidence_colors <- c(
  "Low" = "#fee5d9",
  "Low-Medium" = "#fdd49e", 
  "Medium" = "#fdbb84",
  "Medium-High" = "#a1d99b",
  "High" = "#31a354"
)

# Mark low-coverage cells for hatching pattern
cell_all_pattern <- cell_all %>%
  mutate(low_coverage = confidence %in% c("No Data", "Low", "Low-Medium"))

# Create plot with diagonal hatching for low-coverage areas
p_conf <- ggplot(cell_all_pattern) +
  ggpattern::geom_tile_pattern(
    aes(x = xc, y = yc, fill = confidence,
        pattern = ifelse(low_coverage, "stripe", "none")),
    pattern_fill = "white", pattern_color = "white", 
    width = cell_m, height = cell_m,
        pattern_density = 0.02, pattern_spacing = 0.015,
    pattern_angle = 45, pattern_alpha = 0.5
  ) +
  scale_fill_manual(values = confidence_colors, name = "Data coverage", drop = FALSE) +
  scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
  guides(fill = guide_legend(override.aes = list(
    pattern = c("none", "none", "none", "stripe", "stripe"),
    pattern_fill = "white", pattern_color = "white"
  ))) +
  labs(
    title = "Cetacean sighting coverage + temporal consistency",
    caption = sprintf(
      "Coverage: Low < %.0f, Med %.0f-%.0f, High > %.0f records | Temporal: Low < %.0f, Med %.0f-%.0f, High > %.0f years",
      n_all_q33, n_all_q33, n_all_q67, n_all_q67,
      n_years_q33, n_years_q33, n_years_q67, n_years_q67
    )
  ) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  theme_map()

# Add spatial layers
p_conf <- add_spatial_layers(p_conf, land, osw_wind, study_area, survey_area)

# Save
ggsave(
  file.path(out_dir, "00_data_consistency.png"),
  p_conf, width = map_width, height = map_height, dpi = map_dpi
)

message("  ✓ Saved: 00_data_consistency.png")

# ==============================================================================
# 12. TARGET SPECIES MAPPING FUNCTION
# ==============================================================================

#' Create maps for a single target species/group
#' 
#' Generates two complementary maps:
#' 1. Raw count map - Shows WHERE target species were seen (affected by observer effort)
#' 2. Effort-normalized map - Shows RELATIVE PREFERENCE (controls for effort)
#' 
#' The effort-normalized map calculates:
#'   share_target = (target sightings in cell) / (all cetacean sightings in cell)
#' 
#' This proportion answers: "When observers are in this area and see cetaceans,
#' what percentage are the target species?" 
#' 
#' Higher values suggest:
#' - Target species are preferentially using this habitat, OR
#' - Other cetacean species avoid this area
#' 
#' Lower values suggest:
#' - Other cetacean species dominate encounters, OR
#' - Target species avoid this area
#' 
#' White hatching on normalized maps flags cells where the proportion may be unreliable:
#' - Very few total cetacean sightings (unstable denominator)
#' - Zero target sightings (can't distinguish true absence from insufficient sampling)
#' 
#' @param nm Target name (key in targets list)
#' @param target_regex Compiled regex pattern for matching
#' @param dt_input Data (dataframe or data.table)
#' @param cell_all_conf Grid cells with confidence scores from overall coverage analysis
#' @param use_dt Whether to use data.table operations (faster for large datasets)
#' 
#' @return NULL (saves maps to disk as side effect)
create_target_maps <- function(nm, target_regex, dt_input, cell_all_conf, use_dt = FALSE) {
  
  message("\n=== Processing: ", nm, " ===")
  
  # --- AGGREGATE DATA ---
  # Check if dt_input is actually a data.table object
  if (use_dt && data.table::is.data.table(dt_input)) {
    # Fast data.table aggregation
    cell_sum <- dt_input[, .(
      n_all = sum(w, na.rm = TRUE),
      n_target = sum(w * stringr::str_detect(name_blob, target_regex), na.rm = TRUE)
    ), by = .(cell_id, gx, gy, x0, y0, xc, yc)
        ]
    
    # Get target-specific years
    target_years <- dt_input[
      stringr::str_detect(name_blob, target_regex),
      .(n_target_years = data.table::uniqueN(year, na.rm = TRUE)),
      by = .(cell_id)
    ]
    
    cell_sum <- as_tibble(cell_sum)
    target_years <- as_tibble(target_years)
    
  } else {
    # Standard dplyr aggregation (use this if dt_input is not a data.table)
    cell_sum <- dt_input %>%
      mutate(is_target = stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id, gx, gy, xc, yc) %>%
      summarise(
        n_all = sum(w, na.rm = TRUE),
        n_target = sum(w[is_target], na.rm = TRUE),
        .groups = "drop"
      )
    
    target_years <- dt_input %>%
      filter(stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(n_target_years = n_distinct(year, na.rm = TRUE), .groups = "drop")
  }
  
  # Merge years and calculate effort-normalized share
  cell_sum <- cell_sum %>%
    left_join(target_years, by = "cell_id") %>%
    mutate(
      n_target_years = tidyr::replace_na(n_target_years, 0),
      
      # EFFORT-NORMALIZED SHARE CALCULATION:
      # share_target = (target species count) / (all cetacean count) per grid cell
      # 
      # What this means:
      # - If a cell has 100 total cetacean sightings and 30 are humpback whales,
      #   share_target = 0.30 (30% of encounters were humpback)
      # - This controls for observer effort: cells with more surveys get more sightings,
      #   but the PROPORTION reveals habitat preference
      # 
      # Why mask cells with n_all < min_records_all?
      # - Small denominators create unstable proportions
      # - Example: 1 humpback out of 2 total = 50% (misleading!)
      # - Example: 10 humpback out of 200 total = 5% (more reliable)
      # 
      # NOTE: This is NOT "proportion of all target sightings" (which would sum to 1.0)
      #       It's "proportion of local encounters" (varies by cell, can all be high or all low)
      
      share_target = ifelse(
        n_all >= min_records_all & n_all > 0,  # Only calculate if enough data
        n_target / n_all,                       # Proportion of encounters that are target
        NA_real_                                # Mask unreliable cells
      )
    )
  
  # --- CALCULATE TARGET-SPECIFIC CONFIDENCE ---
  # This combines two dimensions of data quality:
  # 1. Overall dataset coverage (from cell_all_conf): Are there enough cetacean sightings?
  # 2. Target species evidence: Are there enough target species sightings?
  # 
  # The goal is to identify cells where the PROPORTION (share_target) is unreliable,
  # either because:
  # - Too few total cetacean sightings (small denominator = unstable proportion)
  # - Zero or very few target sightings (can't assess habitat preference)
  
  cell_sum <- cell_sum %>%
    # Categorize target species evidence strength
    mutate(
      evidence = case_when(
        n_target == 0 ~ "None",                              # Never seen
        n_target >= 5 & n_target_years >= 3 ~ "High",       # Many sightings, multiple years
        n_target >= 2 & n_target_years >= 2 ~ "Medium",     # Some sightings, some years
        TRUE ~ "Low"                                         # Few sightings or single year
      ),
      evidence = factor(evidence, levels = c("High", "Medium", "Low", "None"))
    ) %>%
    # Merge with overall dataset confidence
    left_join(
      cell_all_conf %>% select(cell_id, confidence),
      by = "cell_id"
    ) %>%
    # Create combined confidence score for the normalized share map
    # RELAXED CRITERIA: Only flag cells that are genuinely problematic for proportion estimates
    mutate(
      target_conf = case_when(
        # Flag if BOTH overall coverage is poor AND target evidence is weak
        confidence %in% c("No Data", "Low") & evidence %in% c("None", "Low") ~ "Low (both)",
        # Flag if zero target sightings (proportion is 0, but is that real absence or just chance?)
        evidence == "None" ~ "No evidence",
        # Otherwise, accept the proportion as reasonable
        TRUE ~ "Acceptable"
      ),
      target_conf = factor(target_conf, levels = c("Acceptable", "No evidence", "Low (both)")),
      # Flag cells for hatching: ONLY if severe data issues
      # This is more permissive than before - allows Medium coverage + Low evidence, etc.
      low_coverage = target_conf %in% c("No evidence", "Low (both)")
    )
  
  # Safe filename
  nm_file <- nm_safe(nm)
  
  # --- MAP 1: RAW COUNTS ---
  p1 <- ggplot(cell_sum) +
    geom_tile(aes(x = xc, y = yc, fill = n_target)) +
    scale_fill_viridis_c(
      option = target_palettes[[nm]], na.value = "grey90", 
      name = "Count", begin = 0.10, end = 0.95
    ) +
    labs(
      title = paste0("Opportunistic sightings: ", target_titles[[nm]]),
      caption = paste0(
        "Raw counts reflect both species distribution AND observer effort.\n",
        "Note: Color scales vary by group - use normalized maps for comparisons."
      )
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
    coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
    theme_map()
  
  p1 <- add_spatial_layers(p1, land, osw_wind, study_area, survey_area)
  
  ggsave(
    file.path(out_dir, paste0(nm_file, "_1_target_records.png")),
    p1, width = map_width, height = map_height, dpi = map_dpi
  )
  
  # --- MAP 2: NORMALIZED SHARE (with hatching for low confidence) ---
  # Get low-confidence cells for outline
  low_conf_cells <- cell_sum %>%
    filter(low_coverage)
  
  p2 <- ggplot(cell_sum) +
    ggpattern::geom_tile_pattern(
      aes(x = xc, y = yc, fill = share_target, 
          pattern = ifelse(low_coverage, "stripe", "none")),
      pattern_fill = "white", pattern_color = "white", 
      width = cell_m, height = cell_m,
      
      pattern_density = 0.03, pattern_spacing = 0.015,
      pattern_angle = 45, pattern_alpha = 1.0, pattern_linewidth = 0.5
    ) +
    scale_fill_viridis_c(
      option = "viridis", na.value = "grey90", 
      limits = c(0, 1), name = "Proportion of\nall sightings",
      direction = 1, begin = 0.1, end = 0.95
    ) +
    scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
    # Add white outline for low-coverage cells
    {if(nrow(low_conf_cells) > 0) {
      geom_tile(data = low_conf_cells, aes(x = xc, y = yc), 
                fill = NA, color = "white", linewidth = 0.6)
    }} +
    labs(
      title = paste0("Effort-normalized sightings: ", target_titles[[nm]]),
      caption = paste0(
        "Each cell shows: (target species sightings) / (all cetacean sightings) = proportion of encounters.\n",
        "This controls for observer effort - cells with more surveys naturally have more sightings.\n",
        "White hatching flags cells with unreliable proportions: either (1) very few total cetacean sightings (poor coverage), or (2) zero target species seen.\n",
        "Cells masked (grey) where total cetacean records < ", min_records_all, "."
      )
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
    coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
    theme_map()
  
  p2 <- add_spatial_layers(p2, land, osw_wind, study_area, survey_area)
  
  ggsave(
    file.path(out_dir, paste0(nm_file, "_2_share.png")),
    p2, width = map_width, height = map_height, dpi = map_dpi
  )
  
  message("  ✓ Saved 2 maps for ", nm)
  
  return(invisible(NULL))
}

# ==============================================================================
# 13. PROCESS ALL TARGETS
# ==============================================================================

message("\n=== CREATING TARGET MAPS ===")

# Pre-compile regex patterns for efficiency
target_patterns <- lapply(targets, function(x) {
  stringr::regex(x, ignore_case = TRUE)
})

if (use_parallel && length(targets) >= 3) {
  # Parallel processing for 3+ targets
  message("Using parallel processing (", future::availableCores() - 1, " workers)")
  future::plan(future::multisession, 
               workers = min(length(targets), future::availableCores() - 1))
  
  furrr::future_walk(names(targets), function(nm) {
    create_target_maps(nm, target_patterns[[nm]], dt, cell_all, use_dt = use_datatable)
  })
  
  future::plan(future::sequential)
  
} else {
  # Sequential processing
  for (nm in names(targets)) {
    create_target_maps(nm, target_patterns[[nm]], dt, cell_all, use_dt = use_datatable)
  }
}

# ==============================================================================
# 14. SPECIES-LEVEL APPENDIX TABLE
# ==============================================================================

message("\n=== CREATING SPECIES SUMMARY TABLE ===")

# Calculate year span
yr_min_data <- min(df$year, na.rm = TRUE)
yr_max_data <- max(df$year, na.rm = TRUE)
n_years_total <- (yr_max_data - yr_min_data + 1)

# Aggregate by species
if (use_datatable && data.table::is.data.table(dt)) {
  sp <- dt[, .(
    n_study_area = sum(w, na.rm = TRUE),
    n_survey_any = sum(w[in_survey == TRUE], na.rm = TRUE),
    n_spring_study = sum(w[is_spring == TRUE], na.rm = TRUE),
    n_spring_survey_any = sum(w[is_spring == TRUE & in_survey == TRUE], na.rm = TRUE),
    years_study_area = data.table::uniqueN(year, na.rm = TRUE),
    years_survey_any = data.table::uniqueN(year[in_survey == TRUE], na.rm = TRUE),
    first_year_osw = suppressWarnings(min(year[in_survey == TRUE], na.rm = TRUE)),
    last_year_osw = suppressWarnings(max(year[in_survey == TRUE], na.rm = TRUE)),
    streak_osw_any = longest_streak(year[in_survey == TRUE])
  ), by = .(species = get(species_field))]
  sp <- as_tibble(sp)
} else {
  sp <- df %>%
    group_by(species = .data[[species_field]]) %>%
    summarise(
      n_study_area = sum(w, na.rm = TRUE),
      n_survey_any = sum(w[in_survey], na.rm = TRUE),
      n_spring_study = sum(w[is_spring], na.rm = TRUE),
      n_spring_survey_any = sum(w[is_spring & in_survey], na.rm = TRUE),
      years_study_area = n_distinct(year),
      years_survey_any = n_distinct(year[in_survey]),
      first_year_osw = suppressWarnings(min(year[in_survey], na.rm = TRUE)),
      last_year_osw = suppressWarnings(max(year[in_survey], na.rm = TRUE)),
      streak_osw_any = longest_streak(year[in_survey]),
      .groups = "drop"
    )
}

# Clean and add derived metrics
sp <- sp %>%
  mutate(
    first_year_osw = ifelse(is.infinite(first_year_osw), NA, first_year_osw),
    last_year_osw = ifelse(is.infinite(last_year_osw), NA, last_year_osw),
    
    # Proportions and percentages for context
    prop_years_survey_any = years_survey_any / n_years_total,
    pct_in_survey_area = ifelse(n_study_area > 0, 
                                round(100 * n_survey_any / n_study_area, 1), 
                                0),
    
    # Consistency based on survey area sightings (NOT entire study area)
    # This describes temporal consistency WITHIN the survey area
    consistency_survey_area = case_when(
      years_survey_any >= 5 ~ "frequent (5+ years in survey area)",
      years_survey_any >= 2 ~ "intermittent (2-4 years in survey area)",
      years_survey_any == 1 ~ "rare (1 year in survey area)",
      TRUE ~ "none in survey area"
    ),
    
    # Simple presence/absence flag for filtering
    present_in_survey_area = years_survey_any > 0
  ) %>%
  arrange(desc(n_survey_any), desc(n_study_area), species)

# Rename columns for clarity (survey area, not OWA)
sp_renamed <- sp %>%
  rename(
    n_survey_area = n_survey_any,
    n_spring_survey_area = n_spring_survey_any,
    years_survey_area_count = years_survey_any,
    first_year_survey_area = first_year_osw,
    last_year_survey_area = last_year_osw,
    streak_survey_area = streak_osw_any,
    prop_years_survey_area = prop_years_survey_any
  )

# Create metadata sheet
metadata <- tribble(
  ~Column, ~Description, ~Data_Type, ~Notes,
  "species", "Species common name (or 'All Seals' for aggregated seal species)", "text", "Based on species_field configuration",
  "n_study_area", "Total sightings across entire study area", "integer", "Includes all locations within study boundary",
  "n_survey_area", "Total sightings within survey area (5km buffer around OWA polygons)", "integer", "Subset of n_study_area",
  "n_spring_study", "Spring sightings in study area (March-May)", "integer", "Seasonal subset of n_study_area",
  "n_spring_survey_area", "Spring sightings in survey area (March-May)", "integer", "Seasonal subset of n_survey_area",
  "years_study_area", "Number of distinct years species was seen anywhere in study area", "integer", paste0("Out of ", n_years_total, " total years (", yr_min_data, "-", yr_max_data, ")"),
  "years_survey_area_count", "Number of distinct years species was seen in survey area", "integer", "Subset of years_study_area",
  "first_year_survey_area", "First year species was observed in survey area", "integer", "NA if never seen in survey area",
  "last_year_survey_area", "Most recent year species was observed in survey area", "integer", "NA if never seen in survey area",
  "streak_survey_area", "Longest consecutive year streak in survey area", "integer", "Maximum consecutive years with sightings",
  "prop_years_survey_area", "Proportion of total years where species was seen in survey area", "decimal", "years_survey_area_count / n_years_total",
  "pct_in_survey_area", "Percentage of total sightings that occurred in survey area", "decimal", "(n_survey_area / n_study_area) × 100. High values indicate survey area preference.",
  "consistency_survey_area", "Categorical label for temporal consistency in survey area", "text", "frequent (5+yrs) | intermittent (2-4yrs) | rare (1yr) | none",
  "present_in_survey_area", "Simple presence/absence flag for survey area", "logical", "TRUE if years_survey_area_count > 0"
)

# Add analysis context to metadata
analysis_info <- tribble(
  ~Parameter, ~Value,
  "Date generated", as.character(Sys.Date()),
  "Study area", basename(study_area_path),
  "Survey area", basename(survey_polygons_path),
  "OWA polygons", basename(owa_polygons_path),
  "Year range", paste0(year_min, "-", year_max),
  "Grid size (km)", as.character(grid_km),
  "Spring months", paste(month.name[spring_months], collapse = ", "),
  "Species aggregation", "Seals → 'All Seals' | Fin+Sei → 'Fin/Sei Whale' | Whale (ns) + Cetacean (ns) → 'Unidentified Whale'",
  "Species excluded from maps", "Beluga Whale, Northern Bottlenose Whale, Beaked whales",
  "Total species", as.character(nrow(sp_renamed)),
  "Species in survey area", as.character(sum(sp_renamed$present_in_survey_area)),
  "Data years available", paste0(yr_min_data, "-", yr_max_data, " (", n_years_total, " years)")
)

# Save as Excel with multiple sheets
writexl::write_xlsx(
  list(
    "Species_Data" = sp_renamed,
    "Column_Metadata" = metadata,
    "Analysis_Info" = analysis_info
  ),
  path = appendix_xlsx
)

message("  ✓ Saved: ", appendix_xlsx)
message("  Species richness: ", nrow(sp_renamed))
message("  Excel file contains 3 sheets: Species_Data, Column_Metadata, Analysis_Info")

# ==============================================================================
# 15. INDIVIDUAL SPECIES MAPS (OPTIONAL)
# ==============================================================================

if (make_species_maps) {
  
  message("\n=== CREATING INDIVIDUAL SPECIES MAPS ===")
  
  #' Create map for a single species
  create_species_map <- function(sp_name, dt_input, use_dt = FALSE) {
    sp_label <- as.character(sp_name)
    sp_file <- nm_safe(sp_label)
    
    # Aggregate by grid cell
    # Check if dt_input is actually a data.table object
    if (use_dt && data.table::is.data.table(dt_input)) {
      cell_sp <- dt_input[get(species_field) == sp_label, .(
        n_sp = sum(w, na.rm = TRUE)
      ), by = .(cell_id, gx, gy, xc, yc)]
      cell_sp <- as_tibble(cell_sp)
    } else {
      # Use dplyr regardless of use_dt setting if not a data.table
      cell_sp <- dt_input %>%
        filter(.data[[species_field]] == sp_label) %>%
        group_by(cell_id, gx, gy, xc, yc) %>%
        summarise(n_sp = sum(w, na.rm = TRUE), .groups = "drop")
    }
    # Also grab the raw point locations for this species (for hollow grey points)
    if (use_dt && data.table::is.data.table(dt_input)) {
      pts_sp <- dt_input[get(species_field) == sp_label, .(x_m, y_m)]
      pts_sp <- as.data.frame(pts_sp)
    } else {
      pts_sp <- dt_input %>%
        filter(.data[[species_field]] == sp_label) %>%
        select(x_m, y_m)
    }
    
    # Optional: drop rows with missing coords (just in case)
    pts_sp <- pts_sp %>% dplyr::filter(!is.na(x_m), !is.na(y_m))
    
    if (nrow(cell_sp) == 0) return(invisible(NULL))
    
    # Create standardized bins for consistent comparison across species
    # Uses configuration variables: species_map_bins and species_map_labels
    cell_sp <- cell_sp %>%
      mutate(
        count_category = cut(
          n_sp,
          breaks = species_map_bins,
          labels = species_map_labels,
          right = TRUE,
          include.lowest = FALSE
        ),
        count_category = factor(count_category, levels = species_map_labels)
      )
    
    # Define consistent colors for categories (light to dark, viridis-inspired)
    # Colors designed to work well for all species while being comparable
    category_colors <- c(
      "1-5" = "#FDE724",      # Light yellow (low counts)
      "6-25" = "#5DC863",     # Light green
      "26-50" = "#21908C",    # Teal (medium)
      "51-100" = "#3B528B",   # Dark blue
      "100+" = "#440154"      # Dark purple (high counts)
    )
    
    # Ensure color vector matches the labels
    names(category_colors) <- species_map_labels
    
    # Add a tiny dummy row for each missing category to force legend display
    # This ensures all categories show in the legend even if not present in data
    all_categories <- data.frame(
      cell_id = paste0("dummy_", species_map_labels),
      gx = rep(NA_real_, length(species_map_labels)),
      gy = rep(NA_real_, length(species_map_labels)),
      xc = rep(NA_real_, length(species_map_labels)),
      yc = rep(NA_real_, length(species_map_labels)),
      n_sp = rep(NA_real_, length(species_map_labels)),
      count_category = factor(species_map_labels, levels = species_map_labels)
    )
    
    # Combine with actual data (dummy rows won't plot because xc/yc are NA)
    cell_sp_with_dummies <- bind_rows(cell_sp, all_categories)
    
    # Create map with categorical bins
    p <- ggplot(cell_sp_with_dummies) +
      geom_tile(aes(x = xc, y = yc, fill = count_category),
                width = cell_m, height = cell_m) +
      scale_fill_manual(
        values = category_colors,
        name = "Sightings",
        na.value = "grey90",
        drop = FALSE  # Show all categories in legend even if not present
      ) +
      geom_point(
        data = pts_sp,
        aes(x = x_m, y = y_m),
        inherit.aes = FALSE,
        shape = 21,
        fill = NA,
        color = "grey40",
        alpha = 0.6,
        size = 1,
        stroke = 0.35
      ) +
      labs(
        title = paste0("Opportunistic sightings: ", sp_label),
        caption = "Standardized bins allow cross-species comparison. Raw counts reflect distribution AND effort."
      ) +
      annotate("text", x = -Inf, y = Inf,
               label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
               hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
      coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
      theme_map()
    
    # Use survey area boundary instead of OWA polygons for species maps
    p <- add_spatial_layers(p, land, osw_wind, study_area, survey_area, use_survey_not_owa = TRUE)
    
    ggsave(
      file.path(species_map_dir, paste0("species_", sp_file, "_raw.png")),
      p, width = map_width, height = map_height, dpi = map_dpi
    )
    
    invisible(NULL)
  }
  
  # Select top N species by OSW count, applying filters
  sp_filtered <- sp %>%
    filter(!is.na(species), species != "")
  
  # Apply OWA filter if enabled
  if (exclude_zero_owa) {
    sp_filtered <- sp_filtered %>%
      filter(n_survey_any > 0)
  }
  
  # Apply beaked whale filter if enabled
  if (exclude_beaked_whales) {
    sp_filtered <- sp_filtered %>%
      filter(!grepl("BEAKED|Mesoplodon|Hyperoodon|Ziphius", species, ignore.case = TRUE))
  }
  
  # Exclude Northern Bottlenose Whale (NBW) - single suspicious record
  sp_filtered <- sp_filtered %>%
    filter(!grepl("NORTHERN\\s+BOTTLENOSE|N\\s+BOTTLENOSE", species, ignore.case = TRUE))
  
  # Exclude Beluga Whale - not relevant for this analysis
  sp_filtered <- sp_filtered %>%
    filter(!grepl("BELUGA", species, ignore.case = TRUE))
  
  # Select top N
  sp_to_map <- sp_filtered %>%
    slice_head(n = top_n_species_maps) %>%
    pull(species)
  
  # Report what was filtered
  n_excluded <- nrow(sp) - nrow(sp_filtered)
  if (n_excluded > 0) {
    message("  Excluded ", n_excluded, " species (", 
            if(exclude_zero_owa) "zero OWA sightings" else "",
            if(exclude_zero_owa && exclude_beaked_whales) " + " else "",
            if(exclude_beaked_whales) "beaked whales" else "",
            ")")
  }
  message("  Creating maps for ", length(sp_to_map), " species")
  
  # Process (parallel if many species)
  if (use_parallel && length(sp_to_map) >= 8) {
    future::plan(future::multisession, 
                 workers = min(length(sp_to_map), future::availableCores() - 1))
    furrr::future_walk(sp_to_map, ~create_species_map(.x, dt, use_dt = use_datatable))
    future::plan(future::sequential)
  } else {
    for (spp in sp_to_map) {
      create_species_map(spp, dt, use_dt = use_datatable)
    }
  }
  
  message("  ✓ Species maps saved to: ", species_map_dir)
}

# Save all species maps as a multi-page PDF
species_pngs <- list.files(species_map_dir, pattern = "\\.png$", full.names = TRUE)
pdf(file.path(out_dir, "SPECIES_MAPS.pdf"), width = 11, height = 8.5)
for (i in seq_along(species_pngs)) {
  if (i > 1) grid::grid.newpage()  # New page before each image EXCEPT the first
  img <- png::readPNG(species_pngs[i])
  grid::grid.raster(img)
}
dev.off()
message("  ✓ Saved species maps PDF: ", file.path(out_dir, "SPECIES_MAPS.pdf"))

# ==============================================================================
# COMPLETION SUMMARY
# ==============================================================================

message("\n" , paste(rep("=", 80), collapse = ""))
message("ANALYSIS COMPLETE")
message(paste(rep("=", 80), collapse = ""))
message("\nOutput directory: ", out_dir)
message("\nMaps created:")
message("  • Coverage map: 00_all_cetacean_records.png")
message("  • Confidence map: 00_data_consistency.png")
message("  • Target maps: ", length(targets), " groups × 2 maps = ", length(targets) * 2, " maps")
if (make_species_maps) {
  message("  • Species maps: ", length(sp_to_map), " individual species")
}
message("\nData summary:")
message("  • Total records: ", nrow(df))
message("  • Years: ", year_min, "-", year_max)
message("  • Grid size: ", grid_km, " km")
message("  • Species richness: ", nrow(sp))
message("\nTarget groups: ", paste(names(targets), collapse = ", "))
message("\n", paste(rep("=", 80), collapse = ""))