# ==============================================================================
# KDE QUALITY ASSURANCE / QUALITY CHECK MAPPING--------
# ==============================================================================
# Purpose: Verify KDE outputs by overlaying quantile shapefiles with source data
#
# Creates QA maps showing:
# - KDE quantile polygons (from shapefiles)
# - Original point data (sightings or PAM stations)
# - Species-specific comparisons
# - Monthly comparisons (for PAM data)
#
# Use this to verify:
# - KDE hotspots align with high-density point clusters
# - Quantile boundaries make sense
# - No data loss or spatial misalignment
# ==============================================================================

options(scipen = 999)

pacman::p_load(sf, tidyverse, terra, ggplot2, viridis, ggspatial, patchwork)

# ==============================================================================
# CONFIGURATION-------
# ==============================================================================

# Projections
UTM20 <- st_crs(32620)

# Input paths
SHAPEFILE_DIR <- "output/shapes/"  # KDE quantile shapefiles
MONTHLY_SHAPEFILE_DIR <- "output/shapes/monthly/"  # Monthly KDE shapefiles
SIGHTINGS_DATA <- "input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv"
PAM_DATA <- "input/2025/baleen_presence_laura_2025.csv"
LAND_PATH <- "/Users/chirp/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"
OSW_WIND_PATH <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"  # OSW wind areas path

# Output directory
QA_OUTPUT_DIR <- "output/FIGS/QA/"
if (!dir.exists(QA_OUTPUT_DIR)) dir.create(QA_OUTPUT_DIR, recursive = TRUE)

# Color palettes
QUANTILE_PALETTE <- "mako"  # for KDE quantiles
POINT_COLOR <- "red"  # for sightings
STATION_COLOR <- "orange"  # for PAM stations

# Month names
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Species name mapping
species_names_pam <- list(
  "Bb" = "Sei Whale",
  "Bm" = "Blue Whale",
  "Bp" = "Fin Whale",
  "Mn" = "Humpback Whale",
  "Ba" = "Minke Whale",
  "Eg" = "North Atlantic Right Whale"
)

# ==============================================================================
# LOAD BASE LAYERS--------
# ==============================================================================

cat("Loading base layers...\n")

# Load land
land <- st_read(LAND_PATH, quiet = TRUE) %>%
  filter(CONTINENT == "North America") %>%
  st_transform(UTM20)

# Load study area
study_area <- st_read(STUDY_AREA_PATH, quiet = TRUE) %>%
  st_transform(UTM20)

# Set up bounding box
bbox <- st_bbox(study_area)
xlims <- c(bbox["xmin"], bbox["xmax"])
ylims <- c(bbox["ymin"], bbox["ymax"])

cat("✓ Base layers loaded\n\n")

# ==============================================================================
# FUNCTION: CREATE QA MAP FOR SIGHTINGS DATA---------
# ==============================================================================

create_sightings_qa_map <- function(shapefile_path, sightings_data, 
                                    species_cd, species_name) {
  
  cat("Creating QA map for:", species_name, "\n")
  
  # Load KDE shapefile
  if (!file.exists(shapefile_path)) {
    cat("  Shapefile not found:", shapefile_path, "\n")
    return(NULL)
  }
  
  kde_shp <- st_read(shapefile_path, quiet = TRUE) %>%
    st_transform(UTM20)
  
  # Debug: Check what columns exist
  cat("  Shapefile columns:", paste(names(kde_shp), collapse = ", "), "\n")
  
  # Try to find the quantile column (handle different possible names)
  quantile_col <- NULL
  possible_names <- c("Quantile", "Quantil", "quant_level", "quantile", "QUANTILE")
  
  for (col_name in possible_names) {
    if (col_name %in% names(kde_shp)) {
      quantile_col <- col_name
      break
    }
  }
  
  if (is.null(quantile_col)) {
    cat("  ERROR: No quantile column found. Available columns:", 
        paste(names(kde_shp), collapse = ", "), "\n")
    return(NULL)
  }
  
  cat("  Using quantile column:", quantile_col, "\n")
  
  # Rename to standard name for easier handling
  kde_shp <- kde_shp %>%
    rename(Quantile = !!quantile_col)
  
  # Filter sightings for this species
  species_sightings <- sightings_data %>%
    filter(species_cd == !!species_cd)
  
  if (nrow(species_sightings) == 0) {
    cat("  No sightings found for species\n")
    return(NULL)
  }
  
  cat("  Sightings:", nrow(species_sightings), "\n")
  
  # Get unique quantiles and create color palette
  unique_quantiles <- unique(kde_shp$Quantile)
  unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
  
  # Handle both numeric and character quantiles
  if (is.numeric(unique_quantiles)) {
    # If numeric, convert to percentage labels
    unique_quantiles <- sort(unique_quantiles)
    kde_shp <- kde_shp %>%
      mutate(Quantile = paste0(Quantile * 100, "%"))
    unique_quantiles <- paste0(unique_quantiles * 100, "%")
  } else {
    # If character, sort by extracting numbers
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    unique_quantiles <- unique_quantiles[order(quantile_nums)]
  }
  
  kde_shp$Quantile <- factor(kde_shp$Quantile, 
                             levels = unique_quantiles, 
                             ordered = TRUE)
  
  n_quantiles <- length(unique_quantiles)
  pal <- rev(get(QUANTILE_PALETTE)(n_quantiles + 1))
  
  # Create map
  p <- ggplot() +
    # KDE quantile polygons
    geom_sf(data = kde_shp, aes(fill = Quantile), 
            color = NA, alpha = 0.7) +
    # Point data overlay
    geom_sf(data = species_sightings, 
            color = POINT_COLOR, size = 0.8, alpha = 0.6) +
    # Land
    geom_sf(data = land, fill = "grey20", color = NA) +
    # Study area outline
    geom_sf(data = study_area, fill = NA, color = "black", 
            linewidth = 1, linetype = "dashed") +
    # Color scale
    scale_fill_manual(values = pal, 
                      name = "KDE Quantile",
                      drop = FALSE) +
    # Coordinate system
    coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
    # Labels
    labs(title = paste("QA Check:", species_name),
         subtitle = paste0("KDE Quantiles vs. ", nrow(species_sightings), 
                           " Sighting Locations")) +
    # Theme
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right"
    ) +
    # Scale bar
    annotation_scale(location = "br", width_hint = 0.25)
  
  # Save
  output_file <- file.path(QA_OUTPUT_DIR, 
                           paste0("QA_sightings_", species_cd, ".png"))
  ggsave(output_file, p, width = 10, height = 8, dpi = 300, bg = "white")
  
  cat("  Saved:", output_file, "\n\n")
  
  return(p)
}

# ==============================================================================
# FUNCTION: CREATE QA MAP FOR PAM DATA (OVERALL)--------
# ==============================================================================

create_pam_qa_map_overall <- function(shapefile_path, pam_data, 
                                      species_code, species_name) {
  
  cat("Creating overall PAM QA map for:", species_name, "\n")
  
  # Load KDE shapefile
  if (!file.exists(shapefile_path)) {
    cat("  Shapefile not found:", shapefile_path, "\n")
    return(NULL)
  }
  
  kde_shp <- st_read(shapefile_path, quiet = TRUE) %>%
    st_transform(UTM20)
  
  # Find the quantile column (handle different possible names)
  quantile_col <- NULL
  possible_names <- c("Quantile", "Quantil", "quant_level", "quantile", "QUANTILE")
  
  for (col_name in possible_names) {
    if (col_name %in% names(kde_shp)) {
      quantile_col <- col_name
      break
    }
  }
  
  if (is.null(quantile_col)) {
    cat("  ERROR: No quantile column found. Available columns:", 
        paste(names(kde_shp), collapse = ", "), "\n")
    return(NULL)
  }
  
  cat("  Using quantile column:", quantile_col, "\n")
  
  # Rename to standard name
  kde_shp <- kde_shp %>%
    rename(Quantile = !!quantile_col)
  
  # Filter PAM stations for this species
  species_pam <- pam_data %>%
    filter(species == species_code)
  
  if (nrow(species_pam) == 0) {
    cat("  No PAM data found for species\n")
    return(NULL)
  }
  
  cat("  Stations:", nrow(species_pam), "\n")
  
  # Get unique quantiles
  unique_quantiles <- unique(kde_shp$Quantile)
  unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
  
  # Handle both numeric and character quantiles
  if (is.numeric(unique_quantiles)) {
    # If numeric, convert to percentage labels
    unique_quantiles <- sort(unique_quantiles)
    kde_shp <- kde_shp %>%
      mutate(Quantile = paste0(Quantile * 100, "%"))
    unique_quantiles <- paste0(unique_quantiles * 100, "%")
  } else {
    # If character, sort by extracting numbers
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    unique_quantiles <- unique_quantiles[order(quantile_nums)]
  }
  
  kde_shp$Quantile <- factor(kde_shp$Quantile, 
                             levels = unique_quantiles, 
                             ordered = TRUE)
  
  n_quantiles <- length(unique_quantiles)
  pal <- rev(get(QUANTILE_PALETTE)(n_quantiles + 1))
  
  # Create map
  p <- ggplot() +
    # KDE quantile polygons
    geom_sf(data = kde_shp, aes(fill = Quantile), 
            color = NA, alpha = 0.7) +
    # PAM stations overlay (sized by detection proportion)
    geom_sf(data = species_pam, aes(size = proportion_det),
            color = STATION_COLOR, alpha = 0.7) +
    # Land
    geom_sf(data = land, fill = "grey20", color = NA) +
    # Study area outline
    geom_sf(data = study_area, fill = NA, color = "black", 
            linewidth = 1, linetype = "dashed") +
    # Color scale
    scale_fill_manual(values = pal, 
                      name = "KDE Quantile",
                      drop = FALSE) +
    scale_size_continuous(name = "Detection\nProportion",
                          range = c(1, 4),
                          limits = c(0, 1),
                          breaks = c(0.2, 0.5, 0.9)) +
    # Coordinate system
    coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
    # Labels
    labs(title = paste("QA Check (PAM):", species_name),
         subtitle = paste0("Overall KDE vs. ", nrow(species_pam), 
                           " Station Locations")) +
    # Theme
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5),
      legend.position = "right"
    ) +
    # Scale bar
    annotation_scale(location = "br", width_hint = 0.25)
  
  # Save
  output_file <- file.path(QA_OUTPUT_DIR, 
                           paste0("QA_pam_overall_", species_code, ".png"))
  ggsave(output_file, p, width = 10, height = 8, dpi = 300, bg = "white")
  
  cat("  Saved:", output_file, "\n\n")
  
  return(p)
}

# ==============================================================================
# FUNCTION: CREATE MONTHLY PAM QA FACET MAP----------
# ==============================================================================

create_pam_qa_map_monthly <- function(shapefile_dir, pam_data, 
                                      species_code, species_name) {
  
  cat("Creating monthly PAM QA map for:", species_name, "\n")
  
  # Find all monthly shapefiles for this species
  pattern <- paste0("KDE_pam_clipped_", species_code, "_month\\d+\\.shp$")
  shapefiles <- list.files(shapefile_dir, pattern = pattern, full.names = TRUE)
  
  if (length(shapefiles) == 0) {
    cat("  No monthly shapefiles found\n")
    return(NULL)
  }
  
  cat("  Found", length(shapefiles), "monthly shapefiles\n")
  
  # Read all monthly shapefiles
  monthly_kde <- lapply(shapefiles, function(f) {
    shp <- st_read(f, quiet = TRUE) %>% st_transform(UTM20)
    month_num <- as.numeric(str_extract(basename(f), "(?<=month)\\d+"))
    shp$month <- month_num
    shp$month_name <- month_names[month_num]
    return(shp)
  })
  
  monthly_kde <- do.call(rbind, monthly_kde)
  
  # Find the quantile column (handle different possible names)
  quantile_col <- NULL
  possible_names <- c("Quantile", "Quantil", "quant_level", "quantile", "QUANTILE")
  
  for (col_name in possible_names) {
    if (col_name %in% names(monthly_kde)) {
      quantile_col <- col_name
      break
    }
  }
  
  if (is.null(quantile_col)) {
    cat("  ERROR: No quantile column found. Available columns:", 
        paste(names(monthly_kde), collapse = ", "), "\n")
    return(NULL)
  }
  
  cat("  Using quantile column:", quantile_col, "\n")
  
  # Rename to standard name
  monthly_kde <- monthly_kde %>%
    rename(Quantile = !!quantile_col)
  
  # Filter PAM data for this species
  species_pam <- pam_data %>%
    filter(species == species_code) %>%
    mutate(month_name = month_names[month])
  
  # Convert month to factor for faceting
  monthly_kde$month_name <- factor(monthly_kde$month_name, 
                                   levels = month_names, ordered = TRUE)
  species_pam$month_name <- factor(species_pam$month_name, 
                                   levels = month_names, ordered = TRUE)
  
  # Get unique quantiles
  unique_quantiles <- unique(monthly_kde$Quantile)
  unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
  
  # Handle both numeric and character quantiles
  if (is.numeric(unique_quantiles)) {
    # If numeric, convert to percentage labels
    unique_quantiles <- sort(unique_quantiles)
    monthly_kde <- monthly_kde %>%
      mutate(Quantile = paste0(Quantile * 100, "%"))
    unique_quantiles <- paste0(unique_quantiles * 100, "%")
  } else {
    # If character, sort by extracting numbers
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    unique_quantiles <- unique_quantiles[order(quantile_nums)]
  }
  
  monthly_kde$Quantile <- factor(monthly_kde$Quantile, 
                                 levels = unique_quantiles, 
                                 ordered = TRUE)
  
  n_quantiles <- length(unique_quantiles)
  pal <- rev(get(QUANTILE_PALETTE)(n_quantiles + 1))
  
  # Create faceted map
  p <- ggplot() +
    # KDE quantile polygons
    geom_sf(data = monthly_kde, aes(fill = Quantile), 
            color = NA, alpha = 0.7) +
    # PAM stations overlay
    geom_sf(data = species_pam, aes(size = proportion_det),
            color = STATION_COLOR, alpha = 0.7) +
    # Land
    geom_sf(data = land, fill = "grey20", color = NA) +
    # Facet by month
    facet_wrap(~ month_name, ncol = 4) +
    # Color scale
    scale_fill_manual(values = pal, 
                      name = "KDE Quantile",
                      drop = FALSE) +
    scale_size_continuous(name = "Detection\nProportion",
                          range = c(0.5, 3),
                          limits = c(0, 1),
                          breaks = c(0.2, 0.5, 0.9)) +
    # Coordinate system
    coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
    # Labels
    labs(title = paste("Monthly QA Check (PAM):", species_name)) +
    # Theme
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 10, face = "bold"),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "right"
    )
  
  # Save
  output_file <- file.path(QA_OUTPUT_DIR, 
                           paste0("QA_pam_monthly_", species_code, ".png"))
  ggsave(output_file, p, width = 16, height = 12, dpi = 300, bg = "white")
  
  cat("  Saved:", output_file, "\n\n")
  
  return(p)
}

# ==============================================================================
# LOAD SOURCE DATA------
# ==============================================================================

cat("Loading source data...\n")

# Load sightings data
sightings_data <- read_csv(SIGHTINGS_DATA, show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(UTM20)

cat("✓ Loaded", nrow(sightings_data), "sightings\n")

# Load PAM data
pam_data_raw <- read_csv(PAM_DATA, show_col_types = FALSE) %>%
  mutate(rec_date = as.Date(rec_date),
         month = month(rec_date)) %>%
  group_by(site, latitude, longitude, species, month) %>%
  summarise(
    total_days = n(),
    detection_days = sum(presence, na.rm = TRUE),
    proportion_det = detection_days / total_days,
    .groups = 'drop'
  ) %>%
  filter(total_days > 0) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

# Clip PAM to study area
pam_data <- st_intersection(pam_data_raw, study_area)

cat("✓ Loaded", nrow(pam_data), "PAM station-month records\n\n")

# ==============================================================================
# RUN QA FOR SIGHTINGS (Example: Top species)------
# ==============================================================================

cat("=== CREATING SIGHTINGS QA MAPS ===\n\n")

# Get top sighting species
top_sightings_species <- sightings_data %>%
  st_drop_geometry() %>%
  filter(cetacean_family %in% c("odontocete", "mysticete")) %>%
  count(species_cd, common_name) %>%
  filter(n >= 50) %>%
  arrange(desc(n)) %>%
  slice_head(n = 5)  # Top 5 species

for (i in 1:nrow(top_sightings_species)) {
  sp_cd <- top_sightings_species$species_cd[i]
  sp_name <- top_sightings_species$common_name[i]
  
  # Find corresponding shapefile
  shp_pattern <- paste0("KDE_combined_species_", gsub(" ", "", sp_name), ".shp")
  shp_files <- list.files(SHAPEFILE_DIR, pattern = shp_pattern, 
                          full.names = TRUE, ignore.case = TRUE)
  
  if (length(shp_files) > 0) {
    create_sightings_qa_map(shp_files[1], sightings_data, sp_cd, sp_name)
  }
}

# ==============================================================================
# RUN QA FOR PAM DATA-----
# ==============================================================================

cat("=== CREATING PAM QA MAPS ===\n\n")

pam_species <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")

for (sp_code in pam_species) {
  sp_name <- species_names_pam[[sp_code]]
  
  # Overall KDE QA
  shp_overall <- file.path(SHAPEFILE_DIR, paste0("KDE_pam_clipped_", sp_code, ".shp"))
  if (file.exists(shp_overall)) {
    create_pam_qa_map_overall(shp_overall, pam_data, sp_code, sp_name)
  }
  
  # Monthly KDE QA
  create_pam_qa_map_monthly(MONTHLY_SHAPEFILE_DIR, pam_data, sp_code, sp_name)
}

# ==============================================================================
# SUMMARY------
# ==============================================================================

cat("\n=== QA MAPPING COMPLETE ===\n")
cat("QA maps saved to:", QA_OUTPUT_DIR, "\n")
cat("\nReview these maps to verify:\n")
cat("  ✓ KDE hotspots align with point clusters\n")
cat("  ✓ Quantile boundaries make sense\n")
cat("  ✓ No spatial misalignment\n")
cat("  ✓ No unexpected data gaps\n")