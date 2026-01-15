# ==============================================================================
# KDE COMPARISON: PAM vs SIGHTINGS-------
# ==============================================================================
# Purpose: Compare KDE outputs from PAM data vs sightings data for species
#          that have both data types
#
# Creates comparison maps showing:
# - Side-by-side PAM and sightings KDE quantiles
# - Original point/station data overlaid
# - Identifies spatial agreement/disagreement
# - INCLUDES GROUPED ALL BALEEN WHALES COMPARISON
# ==============================================================================

options(scipen = 999)

pacman::p_load(sf, tidyverse, terra, ggplot2, viridis, ggspatial, patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Projections
UTM20 <- st_crs(32620)

# Input paths
SHAPEFILE_DIR <- "output/shapes/"
SIGHTINGS_DATA <- "input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv"
PAM_DATA <- "input/2025/baleen_presence_laura_2025.csv"
LAND_PATH <- "/Users/chirp/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"
OSW_WIND_PATH <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"  # OSW wind area path

# Output directory
COMPARE_OUTPUT_DIR <- "output/FIGS/compare_sightings_pam/"
if (!dir.exists(COMPARE_OUTPUT_DIR)) dir.create(COMPARE_OUTPUT_DIR, recursive = TRUE)

# Color palettes
QUANTILE_PALETTE <- "mako"
POINT_COLOR <- "gold3"  # Changed from red to dark yellow for sightings
STATION_COLOR <- "orange"
OSW_COLOR <- "red"  # Color for OSW wind areas
OSW_ALPHA <- 0.3  # Transparency for wind areas

# Species mapping between PAM codes and common names
# Now uses simple base names - the smart finder will handle variations
SPECIES_MAPPING <- list(
  "Bb" = list(
    common_name = "WHALE-SEI",
    pam_name = "Sei Whale",
    sightings_base = "WHALE-SEI",  # Base name for smart finder
    pam_base = "Bb"
  ),
  "Bm" = list(
    common_name = "WHALE-BLUE",
    pam_name = "Blue Whale",
    sightings_base = "WHALE-BLUE",
    pam_base = "Bm"
  ),
  "Bp" = list(
    common_name = "WHALE-FIN",
    pam_name = "Fin Whale",
    sightings_base = "WHALE-FIN",
    pam_base = "Bp"
  ),
  "Mn" = list(
    common_name = "WHALE-HUMPBACK",
    pam_name = "Humpback Whale",
    sightings_base = "WHALE-HUMPBACK",
    pam_base = "Mn"
  ),
  "Ba" = list(
    common_name = "WHALE-MINKE",
    pam_name = "Minke Whale",
    sightings_base = "WHALE-MINKE",
    pam_base = "Ba"
  ),
  "Eg" = list(
    common_name = "WHALE-NORTH ATLANTIC RIGHT",
    pam_name = "North Atlantic Right Whale",
    sightings_base = "WHALE-NORTH_ATLANTIC_RIGHT",
    pam_base = "Eg"
  )
)

# Add grouped baleen whale mapping
GROUPED_BALEEN <- list(
  "ALL_BALEEN" = list(
    common_name = "ALL_BALEEN_MYSTICETES",
    pam_name = "All Baleen Whales",
    sightings_base = "Baleen_Whales",  # Flexible - will match various formats
    pam_base = "Baleen_Whales",
    is_grouped = TRUE
  )
)

# Also check beaked whales (if they have sightings data)
BEAKED_WHALES <- list(
  "Northern Bottlenose" = list(
    sightings_pattern = "KDE_combined_species_WHALE-NORTHERN_BOTTLENOSE",
    common_name = "WHALE-NORTHERN BOTTLENOSE"
  ),
  "Sowerby's Beaked" = list(
    sightings_pattern = "KDE_combined_species_WHALE-SOWERBY'S_BEAKED",
    common_name = "WHALE-SOWERBY'S BEAKED"
  )
)

# ==============================================================================
# LOAD BASE LAYERS
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

# Load OSW wind areas
cat("Loading OSW wind areas...\n")
osw_wind <- st_read(OSW_WIND_PATH, quiet = TRUE) %>%
  st_transform(UTM20)
cat("✓ Loaded", nrow(osw_wind), "OSW wind lease areas\n\n")

# ==============================================================================
# LOAD SOURCE DATA
# ==============================================================================

cat("Loading source data...\n")

# Load sightings data
sightings_data <- read_csv(SIGHTINGS_DATA, show_col_types = FALSE) %>%
  filter(!is.na(lon), !is.na(lat)) %>%
  filter(as.Date(date_utc) >= as.Date("2015-01-01")) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(UTM20)

cat("✓ Loaded", nrow(sightings_data), "sightings (2015+)\n")

# Debug: Check what species column exists and what values it has
cat("\nDEBUG: Sightings data columns:\n")
print(names(sightings_data))
cat("\nDEBUG: Unique species values (first 20):\n")
if ("common_name" %in% names(sightings_data)) {
  print(head(unique(sightings_data$common_name), 20))
}
cat("\n")

# Load PAM data
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

# Clip PAM to study area
pam_data <- st_intersection(pam_data_raw, study_area)

cat("✓ Loaded", nrow(pam_data), "PAM station records\n\n")

# ==============================================================================
# HELPER FUNCTION: SMART SHAPEFILE FINDER
# ==============================================================================

find_kde_shapefile <- function(base_name, shapefile_dir = SHAPEFILE_DIR, prefer_proximity = FALSE) {
  # Smart function to find KDE shapefiles with flexible pattern matching
  # Handles variations in naming (proximity, no proximity, different prefixes)
  
  # Build search patterns from least to most specific
  patterns <- c(
    paste0("KDE_.*", base_name),  # Most flexible - any KDE with this base name
    paste0("KDE_.*species.*", base_name),  # Species-specific
    paste0("KDE_.*combined.*", base_name),  # Combined datasets
    paste0("KDE_.*grouped.*", base_name),  # Grouped analyses
    paste0("KDE_.*mysticetes.*", base_name),  # Mysticetes alternative
    paste0("KDE_pam_", base_name)  # PAM data
  )
  
  # Search for matching files
  all_matches <- c()
  for (pattern in patterns) {
    matches <- list.files(shapefile_dir, 
                          pattern = paste0(pattern, "\\.shp$"),
                          full.names = FALSE,
                          ignore.case = TRUE)
    all_matches <- c(all_matches, matches)
  }
  
  # Remove duplicates
  all_matches <- unique(all_matches)
  
  if (length(all_matches) == 0) {
    return(NULL)
  }
  
  # If prefer_proximity is TRUE, prioritize files with "proximity" in name
  if (prefer_proximity && any(grepl("proximity", all_matches, ignore.case = TRUE))) {
    all_matches <- all_matches[grepl("proximity", all_matches, ignore.case = TRUE)]
  }
  
  # If multiple matches, prefer most recent or most specific
  if (length(all_matches) > 1) {
    # Prioritize by specificity (longer names are usually more specific)
    # But deprioritize ones with "old" or "backup" in name
    all_matches <- all_matches[!grepl("old|backup|test", all_matches, ignore.case = TRUE)]
    
    if (length(all_matches) > 1) {
      # Sort by file modification time (most recent first)
      full_paths <- file.path(shapefile_dir, all_matches)
      file_times <- file.info(full_paths)$mtime
      all_matches <- all_matches[order(file_times, decreasing = TRUE)]
    }
  }
  
  # Return the best match
  if (length(all_matches) > 0) {
    return(file.path(shapefile_dir, all_matches[1]))
  } else {
    return(NULL)
  }
}

# ==============================================================================
# HELPER FUNCTION: LOAD AND PREPARE KDE SHAPEFILE (UPDATED)
# ==============================================================================

load_kde_shapefile <- function(pattern_or_file, shapefile_dir = SHAPEFILE_DIR) {
  
  # Check if this is already a full file path
  if (file.exists(pattern_or_file) && grepl("\\.shp$", pattern_or_file)) {
    shp_file <- pattern_or_file
  } else {
    # It's a pattern - try exact match first
    shp_files <- list.files(shapefile_dir, 
                            pattern = paste0(pattern_or_file, "\\.shp$"), 
                            full.names = TRUE, 
                            ignore.case = TRUE)
    
    if (length(shp_files) == 0) {
      # No exact match found
      return(NULL)
    }
    
    # Take most recent if multiple matches
    if (length(shp_files) > 1) {
      file_times <- file.info(shp_files)$mtime
      shp_file <- shp_files[which.max(file_times)]
      cat("  Note: Multiple matches found, using most recent:", basename(shp_file), "\n")
    } else {
      shp_file <- shp_files[1]
    }
  }
  
  # Load shapefile
  kde_shp <- st_read(shp_file, quiet = TRUE) %>%
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
    cat("  WARNING: No quantile column found in", basename(shp_file), "\n")
    return(NULL)
  }
  
  # Rename to standard name
  kde_shp <- kde_shp %>%
    rename(Quantile = !!quantile_col)
  
  # Standardize quantile format
  unique_quantiles <- unique(kde_shp$Quantile)
  unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
  
  if (is.numeric(unique_quantiles)) {
    kde_shp <- kde_shp %>%
      mutate(Quantile = paste0(Quantile * 100, "%"))
  }
  
  # Convert to ordered factor
  quantile_levels <- unique(kde_shp$Quantile)
  quantile_nums <- as.numeric(gsub("%", "", quantile_levels))
  quantile_levels <- quantile_levels[order(quantile_nums)]
  
  kde_shp$Quantile <- factor(kde_shp$Quantile, 
                             levels = quantile_levels, 
                             ordered = TRUE)
  
  return(kde_shp)
}

# ==============================================================================
# FUNCTION: CREATE COMPARISON MAP (PAM vs SIGHTINGS)
# ==============================================================================

create_comparison_map <- function(species_code, species_info, 
                                  sightings_data, pam_data,
                                  is_grouped = FALSE) {
  
  species_name <- species_info$common_name
  cat("Creating comparison map for:", species_info$pam_name, "\n")
  
  # Use smart finder to locate KDE shapefiles
  cat("  Searching for sightings KDE with base name:", species_info$sightings_base, "\n")
  sightings_file <- find_kde_shapefile(species_info$sightings_base, prefer_proximity = TRUE)
  
  cat("  Searching for PAM KDE with base name:", species_info$pam_base, "\n")
  pam_file <- find_kde_shapefile(species_info$pam_base)
  
  # Load the shapefiles
  sightings_kde <- if (!is.null(sightings_file)) {
    cat("  Found sightings KDE:", basename(sightings_file), "\n")
    load_kde_shapefile(sightings_file)
  } else {
    cat("  No sightings KDE found\n")
    NULL
  }
  
  pam_kde <- if (!is.null(pam_file)) {
    cat("  Found PAM KDE:", basename(pam_file), "\n")
    load_kde_shapefile(pam_file)
  } else {
    cat("  No PAM KDE found\n")
    NULL
  }
  
  if (is.null(sightings_kde) && is.null(pam_kde)) {
    cat("  ERROR: No KDE shapefiles found for", species_info$pam_name, "\n\n")
    return(NULL)
  }
  
  # Filter source data - handle grouped vs individual species differently
  if (is_grouped) {
    # For grouped baleen: filter to all mysticetes
    cat("  Filtering all baleen whale (mysticete) records...\n")
    species_sightings <- sightings_data %>%
      filter(cetacean_family == "mysticete")
    
    # For PAM: get all baleen species
    species_pam <- pam_data  # Already only has baleen species
    
  } else {
    # For individual species: use normal filtering
    # Debug output
    cat("  Looking for species with common_name:", species_name, "\n")
    
    # Try exact match first
    species_sightings <- sightings_data %>%
      filter(common_name == species_name)
    
    # If no exact match, try various alternatives
    if (nrow(species_sightings) == 0) {
      cat("  No exact match, trying alternatives...\n")
      
      # Try with spaces instead of underscores
      species_name_alt1 <- gsub("_", " ", species_name)
      species_sightings <- sightings_data %>%
        filter(common_name == species_name_alt1)
      
      if (nrow(species_sightings) > 0) {
        cat("  Matched with spaces:", species_name_alt1, "\n")
      }
    }
    
    # Still no match? Try without hyphens
    if (nrow(species_sightings) == 0) {
      species_name_alt2 <- gsub("-", "", species_name)
      species_name_alt2 <- gsub("_", "", species_name_alt2)
      
      species_sightings <- sightings_data %>%
        filter(gsub("[-_]", "", common_name) == species_name_alt2)
      
      if (nrow(species_sightings) > 0) {
        cat("  Matched without hyphens/underscores\n")
      }
    }
    
    # Last resort: show what's available that's close
    if (nrow(species_sightings) == 0) {
      cat("  WARNING: Could not match species\n")
      cat("  Available species containing 'WHALE':\n")
      whale_species <- unique(sightings_data$common_name[grepl("WHALE", sightings_data$common_name)])
      print(head(whale_species, 10))
    }
    
    species_pam <- pam_data %>%
      filter(species == species_code)
  }
  
  cat("  Sightings:", nrow(species_sightings), "\n")
  cat("  PAM stations:", nrow(species_pam), "\n")
  
  # Get color palette (use the same for both maps)
  if (!is.null(sightings_kde)) {
    n_quantiles_s <- length(levels(sightings_kde$Quantile))
  } else {
    n_quantiles_s <- 0
  }
  
  if (!is.null(pam_kde)) {
    n_quantiles_p <- length(levels(pam_kde$Quantile))
  } else {
    n_quantiles_p <- 0
  }
  
  n_quantiles <- max(n_quantiles_s, n_quantiles_p)
  pal <- rev(get(QUANTILE_PALETTE)(n_quantiles + 1))
  
  # Base plot theme
  base_theme <- theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 8)
    )
  
  # Create sightings map
  if (!is.null(sightings_kde) && nrow(species_sightings) > 0) {
    # Get date range for sightings (years only)
    sight_dates <- as.Date(species_sightings$date_utc)
    sight_year_range <- paste(format(min(sight_dates), "%Y"), "to", 
                              format(max(sight_dates), "%Y"))
    
    p_sightings <- ggplot() +
      geom_sf(data = sightings_kde, aes(fill = Quantile), 
              color = NA, alpha = 0.7) +
      geom_sf(data = species_sightings, 
              color = POINT_COLOR, size = 1, alpha = 0.6) +
      geom_sf(data = land, fill = "grey60", color = NA) +
      geom_sf(data = osw_wind, fill = NA, color = OSW_COLOR, 
              alpha = OSW_ALPHA, linewidth = 1) +  # OSW wind areas on top
      geom_sf(data = study_area, fill = NA, color = "black", 
              linewidth = 0.8, linetype = "dashed") +
      scale_fill_manual(values = pal, 
                        name = "KDE Quantile",
                        drop = FALSE,
                        guide = "none") +  # Remove legend from sightings map
      coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
      labs(title = "Sightings",
           subtitle = paste0(nrow(species_sightings), " locations")) +
      base_theme +
      annotation_scale(location = "br", width_hint = 0.25)
  } else {
    sight_year_range <- "No data"
    p_sightings <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No sightings data", size = 6) +
      theme_void()
  }
  
  # Create PAM map
  if (!is.null(pam_kde) && nrow(species_pam) > 0) {
    # Separate zero detections from actual detections
    pam_zero <- species_pam %>% filter(proportion_det == 0)
    pam_detections <- species_pam %>% filter(proportion_det > 0)
    
    # Get date range from PAM data (need to access from pam_data_raw which has rec_date)
    # Filter the raw data for this species before it was grouped
    if (!is_grouped) {
      pam_raw_species <- read_csv(PAM_DATA, show_col_types = FALSE) %>%
        filter(species == species_code)
    } else {
      # For grouped, use all species
      pam_raw_species <- read_csv(PAM_DATA, show_col_types = FALSE)
    }
    
    if(nrow(pam_raw_species) > 0) {
      pam_dates <- as.Date(pam_raw_species$rec_date)
      pam_year_range <- paste(format(min(pam_dates, na.rm = TRUE), "%Y"), "to", 
                              format(max(pam_dates, na.rm = TRUE), "%Y"))
    } else {
      pam_year_range <- "No dates available"
    }
    
    # Create a dummy data frame for the "No detections" legend entry
    dummy_zero <- data.frame(
      x = NA,
      y = NA,
      label = "No Validated Detections"
    )
    
    p_pam <- ggplot() +
      geom_sf(data = pam_kde, aes(fill = Quantile), 
              color = NA, alpha = 0.7) +
      # Plot actual detections with combined size and color scale
      {if(nrow(pam_detections) > 0) geom_sf(data = pam_detections, 
                                            aes(size = proportion_det, color = proportion_det),
                                            alpha = 0.8)} +
      geom_sf(data = land, fill = "grey60", color = NA) +
      geom_sf(data = osw_wind, fill = NA, color = OSW_COLOR, 
              alpha = OSW_ALPHA, linewidth = 1) +  # OSW wind areas on top
      # Plot zero detections on top as simple gray points
      {if(nrow(pam_zero) > 0) geom_sf(data = pam_zero, 
                                      color = "gray20",
                                      fill = "gray80",
                                      size = 1.5,
                                      shape = 21)} +
      # Add dummy point for legend (will be invisible but create legend entry)
      geom_point(data = dummy_zero, aes(x = x, y = y, shape = label),
                 color = "gray20", fill = "gray80", size = 1.5) +
      scale_shape_manual(name = "", values = 21, 
                         labels = "No Validated Detections",
                         guide = guide_legend(order = 3, override.aes = list(size = 3))) +
      geom_sf(data = study_area, fill = NA, color = "black", 
              linewidth = 0.8, linetype = "dashed") +
      scale_fill_manual(values = pal, 
                        name = "KDE Quantile",
                        drop = FALSE,
                        guide = guide_legend(order = 1)) +
      scale_size_continuous(name = "Detection\nProportion",
                            range = c(1, 4),
                            limits = c(0, 1),
                            breaks = c(0.2, 0.5, 0.9),
                            guide = guide_legend(order = 2)) +
      scale_color_gradient(low = "orange", high = "darkorange", 
                           name = "Detection\nProportion",
                           limits = c(0, 1),
                           breaks = c(0.2, 0.5, 0.9),
                           guide = guide_legend(order = 2)) +
      coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
      labs(title = "PAM",
           subtitle = paste0(nrow(species_pam), " stations (", 
                             nrow(pam_zero), " with 0 detections)")) +
      base_theme +
      annotation_scale(location = "br", width_hint = 0.25)
  } else {
    pam_year_range <- "No data"
    p_pam <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = "No PAM data", size = 6) +
      theme_void()
  }
  
  # Combine plots side-by-side
  p_combined <- (p_sightings + p_pam) +
    plot_annotation(
      title = paste(species_info$pam_name, " (", sight_year_range, ")", sep = ""),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
  
  # Save
  output_file <- file.path(COMPARE_OUTPUT_DIR, 
                           paste0("compare_", 
                                  gsub(" ", "_", tolower(species_code)), 
                                  ".png"))
  ggsave(output_file, p_combined, width = 16, height = 8, dpi = 300, bg = "white")
  
  cat("  Saved:", output_file, "\n\n")
  
  return(p_combined)
}

# ==============================================================================
# RUN COMPARISONS FOR ALL SPECIES WITH BOTH DATA TYPES
# ==============================================================================

cat("=== CREATING PAM vs SIGHTINGS COMPARISON MAPS ===\n\n")

# Process baleen whales (have both PAM and sightings)
for (species_code in names(SPECIES_MAPPING)) {
  species_info <- SPECIES_MAPPING[[species_code]]
  
  create_comparison_map(species_code, species_info, 
                        sightings_data, pam_data,
                        is_grouped = FALSE)
}

# ==============================================================================
# CREATE GROUPED ALL BALEEN COMPARISON
# ==============================================================================

cat("=== CREATING GROUPED ALL BALEEN WHALES COMPARISON ===\n\n")

for (species_code in names(GROUPED_BALEEN)) {
  species_info <- GROUPED_BALEEN[[species_code]]
  
  create_comparison_map(species_code, species_info, 
                        sightings_data, pam_data,
                        is_grouped = TRUE)
}

# ==============================================================================
# DIAGNOSTIC OUTPUT FOR GROUPED BALEEN
# ==============================================================================

cat("\n=== GROUPED BALEEN KDE DIAGNOSTICS ===\n\n")

# Check if grouped baleen KDEs exist
for (species_code in names(GROUPED_BALEEN)) {
  species_info <- GROUPED_BALEEN[[species_code]]
  
  cat("Checking grouped baleen KDE outputs...\n")
  
  # Use smart finder
  sightings_file <- find_kde_shapefile(species_info$sightings_base)
  pam_file <- find_kde_shapefile(species_info$pam_base)
  
  # Load KDE shapefiles
  sightings_kde <- if (!is.null(sightings_file)) load_kde_shapefile(sightings_file) else NULL
  pam_kde <- if (!is.null(pam_file)) load_kde_shapefile(pam_file) else NULL
  
  if (!is.null(sightings_kde)) {
    cat("✓ Sightings KDE found:", basename(sightings_file), "\n")
    cat("  Total area (all quantiles):", round(sum(st_area(sightings_kde)) / 1e6, 2), "km²\n")
    cat("  Number of polygons:", nrow(sightings_kde), "\n")
    
    # Count sightings used
    all_baleen_sightings <- sightings_data %>%
      filter(cetacean_family == "mysticete")
    cat("  Sightings included:", nrow(all_baleen_sightings), "records\n")
    
    # List species breakdown
    species_breakdown <- all_baleen_sightings %>%
      st_drop_geometry() %>%
      count(common_name, sort = TRUE)
    cat("  Species breakdown:\n")
    print(species_breakdown, n = Inf)
    
  } else {
    cat("✗ Sightings KDE NOT FOUND\n")
    cat("  Searched for base name:", species_info$sightings_base, "\n")
  }
  
  if (!is.null(pam_kde)) {
    cat("\n✓ PAM KDE found:", basename(pam_file), "\n")
    cat("  Total area (all quantiles):", round(sum(st_area(pam_kde)) / 1e6, 2), "km²\n")
    cat("  Number of polygons:", nrow(pam_kde), "\n")
    cat("  PAM stations included:", nrow(pam_data), "deployments\n")
  } else {
    cat("\n✗ PAM KDE NOT FOUND\n")
    cat("  Searched for base name:", species_info$pam_base, "\n")
  }
}

cat("\n=== COMPARISON NOTE ===\n")
cat("If the grouped baleen sightings KDE appears smaller than individual species:\n")
cat("  1. Check if proximity weighting was applied (makes KDE more restrictive)\n")
cat("  2. Check if bw.diggle was used (may produce narrower bandwidth for grouped data)\n")
cat("  3. Individual species KDEs are normalized independently (each maxes at 1.0)\n")
cat("  4. Grouped KDE combines all species, so density is spread across larger area\n")
cat("  5. Consider running grouped baleen with:\n")
cat("     - USE_PROXIMITY_WEIGHTING_FAMILY = FALSE\n")
cat("     - USE_BW_DIGGLE_FAMILY = FALSE (use fixed bandwidth)\n")
cat("\n")

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== COMPARISON MAPPING COMPLETE ===\n")
cat("Comparison maps saved to:", COMPARE_OUTPUT_DIR, "\n")
cat("\nIndividual species comparison maps created\n")
cat("✓ Grouped All Baleen Whales comparison map created\n")
cat("\nReview these maps to assess:\n")
cat("  ✓ Spatial agreement between PAM and sightings\n")
cat("  ✓ Areas where methods detect different hotspots\n")
cat("  ✓ Potential complementarity of data sources\n")
cat("  ✓ Data coverage differences\n")
cat("\nNOTE: If grouped baleen KDE appears smaller than expected:\n")
cat("  - Check diagnostic output above for bandwidth and weighting settings\n")
cat("  - Individual species KDEs are normalized independently (each peaks at 1.0)\n")
cat("  - Grouped KDE spreads density across all baleen whale locations\n")
cat("  - Consider running with fixed bandwidth and no proximity weighting\n")