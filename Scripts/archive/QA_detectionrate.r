# ==============================================================================
# DETECTION-BASED KDE DIAGNOSTIC AND COMPARISON SCRIPT
# ==============================================================================
# Purpose: Evaluate detection-based KDE shapefiles and compare with quantile-based
# ==============================================================================

options(scipen = 999)
pacman::p_load(sf, tidyverse, terra, spatstat, patchwork)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

UTM20 <- st_crs(32620)

# Input paths
PAM_DATA <- "input/2025/baleen_presence_laura_2025.csv"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"
QUANTILE_SHAPEFILE_DIR <- "output/shapes/"
DETECTION_SHAPEFILE_DIR <- "output/shapes/detection_based/"
KDE_RASTER_DIR <- "output/tif/"

# Output
OUTPUT_DIR <- "output/diagnostics/detection_based/"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Species to check
SPECIES_CODE <- "Ba"  # Change to check other species
SPECIES_NAME <- "Minke Whale"

# Stations of interest
STATIONS_TO_CHECK <- c("CS1", "LOC")

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("=== LOADING DATA ===\n\n")

# Load study area
study_area <- st_read(STUDY_AREA_PATH, quiet = TRUE) %>%
  st_transform(UTM20)

# Load PAM data
pam_data_raw <- read_csv(PAM_DATA, show_col_types = FALSE) %>%
  mutate(rec_date = as.Date(rec_date)) %>%
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

# Filter for species
species_pam <- pam_data %>%
  filter(species == SPECIES_CODE)

cat("Total PAM stations for", SPECIES_NAME, ":", nrow(species_pam), "\n\n")

# ==============================================================================
# LOAD SHAPEFILES
# ==============================================================================

cat("=== LOADING KDE SHAPEFILES ===\n\n")

# Load quantile-based shapefile
quantile_shp_pattern <- paste0("KDE_pam_", SPECIES_CODE, "\\.shp$")
quantile_shp_files <- list.files(QUANTILE_SHAPEFILE_DIR, 
                                 pattern = quantile_shp_pattern, 
                                 full.names = TRUE, ignore.case = TRUE)

# Load detection-based shapefile
detection_shp_pattern <- paste0("KDE_pam_detbased_", SPECIES_CODE, "\\.shp$")
detection_shp_files <- list.files(DETECTION_SHAPEFILE_DIR, 
                                  pattern = detection_shp_pattern, 
                                  full.names = TRUE, ignore.case = TRUE)

if (length(quantile_shp_files) == 0) {
  cat("WARNING: No quantile-based shapefile found\n")
  quantile_kde <- NULL
} else {
  cat("Loaded quantile-based shapefile\n")
  quantile_kde <- st_read(quantile_shp_files[1], quiet = TRUE) %>%
    st_transform(UTM20)
  
  # Find quantile column
  q_col <- intersect(names(quantile_kde), c("Quantile", "Quantil", "quantile"))[1]
  if (!is.na(q_col)) {
    quantile_kde <- quantile_kde %>% rename(Level = !!q_col)
  }
}

if (length(detection_shp_files) == 0) {
  cat("WARNING: No detection-based shapefile found\n")
  detection_kde <- NULL
} else {
  cat("Loaded detection-based shapefile\n")
  detection_kde <- st_read(detection_shp_files[1], quiet = TRUE) %>%
    st_transform(UTM20)
  
  # Find threshold column
  t_col <- intersect(names(detection_kde), c("Threshold", "threshold"))[1]
  if (!is.na(t_col)) {
    detection_kde <- detection_kde %>% rename(Level = !!t_col)
  }
}

cat("\n")

# ==============================================================================
# CHECK STATION CONTAINMENT - QUANTILE-BASED
# ==============================================================================

if (!is.null(quantile_kde)) {
  cat("=== QUANTILE-BASED CONTAINMENT ===\n\n")
  
  station_containment_q <- species_pam %>%
    rowwise() %>%
    mutate(
      intersects_kde = list(st_intersects(geometry, quantile_kde, sparse = FALSE)[1,]),
      containing_levels = list(quantile_kde$Level[unlist(intersects_kde)])
    ) %>%
    ungroup()
  
  summary_q <- station_containment_q %>%
    st_drop_geometry() %>%
    mutate(
      n_polygons = lengths(containing_levels),
      level_list = sapply(containing_levels, function(x) {
        if(length(x) == 0) return("NONE")
        paste(sort(unique(as.character(x))), collapse = ", ")
      }),
      in_any_polygon = n_polygons > 0
    ) %>%
    select(site, proportion_det, in_any_polygon, level_list) %>%
    arrange(proportion_det)
  
  cat("Stations in quantile-based contours:", sum(summary_q$in_any_polygon), "/", nrow(summary_q), "\n")
  
  not_in_quantile <- summary_q %>%
    filter(proportion_det > 0, !in_any_polygon)
  
  if (nrow(not_in_quantile) > 0) {
    cat("Stations WITH detections NOT in quantile contours:\n")
    print(not_in_quantile)
  }
  cat("\n")
}

# ==============================================================================
# CHECK STATION CONTAINMENT - DETECTION-BASED
# ==============================================================================

if (!is.null(detection_kde)) {
  cat("=== DETECTION-BASED CONTAINMENT ===\n\n")
  
  station_containment_d <- species_pam %>%
    rowwise() %>%
    mutate(
      intersects_kde = list(st_intersects(geometry, detection_kde, sparse = FALSE)[1,]),
      containing_levels = list(detection_kde$Level[unlist(intersects_kde)])
    ) %>%
    ungroup()
  
  summary_d <- station_containment_d %>%
    st_drop_geometry() %>%
    mutate(
      n_polygons = lengths(containing_levels),
      level_list = sapply(containing_levels, function(x) {
        if(length(x) == 0) return("NONE")
        paste(sort(unique(as.character(x))), collapse = ", ")
      }),
      in_any_polygon = n_polygons > 0
    ) %>%
    select(site, proportion_det, in_any_polygon, level_list) %>%
    arrange(proportion_det)
  
  cat("Stations in detection-based contours:", sum(summary_d$in_any_polygon), "/", nrow(summary_d), "\n")
  
  not_in_detection <- summary_d %>%
    filter(proportion_det > 0, !in_any_polygon)
  
  if (nrow(not_in_detection) > 0) {
    cat("Stations WITH detections NOT in detection contours:\n")
    print(not_in_detection)
  } else {
    cat("✓ ALL stations with detections are in detection-based contours!\n")
  }
  cat("\n")
}

# ==============================================================================
# DETAILED CHECK FOR SPECIFIC STATIONS
# ==============================================================================

cat("=== DETAILED CHECK FOR PROBLEM STATIONS ===\n\n")

for (station_name in STATIONS_TO_CHECK) {
  station_data <- species_pam %>%
    filter(site == station_name)
  
  if (nrow(station_data) == 0) {
    cat("Station", station_name, "not found!\n\n")
    next
  }
  
  cat("Station:", station_name, "\n")
  cat("  Detection proportion:", station_data$proportion_det, "\n")
  cat("  Detection days:", station_data$detection_days, "/", station_data$total_days, "\n")
  
  # Check quantile-based
  if (!is.null(quantile_kde)) {
    cat("  Quantile-based contours:\n")
    for (level in unique(quantile_kde$Level)) {
      level_poly <- quantile_kde %>% filter(Level == level)
      intersects <- st_intersects(station_data, level_poly, sparse = FALSE)[1,1]
      cat("    ", level, ":", ifelse(intersects, "YES", "NO"), "\n")
    }
  }
  
  # Check detection-based
  if (!is.null(detection_kde)) {
    cat("  Detection-based contours:\n")
    for (level in unique(detection_kde$Level)) {
      level_poly <- detection_kde %>% filter(Level == level)
      intersects <- st_intersects(station_data, level_poly, sparse = FALSE)[1,1]
      cat("    ", level, ":", ifelse(intersects, "YES", "NO"), "\n")
    }
  }
  
  cat("\n")
}

# ==============================================================================
# VISUALIZATION 1: SIDE-BY-SIDE COMPARISON
# ==============================================================================

cat("=== CREATING COMPARISON VISUALIZATIONS ===\n\n")

# Base map theme
map_theme <- theme_minimal() +
  theme(
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# Quantile-based map
if (!is.null(quantile_kde)) {
  p_quantile <- ggplot() +
    geom_sf(data = quantile_kde, aes(fill = Level), alpha = 0.5) +
    geom_sf(data = species_pam %>% filter(proportion_det == 0),
            color = "gray20", fill = "gray80", size = 2, shape = 21) +
    geom_sf(data = species_pam %>% filter(proportion_det > 0),
            aes(size = proportion_det, color = proportion_det), alpha = 0.8) +
    geom_sf_label(data = species_pam %>% filter(site %in% STATIONS_TO_CHECK),
                  aes(label = site), nudge_y = 10000, size = 2.5) +
    scale_fill_viridis_d(name = "Quantile") +
    scale_color_gradient(low = "orange", high = "darkorange", 
                         name = "Detection\nProportion") +
    scale_size_continuous(range = c(1, 4), name = "Detection\nProportion") +
    labs(title = "Quantile-based KDE",
         subtitle = paste0(sum(summary_q$in_any_polygon), "/", 
                           nrow(summary_q), " stations in contours")) +
    map_theme
} else {
  p_quantile <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No quantile-based\nshapefile", size = 6) +
    theme_void()
}

# Detection-based map
if (!is.null(detection_kde)) {
  p_detection <- ggplot() +
    geom_sf(data = detection_kde, aes(fill = Thrshld), alpha = 0.5) +
    geom_sf(data = species_pam %>% filter(proportion_det == 0),
            color = "gray20", fill = "gray80", size = 2, shape = 21) +
    geom_sf(data = species_pam %>% filter(proportion_det > 0),
            aes(size = proportion_det, color = proportion_det), alpha = 0.8) +
    geom_sf_label(data = species_pam %>% filter(site %in% STATIONS_TO_CHECK),
                  aes(label = site), nudge_y = 10000, size = 2.5) +
    scale_fill_viridis_d(name = "Detection\nThreshold") +
    scale_color_gradient(low = "orange", high = "darkorange", 
                         name = "Detection\nProportion") +
    scale_size_continuous(range = c(1, 4), name = "Detection\nProportion") +
    labs(title = "Detection-based KDE",
         subtitle = paste0(sum(summary_d$in_any_polygon), "/", 
                           nrow(summary_d), " stations in contours")) +
    map_theme
} else {
  p_detection <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No detection-based\nshapefile", size = 6) +
    theme_void()
}

# Combine
p_combined <- p_quantile + p_detection +
  plot_annotation(
    title = paste(SPECIES_NAME, "- KDE Contour Comparison"),
    subtitle = "Gray points = no detections | Orange points = detection stations"
  )

ggsave(file.path(OUTPUT_DIR, paste0(SPECIES_CODE, "_comparison_side_by_side.png")),
       p_combined, width = 16, height = 8, dpi = 300, bg = "white")

cat("✓ Saved side-by-side comparison\n")

# ==============================================================================
# VISUALIZATION 2: STATION STATUS COMPARISON
# ==============================================================================

if (!is.null(quantile_kde) && !is.null(detection_kde)) {
  
  # Combine containment info
  comparison_df <- data.frame(
    site = species_pam$site,
    proportion_det = species_pam$proportion_det,
    in_quantile = summary_q$in_any_polygon,
    in_detection = summary_d$in_any_polygon
  ) %>%
    mutate(
      status = case_when(
        proportion_det == 0 ~ "No detections",
        in_quantile & in_detection ~ "In both",
        in_detection & !in_quantile ~ "Detection-based ONLY",
        in_quantile & !in_detection ~ "Quantile-based ONLY",
        TRUE ~ "In neither"
      )
    )
  
  # Add geometry
  comparison_sf <- species_pam %>%
    left_join(comparison_df, by = c("site", "proportion_det"))
  
  # Create status map
  p_status <- ggplot() +
    geom_sf(data = quantile_kde, fill = "lightblue", alpha = 0.2) +
    geom_sf(data = detection_kde, fill = "lightyellow", alpha = 0.2) +
    geom_sf(data = comparison_sf,
            aes(color = status, size = proportion_det),
            alpha = 0.8) +
    geom_sf_label(data = comparison_sf %>% filter(site %in% STATIONS_TO_CHECK),
                  aes(label = site), nudge_y = 10000, size = 3) +
    scale_color_manual(
      values = c("No detections" = "gray50",
                 "In both" = "darkgreen",
                 "Detection-based ONLY" = "red",
                 "Quantile-based ONLY" = "blue",
                 "In neither" = "black"),
      name = "Station Status"
    ) +
    scale_size_continuous(range = c(2, 5), name = "Detection\nProportion") +
    labs(title = paste(SPECIES_NAME, "- Station Containment Status"),
         subtitle = "Blue background = quantile contours | Yellow background = detection contours\nRed stations = ONLY in detection-based (FIXED!)") +
    map_theme
  
  ggsave(file.path(OUTPUT_DIR, paste0(SPECIES_CODE, "_station_status.png")),
         p_status, width = 12, height = 10, dpi = 300, bg = "white")
  
  cat("✓ Saved station status map\n")
  
  # Summary table
  cat("\n=== STATION STATUS SUMMARY ===\n\n")
  status_summary <- comparison_df %>%
    count(status) %>%
    arrange(desc(n))
  print(status_summary)
  
  # Save CSV
  write_csv(comparison_df, 
            file.path(OUTPUT_DIR, paste0(SPECIES_CODE, "_station_comparison.csv")))
}

# ==============================================================================
# VISUALIZATION 3: DETECTION PROPORTION VS CONTOUR LEVEL
# ==============================================================================

if (!is.null(detection_kde)) {
  
  # Get the lowest contour level each station falls into
  station_levels <- summary_d %>%
    filter(in_any_polygon) %>%
    mutate(
      lowest_level = sapply(strsplit(level_list, ", "), function(x) x[1])
    )
  
  # Join with detection proportions
  level_analysis <- species_pam %>%
    st_drop_geometry() %>%
    left_join(station_levels %>% select(site, lowest_level), by = "site") %>%
    mutate(
      lowest_level = ifelse(is.na(lowest_level), "Not in contour", lowest_level),
      lowest_level = factor(lowest_level, 
                            levels = c("Not in contour", ">1%", ">2%", ">5%", 
                                       ">10%", ">20%", ">30%"))
    )
  
  p_levels <- ggplot(level_analysis, aes(x = lowest_level, y = proportion_det)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(aes(color = site %in% STATIONS_TO_CHECK), 
                width = 0.2, height = 0, alpha = 0.6, size = 3) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red"),
                       name = "Problem\nStations") +
    labs(title = paste(SPECIES_NAME, "- Detection Proportion by Contour Level"),
         subtitle = "Shows which detection rates fall into which contours",
         x = "Lowest Contour Level",
         y = "Detection Proportion") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(OUTPUT_DIR, paste0(SPECIES_CODE, "_detection_by_level.png")),
         p_levels, width = 10, height = 6, dpi = 300, bg = "white")
  
  cat("✓ Saved detection vs level plot\n")
}

# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== DIAGNOSTIC SUMMARY ===\n\n")
cat("Species:", SPECIES_NAME, "(", SPECIES_CODE, ")\n")
cat("Total stations:", nrow(species_pam), "\n")
cat("Stations with detections:", sum(species_pam$proportion_det > 0), "\n\n")

if (!is.null(quantile_kde)) {
  cat("QUANTILE-BASED:\n")
  cat("  Stations in contours:", sum(summary_q$in_any_polygon), "\n")
  cat("  Stations with detections NOT in contours:", 
      nrow(summary_q %>% filter(proportion_det > 0, !in_any_polygon)), "\n\n")
}

if (!is.null(detection_kde)) {
  cat("DETECTION-BASED:\n")
  cat("  Stations in contours:", sum(summary_d$in_any_polygon), "\n")
  cat("  Stations with detections NOT in contours:", 
      nrow(summary_d %>% filter(proportion_det > 0, !in_any_polygon)), "\n\n")
}

cat("Files saved to:", OUTPUT_DIR, "\n")
cat("  -", paste0(SPECIES_CODE, "_comparison_side_by_side.png\n"))
cat("  -", paste0(SPECIES_CODE, "_station_status.png\n"))
cat("  -", paste0(SPECIES_CODE, "_detection_by_level.png\n"))
cat("  -", paste0(SPECIES_CODE, "_station_comparison.csv\n"))

cat("\n=== DONE ===\n")

