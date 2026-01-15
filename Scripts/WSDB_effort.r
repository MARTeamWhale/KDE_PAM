# Opportunistic cetacean DB mapping (OPTIMIZED VERSION)
# Outputs (per target): 3 PNG maps
#   1) All cetacean records per grid cell (proxy for reporting coverage)
#   2) Target records per grid cell
#   3) Target share = target / all (masked where all < min_records_all)
#
# OPTIMIZATION NOTES:
# - Selective column reading
# - Combined filter operations
# - Pre-compiled regex patterns
# - Grid calculation moved outside loop
# - Optional data.table for large datasets
# - Optional parallel processing for multiple targets

# ---------------------------
# PACKAGE SETUP
# ---------------------------
pkgs <- c("dplyr", "readr", "lubridate", "ggplot2", "stringr", "tidyr", "sf", "ggspatial", "ggpattern")

# Optional but recommended for large datasets
use_datatable <- TRUE
use_parallel <- TRUE  # Only beneficial with 3+ targets

if (use_datatable) pkgs <- c(pkgs, "data.table")
if (use_parallel) pkgs <- c(pkgs, "future", "furrr")

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  message("Installing missing packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
}
invisible(lapply(pkgs, library, character.only = TRUE))

# ---------------------------
# USER EDITS
# ---------------------------
infile <- "input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv"   # path to your CSV
out_dir <- "output/Effort/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Grid size (km): shelf-only overview usually reads well at 25 km
grid_km <- 25

# Mask share map where reporting is too thin
min_records_all <- 5

# Date filter
year_min <- 2000
year_max <- 2025

# Optional spatial filter in lon/lat. Set to NULL to skip.
# bbox_lonlat <- c(xmin = -68, xmax = -54, ymin = 40, ymax = 48)
bbox_lonlat <- NULL

# OWA wind area polygons (sf object or shapefile path)
# Set to path of your OSW wind areas shapefile/geopackage, or NULL to skip
owa_polygons_path <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"
owa_color <- "#FF6B35"  # Orange color for OWA boundaries
owa_alpha <- 0.3
owa_linewidth <- 1

# Study area boundary (optional, set to NULL to skip)
study_area_path <- "shapefiles/studyArea/ESS_study_area_simple.shp"

# Land polygons (optional, for context)
land_path <- "shapefiles/Canada/NE_10mLand.shp"

# Log transform the "all records" map to reduce Gully bias?
log_transform_all <- TRUE

# Weighting:
# Default uses each row as 1 record (recommended for opportunistic reporting).
# If you are sure records_in_group is the right unit, set weight_col <- "records_in_group"
weight_col <- NULL

# Targets for SHELF-ONLY OWA workshop (NARW + likely overlap groups)
# Matching runs against "scientific_name | common_name"
targets <- list(
  NARW = "Eubalaena\\s+glacialis|RIGHT\\s+WHALE|N\\s*ATLANTIC\\s+RIGHT",
  Blue_whale = "Balaenoptera\\s+musculus|BLUE\\s+WHALE",
  Shelf_baleen = "Megaptera\\s+novaeangliae|HUMPBACK|Balaenoptera\\s+acutorostrata|MINKE|Balaenoptera\\s+physalus|FIN\\s+WHALE|Balaenoptera\\s+borealis|SEI\\s+WHALE|Balaenoptera\\s+musculus|BLUE\\s+WHALE|Eubalaena\\s+glacialis|RIGHT\\s+WHALE|RORQUAL|BALEEN",
  Harbour_porpoise = "Phocoena\\s+phocoena|HARBOU?R\\s+PORPOISE",
  Small_odontocetes = "DOLPHIN|Delphin|Stenella|Tursiops|Lagenorhynchus|Grampus|WHITE[-\\s]?SIDED|WHITE[-\\s]?BEAKED|COMMON\\s+DOLPHIN|STRIPED\\s+DOLPHIN|SPOTTED\\s+DOLPHIN|BOTTLENOSE|RISSO|Globicephala|PILOT\\s+WHALE",
  Deep_divers = "Physeter|SPERM\\s+WHALE|Mesoplodon|Hyperoodon|Ziphius|BEAKED\\s+WHALE|SOWERBY|CUVIER|NORTHERN\\s+BOTTLENOSE"
)

# Plain language titles for maps (order must match targets list)
target_titles <- list(
  NARW = "North Atlantic Right Whale",
  Blue_whale = "Blue Whale",
  Shelf_baleen = "Baleen Whales",
  Harbour_porpoise = "Harbour Porpoise",
  Small_odontocetes = "Small Odontocetes",
  Deep_divers = "Deep Divers"
)

# Color palettes for each target group (to signal non-comparability of raw counts)
target_palettes <- list(
  NARW = "inferno",              # Red/orange - signals rare/important species
  Blue_whale = "cividis",        # Blue/yellow - for blue whale
  Shelf_baleen = "magma",        # Dark purple/orange 
  Harbour_porpoise = "turbo",    # Rainbow
  Small_odontocetes = "plasma",  # Purple/pink/orange - abundant species
  Deep_divers = "rocket"         # Purple/pink - distinct from all other palettes
)

# Map output settings
map_width <- 10
map_height <- 8
map_dpi <- 300

# ---------------------------
# READ + CLEAN (OPTIMIZED)
# ---------------------------
message("Reading data from: ", infile)

# Define required columns
required_cols <- c("lon", "lat", "x_m", "y_m", "date_utc", 
                   "scientific_name", "common_name")
optional_cols <- c("is_duplicate", if (!is.null(weight_col)) weight_col else NULL)

# Read only necessary columns
all_cols <- c(required_cols, optional_cols)
df <- readr::read_csv(
  infile, 
  col_select = all_of(all_cols),
  show_col_types = FALSE
)

message("Initial rows: ", nrow(df))

# Combined data cleaning and filtering
df <- df %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    x_m = as.numeric(x_m),
    y_m = as.numeric(y_m),
    date_utc = lubridate::as_date(date_utc),
    year = lubridate::year(date_utc),
    common_name = as.character(common_name),
    scientific_name = as.character(scientific_name),
    # Handle is_duplicate in one step
    is_duplicate = if("is_duplicate" %in% names(.)) {
      as.logical(is_duplicate)
    } else {
      FALSE
    }
  ) %>%
  # Combined filter operation (removed scientific_name filter to keep family-level records)
  filter(
    !is.na(x_m), 
    !is.na(y_m), 
    !is.na(date_utc), 
    year >= year_min, 
    year <= year_max
  )

# Optional bbox filter
if (!is.null(bbox_lonlat)) {
  df <- df %>%
    filter(
      !is.na(lon), !is.na(lat),
      between(lon, bbox_lonlat["xmin"], bbox_lonlat["xmax"]),
      between(lat, bbox_lonlat["ymin"], bbox_lonlat["ymax"])
    )
}

message("Rows after filtering: ", nrow(df))

# Weighting
df <- df %>%
  mutate(
    w = if (!is.null(weight_col) && weight_col %in% names(.)) {
      as.numeric(.data[[weight_col]])
    } else {
      1
    },
    w = ifelse(is.na(w) | w <= 0, 1, w)
  )

# Combined name field for matching
df <- df %>%
  mutate(name_blob = paste0(scientific_name, " | ", common_name))

# ---------------------------
# GRID CALCULATION (ONCE)
# ---------------------------
message("Calculating grid cells...")
cell_m <- grid_km * 1000

df <- df %>%
  mutate(
    gx = floor(x_m / cell_m),
    gy = floor(y_m / cell_m),
    cell_id = paste0(gx, "_", gy),
    x0 = gx * cell_m,
    y0 = gy * cell_m
  )

# ---------------------------
# PRE-COMPILE REGEX PATTERNS
# ---------------------------
message("Compiling regex patterns...")
target_patterns <- lapply(targets, function(x) {
  stringr::regex(x, ignore_case = TRUE)
})

# ---------------------------
# LOAD SPATIAL OVERLAYS
# ---------------------------
# Load OWA polygons and other spatial layers
osw_wind <- NULL
study_area <- NULL
land <- NULL

if (!is.null(owa_polygons_path) && file.exists(owa_polygons_path)) {
  message("Loading OWA wind area polygons...")
  osw_wind <- sf::st_read(owa_polygons_path, quiet = TRUE)
  # Transform to match data CRS if needed (assume UTM zone 20N based on x_m, y_m)
  osw_wind <- sf::st_transform(osw_wind, crs = 32620)
}

if (!is.null(study_area_path) && file.exists(study_area_path)) {
  message("Loading study area boundary...")
  study_area <- sf::st_read(study_area_path, quiet = TRUE)
  study_area <- sf::st_transform(study_area, crs = 32620)
}

if (!is.null(land_path) && file.exists(land_path)) {
  message("Loading land polygons...")
  land <- sf::st_read(land_path, quiet = TRUE)
  land <- sf::st_transform(land, crs = 32620)
}

# ---------------------------
# CONVERT TO DATA.TABLE (OPTIONAL)
# ---------------------------
if (use_datatable) {
  message("Converting to data.table for faster aggregation...")
  dt <- data.table::as.data.table(df)
  data.table::setkey(dt, cell_id)
} else {
  dt <- df
}

# ---------------------------
# SET MAP EXTENT FROM DATA GRID
# ---------------------------
# Use data grid extent instead of study area polygon to avoid cropping data
message("Setting map extent from data grid...")
xlims <- range(dt$x0, na.rm = TRUE) + c(-cell_m/2, cell_m/2)
ylims <- range(dt$y0, na.rm = TRUE) + c(-cell_m/2, cell_m/2)
message(sprintf("Map extent: X [%.0f, %.0f], Y [%.0f, %.0f]", xlims[1], xlims[2], ylims[1], ylims[2]))

# Crop land to map extent if land exists
if (!is.null(land)) {
  message("Cropping land to map extent...")
  bbox_crop <- sf::st_bbox(c(xmin = xlims[1], xmax = xlims[2], 
                             ymin = ylims[1], ymax = ylims[2]), 
                           crs = sf::st_crs(32620))
  land <- sf::st_crop(land, bbox_crop)
}

# ---------------------------
# PLOTTING SETTINGS
# ---------------------------
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

# Cache theme to avoid recreating
cached_theme <- theme_map()

# ---------------------------
# MAPPING FUNCTION
# ---------------------------
create_maps_for_target <- function(nm, target_regex, dt, cell_all_with_conf, use_dt = use_datatable) {
  
  message("Processing target: ", nm)
  
  if (use_dt) {
    # DATA.TABLE APPROACH (FASTER)
    cell_sum <- dt[, .(
      n_all = sum(w, na.rm = TRUE),
      n_target = sum(w * stringr::str_detect(name_blob, target_regex), na.rm = TRUE)
    ), by = .(cell_id, gx, gy, x0, y0)]
    
    cell_sum[, share_target := fifelse(
      n_all >= min_records_all & n_all > 0,
      n_target / n_all,
      NA_real_
    )]
    
    # Calculate target-specific years
    target_years <- dt[
      stringr::str_detect(name_blob, target_regex),
      .(n_target_years = data.table::uniqueN(year, na.rm = TRUE)),
      by = .(cell_id)
    ]
    
    # Convert back to tibble for ggplot
    cell_sum <- as_tibble(cell_sum)
    target_years <- as_tibble(target_years)
    
  } else {
    # DPLYR APPROACH
    cell_sum <- dt %>%
      mutate(is_target = stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id, gx, gy, x0, y0) %>%
      summarise(
        n_all = sum(w, na.rm = TRUE),
        n_target = sum(w[is_target], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(
        share_target = ifelse(
          n_all >= min_records_all & n_all > 0, 
          n_target / n_all, 
          NA_real_
        )
      )
    
    # Calculate target-specific years
    target_years <- dt %>%
      filter(stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(n_target_years = n_distinct(year, na.rm = TRUE), .groups = "drop")
  }
  
  # Merge target years
  cell_sum <- cell_sum %>%
    left_join(target_years, by = "cell_id") %>%
    mutate(n_target_years = tidyr::replace_na(n_target_years, 0))
  
  # Create target-level evidence strength
  cell_sum <- cell_sum %>%
    mutate(
      evidence = case_when(
        n_target == 0 ~ "None",
        n_target >= 5 & n_target_years >= 3 ~ "High",
        n_target >= 2 & n_target_years >= 2 ~ "Medium",
        TRUE ~ "Low"
      ),
      evidence = factor(evidence, levels = c("High", "Medium", "Low", "None"))
    )
  
  # Merge confidence data
  cell_sum <- cell_sum %>%
    left_join(
      cell_all_with_conf %>% select(cell_id, confidence),
      by = "cell_id"
    )
  
  # Combine dataset coverage + target evidence into target-specific confidence
  cell_sum <- cell_sum %>%
    mutate(
      target_conf = case_when(
        confidence %in% c("No Data", "Low") ~ "Low (coverage)",
        evidence == "None" ~ "No evidence",
        confidence == "Low-Medium" & evidence == "Low" ~ "Low",
        confidence %in% c("Medium", "Medium-High") & evidence == "Medium" ~ "Moderate",
        confidence %in% c("Medium-High", "High") & evidence == "High" ~ "High",
        TRUE ~ "Moderate"
      ),
      target_conf = factor(target_conf, levels = c("High", "Moderate", "Low", "Low (coverage)", "No evidence"))
    )
  
  # Generate safe filename
  nm_safe <- stringr::str_replace_all(nm, "[^A-Za-z0-9]+", "_")
  
  # Map 1: Target records
  p1 <- ggplot(cell_sum) +
    geom_tile(aes(x = x0, y = y0, fill = n_target)) +
    scale_fill_viridis_c(option = target_palettes[[nm]], na.value = "grey90", name = "Count",
                         begin = 0.10, end = 0.95) +  # Clip dark and light extremes
    labs(
      title = paste0("Opportunistic sightings: ", target_titles[[nm]]),
      caption = "Raw counts reflect both species distribution AND observer effort. Blank cells indicate no reports (not necessarily species absence).\nNote: Color scales vary by species group - use normalized maps for cross-group comparisons."
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
    coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
    cached_theme
  
  # Add spatial overlays in order: land, then OWA, then study area
  if (!is.null(land)) {
    p1 <- p1 + geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
  }
  if (!is.null(osw_wind)) {
    p1 <- p1 + geom_sf(data = osw_wind, fill = NA, color = owa_color, 
                       alpha = owa_alpha, linewidth = owa_linewidth, inherit.aes = FALSE)
  }
  if (!is.null(study_area)) {
    p1 <- p1 + geom_sf(data = study_area, fill = NA, color = "black", 
                       linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
  }
  
  # Add scale bar
  p1 <- p1 + ggspatial::annotation_scale(location = "br", width_hint = 0.25)
  
  ggsave(
    file.path(out_dir, paste0(nm_safe, "_1_target_records.png")),
    p1,
    width = map_width,
    height = map_height,
    dpi = map_dpi
  )
  
  # Map 2: Normalized share WITH DIAGONAL PATTERN FOR LOW COVERAGE------
  # Create a dataset of low-coverage cells and mark them
  cell_sum_with_coverage <- cell_sum %>%
    mutate(low_coverage = target_conf %in% c("Low (coverage)", "No evidence", "Low"))
  
  # Get low coverage cells for outline
  low_conf_cells <- cell_sum %>%
    filter(target_conf %in% c("Low (coverage)", "No evidence", "Low"))
  
  p2 <- ggplot(cell_sum_with_coverage) +
    ggpattern::geom_tile_pattern(
      aes(x = x0, y = y0, fill = share_target, 
          pattern = ifelse(low_coverage, "stripe", "none")),
      pattern_fill = "white",
      pattern_color = "white", 
      pattern_density = 0.03,  # Fewer lines, more spaced out
      pattern_spacing = 0.015,  # Wider spacing between lines
      pattern_angle = 45,
      pattern_alpha = 1.0,
      pattern_linewidth = 0.5  # Thin line width
    ) +
    scale_fill_viridis_c(option = "viridis", na.value = "grey90", limits = c(0, 1), 
                         name = "Proportion of\nall sightings",
                         direction = 1,  
                         begin = 0.1, end = 0.95) +  # Clip to avoid extreme yellow
    scale_pattern_manual(values = c("none" = "none", "stripe" = "stripe"), guide = "none") +
    # Add white outline for low-coverage cells
    {if(nrow(low_conf_cells) > 0) {
      geom_tile(data = low_conf_cells, aes(x = x0, y = y0), 
                fill = NA, color = "white", linewidth = 0.6)
    }} +
    labs(
      title = paste0("Effort-normalized sightings: ", target_titles[[nm]]),
      caption = paste0("Normalization: divides target species records by ALL cetacean records in each cell, controlling for spatial variation in observer effort.\n",
                       "Higher values indicate areas where this species/group comprises a larger share of cetacean encounters, independent of overall reporting intensity.\n",
                       "Masked where total records < ", min_records_all, " per cell. White diagonal hatching indicates low species-specific data coverage (weak evidence or poor dataset coverage).")
    ) +
    annotate("text", x = -Inf, y = Inf, 
             label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
             hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
    coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
    cached_theme
  
  # Add spatial overlays in order: land, then OWA, then study area
  if (!is.null(land)) {
    p2 <- p2 + geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
  }
  if (!is.null(osw_wind)) {
    p2 <- p2 + geom_sf(data = osw_wind, fill = NA, color = owa_color, 
                       alpha = owa_alpha, linewidth = owa_linewidth, inherit.aes = FALSE)
  }
  if (!is.null(study_area)) {
    p2 <- p2 + geom_sf(data = study_area, fill = NA, color = "black", 
                       linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
  }
  
  # Add scale bar
  p2 <- p2 + ggspatial::annotation_scale(location = "br", width_hint = 0.25)
  
  ggsave(
    file.path(out_dir, paste0(nm_safe, "_2_share.png")),
    p2,
    width = map_width,
    height = map_height,
    dpi = map_dpi
  )
  
  message("  ✓ Saved 2 maps for ", nm)
  
  return(invisible(NULL))
}

# ---------------------------
# CREATE SINGLE "ALL RECORDS" MAP
# ---------------------------
message("\n=== Creating 'all records' coverage map ===")

if (use_datatable) {
  cell_all <- dt[, .(
    n_all = sum(w, na.rm = TRUE)
  ), by = .(cell_id, gx, gy, x0, y0)]
  cell_all <- as_tibble(cell_all)
} else {
  cell_all <- dt %>%
    group_by(cell_id, gx, gy, x0, y0) %>%
    summarise(n_all = sum(w, na.rm = TRUE), .groups = "drop")
}

# ---------------------------
# CALCULATE CONFIDENCE INDEX
# ---------------------------
message("Calculating confidence index based on coverage + temporal repeatability...")

# Calculate number of distinct years per cell
if (use_datatable) {
  cell_years <- dt[, .(
    n_years = data.table::uniqueN(year, na.rm = TRUE)
  ), by = .(cell_id)]
  cell_years <- as_tibble(cell_years)
} else {
  cell_years <- dt %>%
    group_by(cell_id) %>%
    summarise(n_years = n_distinct(year, na.rm = TRUE), .groups = "drop")
}

# Merge years into cell_all
cell_all <- cell_all %>%
  left_join(cell_years, by = "cell_id") %>%
  mutate(n_years = ifelse(is.na(n_years), 0, n_years))

# Define confidence bins based on quantiles
# Use cells with data only (n_all > 0) for quantile calculation
cells_with_data <- cell_all %>% filter(n_all > 0)

# Calculate quantile thresholds for n_all (33rd, 67th percentiles)
n_all_q33 <- quantile(cells_with_data$n_all, 0.33, na.rm = TRUE)
n_all_q67 <- quantile(cells_with_data$n_all, 0.67, na.rm = TRUE)

# Calculate quantile thresholds for n_years (33rd, 67th percentiles)
n_years_q33 <- quantile(cells_with_data$n_years, 0.33, na.rm = TRUE)
n_years_q67 <- quantile(cells_with_data$n_years, 0.67, na.rm = TRUE)

message(sprintf("  Coverage thresholds: Low < %.1f | Medium %.1f-%.1f | High > %.1f", 
                n_all_q33, n_all_q33, n_all_q67, n_all_q67))
message(sprintf("  Temporal thresholds: Low < %.1f yrs | Medium %.1f-%.1f yrs | High > %.1f yrs", 
                n_years_q33, n_years_q33, n_years_q67, n_years_q67))

# Assign confidence categories
# Logic: Both coverage AND repeatability should be considered
# High = both high; Medium = at least one high OR both medium; Low = both low or one very low
cell_all <- cell_all %>%
  mutate(
    coverage_cat = case_when(
      n_all == 0 ~ "None",
      n_all < n_all_q33 ~ "Low",
      n_all < n_all_q67 ~ "Medium",
      TRUE ~ "High"
    ),
    temporal_cat = case_when(
      n_years == 0 ~ "None",
      n_years < n_years_q33 ~ "Low",
      n_years < n_years_q67 ~ "Medium",
      TRUE ~ "High"
    ),
    confidence = case_when(
      n_all == 0 | n_years == 0 ~ "No Data",
      coverage_cat == "High" & temporal_cat == "High" ~ "High",
      coverage_cat == "High" | temporal_cat == "High" ~ "Medium-High",
      coverage_cat == "Medium" & temporal_cat == "Medium" ~ "Medium",
      coverage_cat == "Medium" | temporal_cat == "Medium" ~ "Low-Medium",
      TRUE ~ "Low"
    ),
    # Order confidence levels (reversed so High is at top of legend)
    confidence = factor(confidence, 
                        levels = c("High", "Medium-High", "Medium", 
                                   "Low-Medium", "Low"))
  )

# Apply log transform if requested
if (log_transform_all) {
  cell_all <- cell_all %>%
    mutate(n_all_display = log10(n_all + 1))
  fill_label <- "Records"
  
  # Create breaks at powers of 10 for intuitive reading
  # Log scale: 0, 1, 2, 3, 4 corresponds to: 1, 10, 100, 1000, 10000 records
  log_breaks <- seq(0, ceiling(max(cell_all$n_all_display, na.rm = TRUE)), by = 1)
  actual_values <- 10^log_breaks
  actual_values[1] <- 0  # First break is 0, not 1
  
} else {
  cell_all <- cell_all %>%
    mutate(n_all_display = n_all)
  fill_label <- "Records"
  log_breaks <- NULL
  actual_values <- NULL
}

p_all <- ggplot(cell_all) +
  geom_tile(aes(x = x0, y = y0, fill = n_all_display)) +
  {if (log_transform_all) {
    scale_fill_viridis_c(
      option = "viridis", 
      na.value = "grey90", 
      name = "Count",
      breaks = log_breaks,
      labels = actual_values,
      begin = 0.15, end = 0.95
    )
  } else {
    scale_fill_viridis_c(
      option = "viridis", 
      na.value = "grey90", 
      name = "Count",
      begin = 0.15, end = 0.95
    )
  }} +
  labs(
    title = "Density of all cetacean sightings",
    # caption = "Higher values reflect where reports occur, not necessarily whale abundance."
  ) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid", 
                          if(log_transform_all) " | log scale" else ""),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  cached_theme

# Add spatial overlays in order: land, then OWA, then study area
if (!is.null(land)) {
  p_all <- p_all + geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
}
if (!is.null(osw_wind)) {
  p_all <- p_all + geom_sf(data = osw_wind, fill = NA, color = owa_color, 
                           alpha = owa_alpha, linewidth = owa_linewidth, inherit.aes = FALSE)
}
if (!is.null(study_area)) {
  p_all <- p_all + geom_sf(data = study_area, fill = NA, color = "black", 
                           linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
}

# Add scale bar
p_all <- p_all + ggspatial::annotation_scale(location = "br", width_hint = 0.25)

ggsave(
  file.path(out_dir, "00_all_cetacean_records.png"),
  p_all,
  width = map_width,
  height = map_height,
  dpi = map_dpi
)

message("  ✓ Saved coverage map")

# ---------------------------
# CREATE CONFIDENCE MAP
# ---------------------------
message("Creating confidence map...")

# Create labels with actual threshold values
confidence_subtitle <- sprintf(
  paste0("Based on reporting coverage (Low < %.0f, Med %.0f-%.0f, High > %.0f records) ",
         "and temporal consistency (Low < %.0f, Med %.0f-%.0f, High > %.0f years)"),
  n_all_q33, n_all_q33, n_all_q67, n_all_q67,
  n_years_q33, n_years_q33, n_years_q67, n_years_q67
)

# Define confidence color palette (grey to green - green indicates good confidence)
confidence_colors <- c(
  
  "Low" = "#fee5d9",
  "Low-Medium" = "#fdd49e", 
  "Medium" = "#fdbb84",
  "Medium-High" = "#a1d99b",
  "High" = "#31a354"
)

# Mark low coverage cells for hatching
cell_all_with_pattern <- cell_all %>%
  mutate(low_coverage = confidence %in% c("No Data", "Low", "Low-Medium"))

p_confidence <- ggplot(cell_all_with_pattern) +
  ggpattern::geom_tile_pattern(
    aes(x = x0, y = y0, fill = confidence,
        pattern = ifelse(low_coverage, "stripe", "none")),
    pattern_fill = "white",
    pattern_color = "white", 
    pattern_density = 0.02, 
    pattern_spacing = 0.015,
    pattern_angle = 45,
    pattern_alpha = .5
  ) +
  scale_fill_manual(
    values = confidence_colors,
    name = "Data coverage",
    drop = FALSE
  ) +
  scale_pattern_manual(
    values = c("none" = "none", "stripe" = "stripe"), 
    guide = "none"
  ) +
  guides(fill = guide_legend(override.aes = list(
    pattern = c("none", "none", "none", "stripe", "stripe"),  # Order: High, Med-High, Med, Low-Med, Low
    pattern_fill = "white",
    pattern_color = "white"
  ))) +
  labs(
    title = "Cetacean sighting coverage + temporal consistency",
    caption =  confidence_subtitle) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0(year_min, "–", year_max, " | ", grid_km, " km grid"),
           hjust = -0.1, vjust = 1.5, size = 3.5, color = "grey30") +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  cached_theme

# Add spatial overlays
if (!is.null(land)) {
  p_confidence <- p_confidence + geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
}
if (!is.null(osw_wind)) {
  p_confidence <- p_confidence + geom_sf(data = osw_wind, fill = NA, color = owa_color, 
                                         alpha = owa_alpha, linewidth = owa_linewidth, inherit.aes = FALSE)
}
if (!is.null(study_area)) {
  p_confidence <- p_confidence + geom_sf(data = study_area, fill = NA, color = "black", 
                                         linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE)
}

# Add scale bar
p_confidence <- p_confidence + ggspatial::annotation_scale(location = "br", width_hint = 0.25)

ggsave(
  file.path(out_dir, "00_data_consistency.png"),
  p_confidence,
  width = map_width,
  height = map_height,
  dpi = map_dpi
)

message("  ✓ Saved data consistency map")

# ---------------------------
# PROCESS ALL TARGETS
# ---------------------------
message("\n=== Creating target-specific maps ===")

if (use_parallel && length(targets) >= 3) {
  # PARALLEL PROCESSING
  message("Using parallel processing with ", future::availableCores() - 1, " workers")
  future::plan(future::multisession, workers = min(length(targets), future::availableCores() - 1))
  
  furrr::future_walk(names(targets), function(nm) {
    create_maps_for_target(nm, target_patterns[[nm]], dt, cell_all, use_dt = use_datatable)
  })
  
  future::plan(future::sequential)
  
} else {
  # SEQUENTIAL PROCESSING
  for (nm in names(targets)) {
    create_maps_for_target(nm, target_patterns[[nm]], dt, cell_all, use_dt = use_datatable)
  }
}

message("\n=== COMPLETE ===")
message("All maps saved to: ", out_dir)
message("Total maps created: ", 2 + (length(targets) * 2), " (1 coverage + 1 data consistency + ", length(targets), " targets × 2 maps)")
message("Target groups: ", paste(names(targets), collapse = ", "))


