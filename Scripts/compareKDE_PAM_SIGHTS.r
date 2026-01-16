# ==============================================================================
# Author: Laura Joan Feyrer
# Date: January 13, 2026
# Description: KDE COMPARISON: PAM vs SIGHTINGS - REVISED WITH COMMON NAMES
# 
# Changes from previous version (2026-01-13):
# - MAJOR REVISION: Uses full common names throughout instead of two-letter codes
# - Fixed beaked whale PAM matching (was matching too many files)
# - Ha = Northern Bottlenose Whale, Mb = Sowerby's Beaked Whale
# - Zc = Goose Beaked Whale (NOT Cuvier's), MmMe = True's/Gervais' Beaked Whale
# - More explicit PAM base names to avoid false matches
# ==============================================================================

options(scipen = 999)

pacman::p_load(sf, tidyverse, terra, ggplot2, viridis, ggspatial, patchwork, rlang)

# ==============================================================================
# CONFIGURATION--------
# ==============================================================================

CREATE_SIGHTINGS_ONLY_MAPS <- TRUE

# Projections
UTM20 <- st_crs(32620)

# Input paths
SHAPEFILE_DIR <- "output/shapes/baseline_match_v3/multi"
SIGHTINGS_DATA <- "input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv"
BALEEN_PAM_DATA <- "input/2025/baleen_presence_laura_2025.csv"
BEAKED_PAM_DATA <- "input/2025/beaked_PAM/beaked_pam_results_2026-01-12.csv"
LAND_PATH <- "/Users/chirp/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"
OSW_WIND_PATH <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"

# Output directory
COMPARE_OUTPUT_DIR <- "output/FIGS/compare_sightings_pam/"
if (!dir.exists(COMPARE_OUTPUT_DIR)) dir.create(COMPARE_OUTPUT_DIR, recursive = TRUE)

# Color palettes
QUANTILE_PALETTE <- "mako"
POINT_COLOR <- "gray40"
STATION_COLOR <- "orange"
OSW_COLOR <- "red"
OSW_ALPHA <- 0.3

# ==============================================================================
# SPECIES MAPPING - BALEEN WHALES (PAM + SIGHTINGS)
# ==============================================================================

BALEEN_SPECIES_MAPPING <- list(
  "Sei_Whale" = list(
    pam_code = "Bb",
    common_name = "WHALE-SEI",
    display_name = "Sei Whale",
    sightings_base = "WHALE-SEI",
    pam_base = "pam_baleen_Sei_Whale_Bb"  # Matches "KDE_pam_Bb_..."
  ),
  "Blue_Whale" = list(
    pam_code = "Bm",
    common_name = "WHALE-BLUE",
    display_name = "Blue Whale",
    sightings_base = "WHALE-BLUE",
    pam_base = "pam_baleen_Blue_Whale_Bm"
  ),
  "Fin_Whale" = list(
    pam_code = "Bp",
    common_name = "WHALE-FIN",
    display_name = "Fin Whale",
    sightings_base = "WHALE-FIN",
    pam_base = "pam_baleen_Fin_Whale_Bp"
  ),
  "Humpback_Whale" = list(
    pam_code = "Mn",
    common_name = "WHALE-HUMPBACK",
    display_name = "Humpback Whale",
    sightings_base = "WHALE-HUMPBACK",
    pam_base = "pam_baleen_Humpback_Whale_Mn"
  ),
  "Minke_Whale" = list(
    pam_code = "Ba",
    common_name = "WHALE-MINKE",
    display_name = "Minke Whale",
    sightings_base = "WHALE-MINKE",
    pam_base = "pam_baleen_Minke_Whale_Ba"
  ),
  "North_Atlantic_Right_Whale" = list(
    pam_code = "Eg",
    common_name = "WHALE-NORTH ATLANTIC RIGHT",
    display_name = "North Atlantic Right Whale",
    sightings_base = "WHALE-NORTH_ATLANTIC_RIGHT",
    pam_base = "pam_baleen_North_Atlantic_Right_Whale_Eg"
  )
)

# ==============================================================================
# SPECIES MAPPING - BEAKED WHALES
# ==============================================================================

BEAKED_SPECIES_MAPPING <- list(
  "Northern_Bottlenose_Whale" = list(
    pam_code = "Ha",
    common_name = "WHALE-NORTHERN BOTTLENOSE",
    display_name = "Northern Bottlenose Whale",
    sightings_base = "WHALE-NORTHERN_BOTTLENOSE",
    pam_base = "pam_beaked_Northern_Bottlenose_Whale_Ha",  # Very specific to avoid false matches
    has_sightings = TRUE
  ),
  "Sowerbys_Beaked_Whale" = list(
    pam_code = "Mb",
    common_name = "WHALE-SOWERBY'S BEAKED",
    display_name = "Sowerby's Beaked Whale",
    sightings_base = "WHALE-SOWERBY'S_BEAKED",
    pam_base = "pam_beaked_Sowerbys_Beaked_Whale_Mb",
    has_sightings = TRUE
  ),
  "Goose_Beaked_Whale" = list(
    pam_code = "Zc",
    common_name = "WHALE-GOOSE BEAKED",  # NOT Cuvier's
    display_name = "Goose Beaked Whale",
    sightings_base = "WHALE-GOOSE_BEAKED",
    pam_base = "pam_beaked_Goose_Beaked_Whale_Zc",
    has_sightings = FALSE,
    alternate_sightings_names = c("WHALE- CUVIER'S BEAKED", "WHALE-CUVIER'S BEAKED")  # Search but don't display
  ),
  "Trues_Gervais_Beaked_Whale" = list(
    pam_code = "MmMe",
    common_name = "WHALE-TRUE'S/GERVAIS' BEAKED",
    display_name = "True's/Gervais' Beaked Whale",
    sightings_base = "WHALE-TRUE'S_BEAKED",
    pam_base = "pam_beaked_Trues_Gervais_Beaked_Whale_MmMe",
    has_sightings = FALSE
  )
)

# ==============================================================================
# GROUPED SPECIES MAPPING
# ==============================================================================

GROUPED_BALEEN <- list(
  "All_Baleen_Whales" = list(
    common_name = "ALL_BALEEN_MYSTICETES",
    display_name = "All Baleen Whales",
    sightings_base = "KDE_sightings_mysticetes",
    pam_base = "grouped_baleen_pam_All_Baleen_Whales_(PAM)",
    is_grouped = TRUE
  )
)

GROUPED_BEAKED <- list(
  "All_Beaked_Whales" = list(
    common_name = "ALL_BEAKED_WHALES",
    display_name = "All Beaked Whales",
    sightings_base = NA_character_,
    pam_base = "grouped_beaked_pam_All_Beaked_Whales_(PAM)",
    is_grouped = TRUE,
    pam_only = TRUE
  )
)

GROUPED_ODONTOCETES_SIGHTINGS <- list(
  "All_Odontocetes" = list(
    common_name = "ODONTOCETES",
    display_name = "All Odontocetes (Toothed Whales)",
    sightings_base = "sightings_odontocetes",
    pam_base = NA_character_,
    is_grouped = TRUE,
    pam_only = FALSE
  )
)


# ==============================================================================
# LOAD BASE LAYERS
# ==============================================================================

cat("Loading base layers...\n")

land <- st_read(LAND_PATH, quiet = TRUE) %>%
  filter(CONTINENT == "North America") %>%
  st_transform(UTM20)

study_area <- st_read(STUDY_AREA_PATH, quiet = TRUE) %>%
  st_transform(UTM20)

bbox <- st_bbox(study_area)
xlims <- c(bbox["xmin"], bbox["xmax"])
ylims <- c(bbox["ymin"], bbox["ymax"])

cat("✓ Base layers loaded\n\n")

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

# Load BALEEN PAM data
baleen_pam_raw <- read_csv(BALEEN_PAM_DATA, show_col_types = FALSE) %>%
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

baleen_pam_data <- st_intersection(baleen_pam_raw, study_area)
cat("✓ Loaded", nrow(baleen_pam_data), "baleen PAM station records\n")

# Load BEAKED WHALE PAM data
cat("Loading beaked whale PAM data...\n")
beaked_pam_raw <- read_csv(BEAKED_PAM_DATA, show_col_types = FALSE)

beaked_pam_processed <- beaked_pam_raw %>%
  group_by(deployment, station, latitude, longitude, species) %>%
  summarise(
    detection_days = sum(presence, na.rm = TRUE),
    effort_days = n_distinct(rec_date),
    proportion_det = detection_days / effort_days,
    .groups = "drop"
  ) %>%
  filter(detection_days > 0) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

beaked_pam_data <- st_intersection(beaked_pam_processed, study_area)

cat("✓ Loaded", nrow(beaked_pam_data), "beaked whale PAM station records (with detections)\n")

beaked_species_summary <- beaked_pam_data %>%
  st_drop_geometry() %>%
  group_by(species) %>%
  summarise(
    n_deployments = n(),
    total_det_days = sum(detection_days),
    mean_proportion = mean(proportion_det),
    .groups = "drop"
  ) %>%
  arrange(desc(total_det_days))

cat("\nBeaked whale species in PAM data:\n")
print(beaked_species_summary, n = Inf)
cat("\n")

pam_data <- baleen_pam_data  # Keep for backward compatibility

# ==============================================================================
# HELPER FUNCTIONS-------
# ==============================================================================

extract_kde_params <- function(filename) {
  fname <- basename(filename)
  
  has_proximity <- grepl("proximity", fname, ignore.case = TRUE)
  
  bw <- NA_real_
  m <- regexec("bw[_-]?(\\d+)", fname, ignore.case = TRUE)
  hit <- regmatches(fname, m)[[1]]
  if (length(hit) > 1) bw <- as.numeric(hit[2])
  
  is_diggle <- grepl("diggle", fname, ignore.case = TRUE)
  is_adaptive <- grepl("adaptive", fname, ignore.case = TRUE)
  is_fixed <- grepl("fixed", fname, ignore.case = TRUE)
  
  method <- dplyr::case_when(
    is_diggle ~ "diggle",
    is_adaptive ~ "adaptive",
    is_fixed ~ "fixed",
    TRUE ~ "unknown"
  )
  
  list(
    has_proximity = has_proximity,
    bandwidth = bw,
    method = method
  )
}

find_kde_shapefiles <- function(base_name,
                                 shapefile_dir = SHAPEFILE_DIR,
                                 prefer_proximity = FALSE) {
  
  files <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE, ignore.case = TRUE)
  
  # Escape regex special characters so base_name is treated literally
  base_esc <- gsub("([][{}()+*^$|\\\\.?])", "\\\\\\1", base_name)
  
  # Match base_name as an underscore-delimited token (or start/end of string)
  token_pat <- paste0("(^|_)", base_esc, "(_|$)")
  
  hits <- files[
    grepl("^KDE_", basename(files), ignore.case = TRUE) &
      grepl(token_pat, basename(files), ignore.case = TRUE)
  ]
  
  if (length(hits) == 0) return(NULL)
  
  hits <- hits[!grepl("old|backup|test", basename(hits), ignore.case = TRUE)]
  
  df <- tibble::tibble(
    path = hits,
    file = basename(hits),
    mtime = file.info(hits)$mtime,
    params = purrr::map(hits, extract_kde_params)
  ) |>
    dplyr::mutate(
      has_proximity = purrr::map_lgl(params, ~ .x$has_proximity),
      bandwidth = purrr::map_dbl(params, ~ .x$bandwidth %||% NA_real_),
      method = purrr::map_chr(params, ~ .x$method %||% "unknown")
    )
  
  if (prefer_proximity) {
    df <- df |>
      dplyr::mutate(priority = dplyr::if_else(has_proximity, 0L, 1L)) |>
      dplyr::arrange(priority, dplyr::desc(mtime))
  } else {
    df <- df |>
      dplyr::arrange(dplyr::desc(mtime))
  }
  
  df
}


load_kde_shapefile <- function(pattern_or_file, shapefile_dir = SHAPEFILE_DIR) {
  
  if (file.exists(pattern_or_file) && grepl("\\.shp$", pattern_or_file)) {
    shp_file <- pattern_or_file
  } else {
    shp_files <- list.files(shapefile_dir, 
                            pattern = paste0(pattern_or_file, "\\.shp$"), 
                            full.names = TRUE, 
                            ignore.case = TRUE)
    
    if (length(shp_files) == 0) {
      return(NULL)
    }
    
    if (length(shp_files) > 1) {
      file_times <- file.info(shp_files)$mtime
      shp_file <- shp_files[which.max(file_times)]
      cat("  Note: Multiple matches found, using most recent:", basename(shp_file), "\n")
    } else {
      shp_file <- shp_files[1]
    }
  }
  
  kde_shp <- st_read(shp_file, quiet = TRUE) %>%
    st_transform(UTM20)
  
  quantile_col <- NULL
  possible_names <- c("Quantile", "Quantil", "quant_level", "quantil", "QUANTILE")
  
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
  
  kde_shp <- kde_shp %>%
    rename(Quantile = !!quantile_col)
  
  unique_quantiles <- unique(kde_shp$Quantile)
  unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
  
  if (is.numeric(unique_quantiles)) {
    kde_shp <- kde_shp %>%
      mutate(Quantile = paste0(Quantile * 100, "%"))
  }
  
  quantile_levels <- unique(kde_shp$Quantile)
  quantile_nums <- as.numeric(gsub("%", "", quantile_levels))
  quantile_levels <- quantile_levels[order(quantile_nums)]
  
  kde_shp$Quantile <- factor(kde_shp$Quantile, 
                             levels = quantile_levels, 
                             ordered = TRUE)
  
  return(kde_shp)
}

# ==============================================================================
# FUNCTION: CREATE COMPARISON MAP------
# ==============================================================================
create_comparison_map <- function(species_key, species_info, 
                                  sightings_data, pam_data,
                                  is_grouped = FALSE,
                                  pam_only = FALSE) {
  
  species_name <- species_info$common_name
  cat("Creating", ifelse(pam_only, "PAM-only", "comparison"), "map for:", species_info$display_name, "\n")
  
  # Find KDE files
  has_sightings_kde <- FALSE
  sightings_df <- NULL
  
  if (!pam_only && !is.na(species_info$sightings_base)) {
    cat("  Searching for sightings KDE variants with base name:", species_info$sightings_base, "\n")
    sightings_df <- find_kde_shapefiles(species_info$sightings_base, prefer_proximity = TRUE)
    has_sightings_kde <- !is.null(sightings_df) && nrow(sightings_df) > 0
  }
  
  cat("  Searching for PAM KDE variants with base name:", species_info$pam_base, "\n")
  pam_df <- find_kde_shapefiles(species_info$pam_base, prefer_proximity = FALSE)
  
  has_pam_kde <- !is.null(pam_df) && nrow(pam_df) > 0
  
  if (!has_sightings_kde && !has_pam_kde) {
    cat("  ERROR: No KDE shapefiles found for", species_info$display_name, "\n\n")
    return(NULL)
  }
  
  if (!has_pam_kde && !CREATE_SIGHTINGS_ONLY_MAPS) {
    cat("  SKIPPING: No PAM data and CREATE_SIGHTINGS_ONLY_MAPS is FALSE\n\n")
    return(NULL)
  }
  
  pam_ref_file <- NA_character_
  pam_ref_kde <- NULL
  if (has_pam_kde) {
    pam_ref_file <- pam_df$path[1]
    cat("  Using PAM reference KDE:", basename(pam_ref_file), "\n")
    pam_ref_kde <- load_kde_shapefile(pam_ref_file)
  } else {
    cat("  No PAM KDE found")
    if (CREATE_SIGHTINGS_ONLY_MAPS) {
      cat(" (will create sightings-only maps)\n")
    } else {
      cat("\n")
    }
  }
  
  MAX_SIGHTINGS_VARIANTS <- 12
  if (has_sightings_kde && nrow(sightings_df) > MAX_SIGHTINGS_VARIANTS) {
    sightings_df <- sightings_df[seq_len(MAX_SIGHTINGS_VARIANTS), , drop = FALSE]
    cat("  Note: capped sightings variants to", MAX_SIGHTINGS_VARIANTS, "\n")
  }
  
  # Filter source data
  if (is_grouped && !pam_only) {
    cat("  Filtering all baleen whale (mysticete) records...\n")
    species_sightings <- sightings_data %>% filter(cetacean_family == "mysticete")
    species_pam <- pam_data
  } else if (!is.null(species_info$is_grouped) && species_info$is_grouped && !pam_only) {
    cat("  Filtering all toothed whale (odontocete) records...\n")
    species_sightings <- sightings_data %>% filter(cetacean_family == "odontocete")
    species_pam <- data.frame()
  } else if (pam_only) {
    cat("  PAM-only species - no sightings data\n")
    species_sightings <- data.frame()
    
    if (isTRUE(is_grouped) || isTRUE(species_info$is_grouped)) {
      # Grouped PAM-only: keep all PAM stations passed into the function
      species_pam <- pam_data
    } else if (!is.null(species_info$pam_code)) {
      species_pam <- pam_data %>% dplyr::filter(species == species_info$pam_code)
    } else {
      species_pam <- data.frame()
    }
  } else {
    cat("  Looking for species with common_name:", species_name, "\n")
    
    # Try exact match first
    species_sightings <- sightings_data %>% filter(common_name == species_name)
    
    # Check for alternate names (e.g., Cuvier's for Goose Beaked)
    if (nrow(species_sightings) == 0 && !is.null(species_info$alternate_sightings_names)) {
      for (alt_name in species_info$alternate_sightings_names) {
        cat("  Trying alternate name:", alt_name, "\n")
        species_sightings <- sightings_data %>% filter(common_name == alt_name)
        if (nrow(species_sightings) > 0) {
          cat("  Found using alternate name\n")
          break
        }
      }
    }
    
    # Try alternative formats
    if (nrow(species_sightings) == 0) {
      cat("  No exact match, trying alternatives...\n")
      species_name_alt1 <- gsub("_", " ", species_name)
      species_sightings <- sightings_data %>% filter(common_name == species_name_alt1)
      if (nrow(species_sightings) > 0) cat("  Matched with spaces:", species_name_alt1, "\n")
    }
    
    if (nrow(species_sightings) == 0) {
      species_name_alt2 <- gsub("-", " ", species_name)
      species_name_alt2 <- gsub("_", " ", species_name_alt2)
      species_sightings <- sightings_data %>% filter(common_name == species_name_alt2)
      if (nrow(species_sightings) > 0) cat("  Matched with all spaces:", species_name_alt2, "\n")
    }
    
    if (nrow(species_sightings) == 0 && !pam_only) {
      cat("  WARNING: Could not match species in sightings data\n")
      cat("  Searched for common_name:", species_name, "\n\n")
    }
    
    # Get PAM data using PAM code
    if (!is.null(species_info$pam_code)) {
      species_pam <- pam_data %>% filter(species == species_info$pam_code)
    } else {
      species_pam <- data.frame()
    }
  }
  
  # Count deployments
  deployments <- if (nrow(species_pam) > 0 && "site" %in% names(species_pam)) {
    nrow(species_pam %>% group_by(site) %>% summarise(n()))
  } else if (nrow(species_pam) > 0 && "deployment" %in% names(species_pam)) {
    nrow(species_pam %>% group_by(deployment) %>% summarise(n()))
  } else {
    0
  }
  
  cat("  Sightings:", nrow(species_sightings), "\n")
  cat("  PAM deployments:", deployments, "\n")
  
  if (nrow(species_sightings) == 0 && deployments == 0) {
    cat("  ERROR: No data found for either sightings or PAM - skipping\n\n")
    return(NULL)
  }
  
  # Helper for tagging
  make_tag <- function(file) {
    if (is.na(file) || !file.exists(file)) return("none")
    p <- extract_kde_params(file)
    prox <- ifelse(isTRUE(p$has_proximity), "prox", "unw")
    meth <- ifelse(is.null(p$method) || is.na(p$method), "unknown", p$method)
    bw <- ifelse(is.na(p$bandwidth), "bwNA", paste0("bw", p$bandwidth))
    paste(prox, meth, bw, sep = "_")
  }
  
  pam_tag <- make_tag(pam_ref_file)
  
  if (!has_sightings_kde || pam_only) {
    sightings_files <- character(0)
  } else {
    sightings_files <- sightings_df$path
  }
  
  plot_list <- list()
  
  # Loop over sightings KDE variants
  n_variants <- max(1, length(sightings_files))
  
  for (i in seq_len(n_variants)) {
    
    if (length(sightings_files) > 0 && i <= length(sightings_files)) {
      sightings_file <- sightings_files[i]
      s_tag <- make_tag(sightings_file)
      cat("  Variant", i, "sightings KDE:", basename(sightings_file), "\n")
      sightings_kde <- load_kde_shapefile(sightings_file)
    } else {
      sightings_file <- NA_character_
      s_tag <- "none"
      sightings_kde <- NULL
    }
    
    pam_kde <- pam_ref_kde
    
    if (is.null(sightings_kde) && is.null(pam_kde)) {
      cat("  WARNING: could not load KDEs for variant", i, "\n")
      next
    }
    
    # Color palette
    n_quantiles_s <- if (!is.null(sightings_kde)) length(levels(sightings_kde$Quantile)) else 0
    n_quantiles_p <- if (!is.null(pam_kde)) length(levels(pam_kde$Quantile)) else 0
    n_quantiles <- max(n_quantiles_s, n_quantiles_p)
    pal <- rev(get(QUANTILE_PALETTE)(n_quantiles + 1))
    
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
    
    # Sightings map
    if (!is.null(sightings_kde) && nrow(species_sightings) > 0) {
      
      sight_dates <- as.Date(species_sightings$date_utc)
      sight_year_range <- paste(format(min(sight_dates), "%Y"), "to", format(max(sight_dates), "%Y"))
      
      sightings_params <- extract_kde_params(sightings_file)
      bw_text <- if (!is.na(sightings_params$bandwidth)) {
        paste0("BW: ", round(sightings_params$bandwidth/1000, 2), "km")
      } else {
        "BW: unknown"
      }
      weight_text <- ifelse(isTRUE(sightings_params$has_proximity), "Proximity weighted", "Unweighted")
      method_text <- ifelse(is.null(sightings_params$method), "unknown", sightings_params$method)
      param_text <- paste(bw_text, "|", method_text, "|", weight_text)
      
      p_sightings <- ggplot() +
        geom_sf(data = species_sightings, color = "grey40", fill = NULL, size = 1, shape = 21, alpha = 0.6) +
        
        geom_sf(data = sightings_kde, aes(fill = Quantile), color = NA, alpha = .6) +
        
        geom_sf(data = land, fill = "grey60", color = NA) +
        geom_sf(data = osw_wind, fill = NA, color = OSW_COLOR, alpha = OSW_ALPHA, linewidth = 1) +
        geom_sf(data = study_area, fill = NA, color = "black", linewidth = 0.8, linetype = "dashed") +
        scale_fill_manual(values = pal, name = "KDE Quantile", drop = FALSE, guide = "none") +
        coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
        labs(title = "Sightings",
             subtitle = paste0(nrow(species_sightings), " locations | ", param_text)) +
        base_theme +
        annotation_scale(location = "br", width_hint = 0.25)
      
    } else {
      sight_year_range <- "No data"
      p_sightings <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "No sightings data", size = 6) +
        theme_void()
    }
    
    # PAM map
    if (!is.null(pam_kde) && deployments > 0) {
      
      pam_zero <- species_pam %>% filter(proportion_det == 0)
      pam_detections <- species_pam %>% filter(proportion_det > 0)
      
      pam_params <- if (!is.na(pam_ref_file)) extract_kde_params(pam_ref_file) else list(bandwidth = NA_real_, method = "unknown")
      bw_text <- if (!is.na(pam_params$bandwidth)) paste0("BW: ", round(pam_params$bandwidth/1000, 2), "km") else "BW: unknown"
      method_text <- ifelse(is.null(pam_params$method), "unknown", pam_params$method)
      
      if (is_grouped) {
        subtitle_text <- paste0(deployments, " deployments with detections | ", bw_text, " | ", method_text,
                                " | Weighted by % days detected")
      } else {
        subtitle_text <- paste0(deployments, " deployments (", nrow(pam_zero), " with 0 detections) | ",
                                bw_text, " | ", method_text, " | Weighted by % days detected")
      }
      
      p_pam <- ggplot() +
        geom_sf(data = pam_kde, aes(fill = Quantile), color = NA, alpha = 0.6) +
        {if (nrow(pam_detections) > 0) geom_sf(data = pam_detections,
                                               aes(size = proportion_det, color = proportion_det),
                                               alpha = 0.8)} +
        geom_sf(data = land, fill = "grey60", color = NA) +
        geom_sf(data = osw_wind, fill = NA, color = OSW_COLOR, alpha = OSW_ALPHA, linewidth = 1)
      
      if (!is_grouped && nrow(pam_zero) > 0) {
        dummy_zero <- data.frame(x = NA, y = NA, label = "No Validated Detections")
        p_pam <- p_pam +
          geom_sf(data = pam_zero, color = "gray20", fill = "gray80", size = 1.5, shape = 21) +
          geom_point(data = dummy_zero, aes(x = x, y = y, shape = label),
                     color = "gray20", fill = "gray80", size = 1.5) +
          scale_shape_manual(name = "", values = 21,
                             labels = "No Validated Detections",
                             guide = guide_legend(order = 3, override.aes = list(size = 3)))
      }
      
      p_pam <- p_pam +
        geom_sf(data = study_area, fill = NA, color = "black", linewidth = 0.8, linetype = "dashed") +
        scale_fill_manual(values = pal, name = "KDE Quantile", drop = FALSE, guide = guide_legend(order = 1)) +
        scale_size_continuous(name = "Detection\nProportion",
                              range = c(1, 4), limits = c(0, 1),
                              breaks = c(0.2, 0.5, 0.9),
                              guide = guide_legend(order = 2)) +
        scale_color_gradient(low = "orange", high = "darkorange",
                             name = "Detection\nProportion",
                             limits = c(0, 1),
                             breaks = c(0.2, 0.5, 0.9),
                             guide = guide_legend(order = 2)) +
        coord_sf(xlim = xlims, ylim = ylims, crs = UTM20, expand = FALSE) +
        labs(title = "PAM", subtitle = subtitle_text) +
        base_theme +
        annotation_scale(location = "br", width_hint = 0.25)
      
    } else {
      p_pam <- NULL
    }
    
    # Combine plots
    if (!is.null(pam_kde) && deployments > 0 && !pam_only) {
      plot_title <- paste0(species_info$display_name, " 2015-2024 (PAM vs Sightings)")
      
      p_combined <- (p_sightings + p_pam) +
        plot_annotation(
          title = plot_title,
          theme = theme(plot.title = element_text(hjust = 0.5))
        )
      
      out_name <- paste0(
        "compare_",
        gsub("[^A-Za-z0-9]+", "_", tolower(species_key)),
        "__S_", s_tag,
        "__P_", pam_tag,
        ".png"
      )
      
      ggsave(file.path(COMPARE_OUTPUT_DIR, out_name),
             p_combined, width = 16, height = 8, dpi = 300, bg = "white")
      
    } else if (pam_only && !is.null(pam_kde) && deployments > 0) {
      plot_title <- paste0(species_info$display_name, " 2015-2024 (PAM Only)")
      
      p_combined <- p_pam +
        labs(title = plot_title) +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      
      out_name <- paste0(
        "pam_only_",
        gsub("[^A-Za-z0-9]+", "_", tolower(species_key)),
        "__P_", pam_tag,
        ".png"
      )
      
      ggsave(file.path(COMPARE_OUTPUT_DIR, out_name),
             p_combined, width = 10, height = 8, dpi = 300, bg = "white")
      
    } else {
      plot_title <- paste0(species_info$display_name, " 2015-2024 (Sightings)")
      
      p_combined <- p_sightings +
        labs(title = plot_title) +
        theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
      
      out_name <- paste0(
        "sightings_only_",
        gsub("[^A-Za-z0-9]+", "_", tolower(species_key)),
        "__S_", s_tag,
        ".png"
      )
      
      ggsave(file.path(COMPARE_OUTPUT_DIR, out_name),
             p_combined, width = 10, height = 8, dpi = 300, bg = "white")
    }
    
    cat("  Saved:", out_name, "\n\n")
    
    plot_list[[length(plot_list) + 1]] <- p_combined
  }
  
  if (length(plot_list) == 0) {
    cat("  WARNING: No comparison plots generated for", species_info$display_name, "\n\n")
    return(NULL)
  }
  
  if (length(plot_list) > 1) {
    return(patchwork::wrap_plots(plot_list, ncol = 1))
  } else {
    return(plot_list[[1]])
  }
}

# ==============================================================================
# HELPER FUNCTION: DISCOVER SIGHTINGS-ONLY SPECIES
# ==============================================================================

discover_sightings_only_species <- function(shapefile_dir = SHAPEFILE_DIR,
                                            exclude_species = character()) {
  
  cat("Discovering sightings-only species from KDE files...\n")
  
  all_kde_files <- list.files(shapefile_dir, pattern = "^KDE_.*\\.shp$",
                              full.names = FALSE, ignore.case = TRUE)
  
  # NEW: match your current naming scheme
  sightings_kdes <- all_kde_files[grepl("^KDE_sightings_species_proximity_", all_kde_files, ignore.case = TRUE)]
  
  if (length(sightings_kdes) == 0) {
    cat("  No sightings KDE files found\n")
    return(list())
  }
  
  species_list <- list()
  
  for (file in sightings_kdes) {
    
    temp <- sub("^KDE_sightings_species_proximity_", "", file, ignore.case = TRUE)
    temp <- sub("\\.shp$", "", temp, ignore.case = TRUE)
    
    # Strip method / bandwidth suffixes and anything after
    temp <- sub("_(fixed|adaptive|diggle).*", "", temp, ignore.case = TRUE)
    temp <- sub("_bw\\d+.*", "", temp, ignore.case = TRUE)
    
    species_base <- temp
    species_base <- sub("_+$", "", species_base)
    
    if (species_base %in% exclude_species) next
    
    if (!species_base %in% names(species_list)) {
      display_name <- species_base
      display_name <- gsub("^WHALE-", "", display_name)
      display_name <- gsub("^DOLPHINS-", "", display_name)
      display_name <- gsub("^PORPOISE-", "", display_name)
      display_name <- gsub("_", " ", display_name)
      display_name <- gsub("-", " ", display_name)
      display_name <- tools::toTitleCase(tolower(display_name))
      
      species_list[[species_base]] <- list(
        common_name = species_base,
        display_name = display_name,
        sightings_base = species_base,
        pam_base = NA_character_,
        is_sightings_only = TRUE,
        is_grouped = FALSE
      )
      
      cat("  Found:", display_name, "\n")
    }
  }
  
  cat("  Total sightings-only species found:", length(species_list), "\n\n")
  species_list
}


# ==============================================================================
# RUN COMPARISONS
# ==============================================================================

cat("=== CREATING PAM vs SIGHTINGS COMPARISON MAPS ===\n")
cat("CREATE_SIGHTINGS_ONLY_MAPS =", CREATE_SIGHTINGS_ONLY_MAPS, "\n\n")

# Process baleen whales
cat("Processing baleen whales (with PAM data)...\n\n")
for (species_key in names(BALEEN_SPECIES_MAPPING)) {
  species_info <- BALEEN_SPECIES_MAPPING[[species_key]]
  
  create_comparison_map(species_key, species_info, 
                        sightings_data, baleen_pam_data,
                        is_grouped = FALSE)
}

# Process beaked whales with sightings
cat("\n=== PROCESSING BEAKED WHALES WITH SIGHTINGS ===\n\n")
for (species_key in names(BEAKED_SPECIES_MAPPING)) {
  species_info <- BEAKED_SPECIES_MAPPING[[species_key]]
  
  if (species_info$has_sightings) {
    create_comparison_map(species_key, species_info, 
                          sightings_data, beaked_pam_data,
                          is_grouped = FALSE,
                          pam_only = FALSE)
  }
}

# Process beaked whales PAM-only
cat("\n=== PROCESSING PAM-ONLY BEAKED WHALES ===\n\n")
for (species_key in names(BEAKED_SPECIES_MAPPING)) {
  species_info <- BEAKED_SPECIES_MAPPING[[species_key]]
  
  if (!species_info$has_sightings) {
    create_comparison_map(species_key, species_info, 
                          sightings_data, beaked_pam_data,
                          is_grouped = FALSE,
                          pam_only = TRUE)
  }
}

# Process sightings-only species
if (CREATE_SIGHTINGS_ONLY_MAPS) {
  cat("\n=== PROCESSING SIGHTINGS-ONLY SPECIES ===\n\n")
  
  baleen_species_names <- sapply(BALEEN_SPECIES_MAPPING, function(x) x$sightings_base)
  beaked_with_sightings <- sapply(BEAKED_SPECIES_MAPPING[sapply(BEAKED_SPECIES_MAPPING, function(x) x$has_sightings)], 
                                  function(x) x$sightings_base)
  
  exclude_species <- c(baleen_species_names, beaked_with_sightings)
  
  sightings_only_species <- discover_sightings_only_species(exclude_species = exclude_species)
  
  if (length(sightings_only_species) > 0) {
    cat("Processing", length(sightings_only_species), "sightings-only species...\n\n")
    
    for (species_key in names(sightings_only_species)) {
      species_info <- sightings_only_species[[species_key]]
      
      create_comparison_map(species_key, species_info, 
                            sightings_data, pam_data,
                            is_grouped = FALSE)
    }
  }
}

cat("=== CREATING GROUPED ODONTOCETES SIGHTINGS MAP ===\n\n")
for (species_key in names(GROUPED_ODONTOCETES_SIGHTINGS)) {
  species_info <- GROUPED_ODONTOCETES_SIGHTINGS[[species_key]]
  create_comparison_map(species_key, species_info,
                        sightings_data, pam_data,
                        is_grouped = TRUE,
                        pam_only = FALSE)
}


# Grouped baleen
cat("=== CREATING GROUPED ALL BALEEN WHALES COMPARISON ===\n\n")
for (species_key in names(GROUPED_BALEEN)) {
  species_info <- GROUPED_BALEEN[[species_key]]
  
  create_comparison_map(species_key, species_info, 
                        sightings_data, baleen_pam_data,
                        is_grouped = TRUE)
}

# Grouped beaked
cat("=== CREATING GROUPED ALL BEAKED WHALES PAM MAP ===\n\n")
for (species_key in names(GROUPED_BEAKED)) {
  species_info <- GROUPED_BEAKED[[species_key]]
  
  create_comparison_map(species_key, species_info, 
                        sightings_data, beaked_pam_data,
                        is_grouped = TRUE,
                        pam_only = TRUE)
}

cat("\n=== COMPARISON MAPPING COMPLETE ===\n")
cat("Comparison maps saved to:", COMPARE_OUTPUT_DIR, "\n")