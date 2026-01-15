# ==============================================================================
# MARINE SPECIES SIGHTINGS DATA PROCESSOR
# ==============================================================================
# Purpose: Combine aerial and WSDB sightings with standardized schema and deduplication
# 
# Deduplication rule: Maximum 1 sighting per species per day per 1km grid cell
# 
# Key outputs:
#   1. Combined dataset with standardized fields
#   2. Deduplicated dataset 
#   3. Deduplication summary statistics
# ==============================================================================

# Load required packages
pacman::p_load(dplyr, stringr, lubridate, readr, sf, tidyr)

# ==============================================================================
# CONFIGURATION
# ==============================================================================

# Input file paths (update as needed)
AERIAL_PATH <- "input/2025/aerial/deduplicated_sightings.csv"
WSDB_PATH   <- "input/2025/sights/WSDB_DFO_Sep2025.csv"
SPECIES_CODES_PATH <- "input/speciesCodes.csv"
SPECIES_LOOKUP_PATH <- "input/species_lookup_aerial_to_wsdb.csv"

# Output paths
COMBINED_OUTPUT     <- "output/combined_sightings_standardized.csv"
DEDUP_OUTPUT        <- "output/combined_dedup_1km_day.csv"
SUMMARY_OUTPUT      <- "output/deduplication_summary.csv"

# Spatial configuration
# EPSG:32620 = UTM Zone 20N (appropriate for Atlantic Canada/Northeast US)
PROJ_CRS <- 32620
CELL_SIZE_M <- 1000  # 1 km grid cells

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Parse WSDB datetime format (handles unusual formatting)
parse_wsdb_datetime <- function(x) {
  # WSDB format: "2000-01-12 18:00.00.000000000"
  # Convert "HH:MM.SS.nnn" to "HH:MM:SS"
  x_clean <- str_replace(as.character(x), "(\\d{2}:\\d{2})\\.(\\d{2}).*$", "\\1:\\2")
  suppressWarnings(ymd_hms(x_clean, tz = "UTC"))
}

#' Parse time in HHMM format to HH:MM:SS
parse_hhmm_time <- function(hhmm) {
  x <- str_pad(str_trim(as.character(hhmm)), width = 4, side = "left", pad = "0")
  x[x == "NA" | is.na(x)] <- NA_character_
  
  hh <- str_sub(x, 1, 2)
  mm <- str_sub(x, 3, 4)
  out <- paste0(hh, ":", mm, ":00")
  out[is.na(x)] <- NA_character_
  
  return(out)
}

#' Create 1km grid cell ID from projected coordinates
create_grid_cell_id <- function(x_m, y_m, cell_size = CELL_SIZE_M) {
  paste0(floor(x_m / cell_size), "_", floor(y_m / cell_size))
}

# ==============================================================================
# 1. READ REFERENCE DATA---------
# ==============================================================================

message("Loading reference data...")

# Load WSDB species code reference
wsdb_species_ref <- read_csv(SPECIES_CODES_PATH, show_col_types = FALSE) %>%
  rename_with(tolower) %>%
  rename_with(~str_replace_all(., "commonname", "common_name")) %>%
  rename_with(~str_replace_all(., "scientif", "scientific_name")) %>%
  mutate(species_cd = as.character(species_cd))

# Load species lookup (aerial -> WSDB mapping)
species_lookup_raw <- read_csv(SPECIES_LOOKUP_PATH, show_col_types = FALSE)

species_lookup <- species_lookup_raw %>%
  rename_with(~str_replace(., "sp_code_aerial", "sp_code_raw")) %>%
  select(sp_code_raw, wsdb_species_cd) %>%
  mutate(
    sp_code_raw = as.character(sp_code_raw),
    wsdb_species_cd = as.character(wsdb_species_cd)
  ) %>%
  filter(!is.na(wsdb_species_cd), wsdb_species_cd != "", wsdb_species_cd != "NA")

# Check if we have valid mappings
if (nrow(species_lookup) == 0) {
  stop(paste0(
    "No valid species mappings found in: ", SPECIES_LOOKUP_PATH, "\n",
    "Please ensure the wsdb_species_cd column is filled with valid codes.\n",
    "Run the quick fix script to create a properly formatted lookup file."
  ))
}

message(sprintf("  ✓ Loaded %d species mappings", nrow(species_lookup)))
# ==============================================================================
# 2. READ AND STANDARDIZE AERIAL DATA----------
# ==============================================================================

message("Reading aerial sightings data...")

aerial_raw <- read_csv(
  AERIAL_PATH,
  locale = locale(encoding = "Windows-1252"),
  show_col_types = FALSE
)

# Standardize to common schema
aerial_std <- aerial_raw %>%
  transmute(
    source = "aerial",
    source_id = coalesce(as.character(ID), as.character(`...1`)),
    sp_code_raw = as.character(Sp_code),
    sp_name_raw = as.character(Sp_name_clean),
    lon = as.numeric(Long_sight),
    lat = as.numeric(Lat_sight),
    datetime_utc = suppressWarnings(ymd_hms(Date_time_utc, tz = "UTC"))
  ) %>%
  filter(!is.na(lon), !is.na(lat), !is.na(datetime_utc), !is.na(sp_code_raw))

message(sprintf("  ✓ Loaded %d aerial records", nrow(aerial_std)))

# Map aerial codes to WSDB codes
aerial_mapped <- aerial_std %>%
  left_join(species_lookup, by = "sp_code_raw") %>%
  mutate(species_cd = wsdb_species_cd) %>%
  filter(!is.na(species_cd)) %>%
  select(source, source_id, species_cd, lon, lat, datetime_utc)

message(sprintf("  ✓ Mapped %d aerial records (excluded %d unmapped)", 
                nrow(aerial_mapped), nrow(aerial_std) - nrow(aerial_mapped)))

# ==============================================================================
# 3. READ AND STANDARDIZE WSDB DATA-------
# ==============================================================================

message("Reading WSDB sightings data...")

# Read with explicit column types to avoid parsing issues
wsdb_raw <- read_csv(
  WSDB_PATH, 
  show_col_types = FALSE,
  col_types = cols(.default = col_character())
)

# Standardize to common schema
wsdb_std <- wsdb_raw %>%
  transmute(
    source = "wsdb",
    source_id = as.character(WS_ID),
    species_cd = as.character(SPECIES_CD),
    lon = as.numeric(LONGITUDE),
    lat = as.numeric(LATITUDE),
    datetime_utc = parse_wsdb_datetime(WS_DATETIME_UTC),
    date_local = suppressWarnings(ymd(WS_DATE, tz = "UTC")),
    time_local = parse_hhmm_time(WS_TIME)
  ) %>%
  # Fall back to date + time if datetime_utc is missing
  mutate(
    datetime_utc = if_else(
      is.na(datetime_utc) & !is.na(date_local) & !is.na(time_local),
      suppressWarnings(ymd_hms(paste(date_local, time_local), tz = "UTC")),
      datetime_utc
    )
  ) %>%
  select(source, source_id, species_cd, lon, lat, datetime_utc) %>%
  filter(!is.na(lon), !is.na(lat), !is.na(datetime_utc), !is.na(species_cd))

message(sprintf("  ✓ Loaded %d WSDB records", nrow(wsdb_std)))

# ==============================================================================
# 4. COMBINE DATASETS---------
# ==============================================================================

message("Combining datasets...")

sightings_combined <- bind_rows(wsdb_std, aerial_mapped) %>%
  mutate(
    datetime_utc = as.POSIXct(datetime_utc, tz = "UTC"),
    date_utc = as.Date(datetime_utc, tz = "UTC")
  ) %>%
  # Add human-readable species info from reference
  left_join(
    wsdb_species_ref %>% select(species_cd, common_name, scientific_name),
    by = "species_cd"
  )

message(sprintf("  ✓ Combined dataset: %d total records", nrow(sightings_combined)))

# Save combined standardized dataset
write_csv(sightings_combined, COMBINED_OUTPUT)
message(sprintf("  ✓ Saved: %s", COMBINED_OUTPUT))

# ==============================================================================
# 5. CREATE SPATIAL GRID CELLS (1 KM)------------
# ==============================================================================

message("Creating 1km grid cells...")

# Convert to spatial object and project to UTM
sightings_sf <- st_as_sf(
  sightings_combined,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
) %>%
  st_transform(PROJ_CRS)

# Extract projected coordinates and create grid cells
xy_coords <- st_coordinates(sightings_sf)

sightings_sf <- sightings_sf %>%
  mutate(
    x_m = xy_coords[, 1],
    y_m = xy_coords[, 2],
    grid_cell_1km = create_grid_cell_id(x_m, y_m)
  )

message(sprintf("  ✓ Created grid cells across %d unique cells", 
                n_distinct(sightings_sf$grid_cell_1km)))

# ==============================================================================
# 6. DEDUPLICATION------------
# ==============================================================================

message("Applying deduplication rules...")
message("  Rule: Max 1 sighting per species per day per 1km cell")

# Create deduplication key for auditability
sightings_with_key <- sightings_sf %>%
  st_drop_geometry() %>%
  mutate(
    dedupe_key = paste(species_cd, date_utc, grid_cell_1km, sep = "|")
  )

pre_dedup_count <- nrow(sightings_with_key)

# Deduplicate: keep earliest datetime in each group
sightings_dedup <- sightings_with_key %>%
  arrange(species_cd, date_utc, grid_cell_1km, datetime_utc) %>%
  group_by(species_cd, date_utc, grid_cell_1km) %>%
  mutate(
    records_in_group = n(),
    is_duplicate = row_number() > 1
  ) %>%
  slice(1) %>%  # Keep only first record
  ungroup()

post_dedup_count <- nrow(sightings_dedup)
records_removed <- pre_dedup_count - post_dedup_count

message(sprintf("  ✓ Deduplicated: %d → %d records", pre_dedup_count, post_dedup_count))
message(sprintf("  ✗ Removed: %d duplicate records (%.1f%%)", 
                records_removed, 100 * records_removed / pre_dedup_count))

# ==============================================================================
# 7. DEDUPLICATION SUMMARY---------
# ==============================================================================

dedup_summary <- sightings_with_key %>%
  group_by(species_cd, date_utc, grid_cell_1km) %>%
  summarise(
    n_records = n(),
    sources = paste(unique(source), collapse = ", "),
    .groups = "drop"
  ) %>%
  filter(n_records > 1) %>%
  left_join(
    wsdb_species_ref %>% select(species_cd, common_name),
    by = "species_cd"
  ) %>%
  arrange(desc(n_records))

if (nrow(dedup_summary) > 0) {
  message(sprintf("\nDuplicate groups found: %d", nrow(dedup_summary)))
  message("Top duplicate groups:")
  print(head(dedup_summary, 10))
  
  write_csv(dedup_summary, SUMMARY_OUTPUT)
  message(sprintf("  ✓ Saved deduplication summary: %s", SUMMARY_OUTPUT))
}

# ==============================================================================
# 8. SAVE FINAL OUTPUT-----------
# ==============================================================================

write_csv(sightings_dedup, DEDUP_OUTPUT)
message(sprintf("\n✓ COMPLETE - Saved deduplicated dataset: %s", DEDUP_OUTPUT))

# Summary statistics by source
source_summary <- sightings_dedup %>%
  group_by(source) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species_cd),
    date_range = paste(min(date_utc), "to", max(date_utc))
  )

message("\n==============================================================================")
message("FINAL SUMMARY")
message("==============================================================================")
print(source_summary)
message("==============================================================================")