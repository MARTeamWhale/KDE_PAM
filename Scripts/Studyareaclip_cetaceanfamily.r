#depends on dedup sightings

#will give some additional infos on sightings data in prep for KDE analysis:

# The script will print:
#   
      # Complete species counts for all species
      # List of species meeting 25 record threshold
      # List of species meeting 50 record threshold
      # Breakdown by cetacean family (odontocetes vs mysticetes vs other)
      # Species lists for each family



# ==============================================================================
# CLIP SIGHTINGS DATA TO STUDY AREA------
# ==============================================================================
# Clips deduplicated sightings to study area boundary before analysis
# Adds spatial clipping step before species summaries and family classification

library(dplyr)
library(readr)
library(sf)

# ==============================================================================
# CONFIGURATION------
# ==============================================================================

# Input paths
DEDUP_DATA_PATH <- "output/combined_dedup_1km_day.csv"
STUDY_AREA_PATH <- "shapefiles/studyArea/ESS_study_area_simple.shp"

# Output paths
CLIPPED_OUTPUT <- "output/combined_dedup_1km_day_clipped.csv"

# ==============================================================================
# 1. LOAD DATA-------
# ==============================================================================

message("Loading deduplicated sightings data...")
sightings_dedup <- read_csv(DEDUP_DATA_PATH, show_col_types = FALSE)
message(sprintf("  ✓ Loaded %d records", nrow(sightings_dedup)))

message("Loading study area shapefile...")
study_area <- st_read(STUDY_AREA_PATH, quiet = TRUE)
message(sprintf("  ✓ Loaded study area (CRS: %s)", st_crs(study_area)$input))

# ==============================================================================
# 2. CONVERT SIGHTINGS TO SPATIAL AND CLIP------
# ==============================================================================

message("\nClipping sightings to study area...")

# Convert sightings to sf object
sightings_sf <- st_as_sf(
  sightings_dedup,
  coords = c("lon", "lat"),
  crs = 4326,
  remove = FALSE
)

# Transform study area to match sightings CRS if needed
if (st_crs(study_area) != st_crs(sightings_sf)) {
  message("  Transforming study area to WGS84...")
  study_area <- st_transform(study_area, 4326)
}

# Spatial clip using st_intersection
sightings_clipped_sf <- st_intersection(sightings_sf, study_area)

# Convert back to regular data frame (drop geometry)
sightings_clipped <- sightings_clipped_sf %>%
  st_drop_geometry()

records_before <- nrow(sightings_dedup)
records_after <- nrow(sightings_clipped)
records_removed <- records_before - records_after

message(sprintf("  ✓ Clipped: %d → %d records", records_before, records_after))
message(sprintf("  ✗ Removed %d records outside study area (%.1f%%)", 
                records_removed, 100 * records_removed / records_before))

# Save clipped data
write_csv(sightings_clipped, CLIPPED_OUTPUT)
message(sprintf("\n✓ Saved clipped dataset: %s", CLIPPED_OUTPUT))

# ==============================================================================
# 3. SPECIES SUMMARIES AND KDE THRESHOLD ANALYSIS------
# ==============================================================================

message("\n==============================================================================")
message("SPECIES ANALYSIS (CLIPPED DATA)")
message("==============================================================================")

# Count records per species
species_counts <- sightings_clipped %>%
  group_by(species_cd, common_name, scientific_name) %>%
  summarise(
    n_records = n(),
    n_aerial = sum(source == "aerial"),
    n_wsdb = sum(source == "wsdb"),
    date_range = paste(min(date_utc), "to", max(date_utc)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_records))

message("\nSpecies record counts (all species):")
print(species_counts, n = Inf)

# Save species counts
write_csv(species_counts, "output/species_record_counts_clipped.csv")
message("\n✓ Saved: species_record_counts_clipped.csv")

# KDE threshold analysis
threshold_50 <- species_counts %>% filter(n_records >= 50)

message("\n==============================================================================")
message("KDE THRESHOLD ANALYSIS")
message("==============================================================================")

message(sprintf("\nSpecies meeting threshold of 50+ records: %d species", nrow(threshold_50)))
print(threshold_50 %>% select(species_cd, common_name, n_records), n = Inf)

# Save threshold summaries
write_csv(threshold_50, "output/species_kde_threshold_50_clipped.csv")
message("\n✓ Saved threshold summary")

# ==============================================================================
# 4. CETACEAN FAMILY CLASSIFICATION-------
# ==============================================================================

message("\n==============================================================================")
message("CETACEAN FAMILY CLASSIFICATION")
message("==============================================================================")

# Define odontocetes (toothed whales) and mysticetes (baleen whales)-----
odontocetes <- c(
  "7020",  # Sperm Whale
  "7019",  # Pygmy Sperm Whale
  "921",   # Atlantic Pilot Whale
  "7031",  # Long-finned Pilot Whale
  "922",   # Northern Bottlenose Whale
  "923",   # Sowerby's Beaked Whale
  "924",   # Beaked Whale (NS)
  "925",   # Cuvier's Beaked Whale
  "7039",  # Blainville's Beaked Whale
  "7041",  # Mesoplodont (NS)
  "931",   # Atlantic Bottlenose Dolphin
  "932",   # White-beaked Dolphin
  "933",   # Atlantic White-sided Dolphin
  "934",   # Common Dolphin
  "935",   # Risso's Dolphin
  "936",   # Striped Dolphin
  "937",   # Atlantic Spotted Dolphin
  "938",   # Fraser's Dolphin
  "7034",  # Pacific White-sided Dolphin
  "7025",  # Harbour Porpoise
  "7035",  # Dall's Porpoise
  "7028",  # Killer Whale
  "7029",  # Beluga
  "7037",  # False Killer Whale
  "7038",  # Spinner Dolphin
  "930"   # Dolphins/Porpoise (NS)

)

mysticetes <- c(
  "7021",  # Fin Whale
  "7022",  # Minke Whale
  "7023",  # North Atlantic Right Whale
  "7024",  # Humpback Whale
  "7026",  # Blue Whale
  "7027",  # Sei Whale
  "7030",  # Baleen Whale (NS)
  "7032",  # Bowhead
  "7033",  # Grey Whale
  "7040"   # Fin/Sei
)

# Add family classification-------
sightings_classified <- sightings_clipped %>%
  mutate(
    cetacean_family = case_when(
      species_cd %in% odontocetes ~ "odontocete",
      species_cd %in% mysticetes ~ "mysticete",
      TRUE ~ "other"  # seals, turtles, fish, birds
    )
  )

# Summary by family------
family_summary <- sightings_classified %>%
  group_by(cetacean_family) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species_cd),
    .groups = "drop"
  )

message("\nRecords by cetacean family:")
print(family_summary)

# Odontocete species summary-----
odontocete_summary <- sightings_classified %>%
  filter(cetacean_family == "odontocete") %>%
  group_by(species_cd, common_name) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  arrange(desc(n_records))

message("\nOdontocete (toothed whale) species:")
print(odontocete_summary, n = Inf)

# Mysticete species summary----
mysticete_summary <- sightings_classified %>%
  filter(cetacean_family == "mysticete") %>%
  group_by(species_cd, common_name) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  arrange(desc(n_records))

message("\nMysticete (baleen whale) species:")
print(mysticete_summary, n = Inf)

# Save classified dataset------
write_csv(sightings_classified, "output/combined/combined_dedup_1km_day_clipped_classified.csv")
message("\n✓ Saved classified dataset: combined_dedup_1km_day_clipped_classified.csv")

# Save family-specific datasets for KDE analysis
odontocete_data <- sightings_classified %>% filter(cetacean_family == "odontocete")
mysticete_data <- sightings_classified %>% filter(cetacean_family == "mysticete")

write_csv(odontocete_data, "input/2025/combined_sights/odontocete_sightings_dedup_clipped.csv")
write_csv(mysticete_data, "input/2025/combined_sights/mysticete_sightings_dedup_clipped.csv")

message(sprintf("\n✓ Saved odontocete dataset: %d records", nrow(odontocete_data)))
message(sprintf("✓ Saved mysticete dataset: %d records", nrow(mysticete_data)))

# Summary statistics by source
source_summary <- sightings_clipped %>%
  group_by(source) %>%
  summarise(
    n_records = n(),
    n_species = n_distinct(species_cd),
    date_range = paste(min(date_utc), "to", max(date_utc))
  )

message("\n==============================================================================")
message("FINAL SUMMARY BY DATA SOURCE (CLIPPED)")
message("==============================================================================")
print(source_summary)
message("==============================================================================")

# Create simple comparison summary
comparison <- data.frame(
  dataset = c("Original", "Clipped", "Removed"),
  n_records = c(records_before, records_after, records_removed),
  percentage = c(100, 100 * records_after / records_before, 
                 100 * records_removed / records_before)
)

message("\n==============================================================================")
message("CLIPPING SUMMARY")
message("==============================================================================")
print(comparison)
message("==============================================================================")