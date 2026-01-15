# look at season in baleen presence days

# Load necessary libraries
options(scipen = 999)

# libraries------
pacman::p_load(terra, spatstat, sf, dplyr, patchwork, ggplot2, lubridate, grafify)

filepath = "input/2025/baleen_presence_laura_2025.csv"

# projections------
UTM20 <- "EPSG:32620" # CODE FOR UTM Zone 20

# Read baleen whale PAM season data------
baleen_DOY <- read.csv(filepath)
baleen_DOY %>% group_by(species) %>% summarise(count = n())

whale_data <- baleen_DOY %>%
  mutate(UTC = as.POSIXct(rec_date, format = "%Y-%m-%d")) %>%
  mutate(Date = format(as_date(rec_date), "%Y-%m-%d")) %>%
  mutate(
    month = month(UTC),
    Season = case_when(
      month %in% c(12, 1, 2) ~ "Winter",
      month %in% c(3, 4, 5) ~ "Spring",
      month %in% c(6, 7, 8) ~ "Summer",
      month %in% c(9, 10, 11) ~ "Fall",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(site, Season) %>%
  mutate(
    season_days = n_distinct(Date)
  ) %>%
  ungroup()

baleen_PA_season <- whale_data %>%
  group_by(site, Season, species) %>%
  summarise(
    efort_days = first(season_days),
    detection_days = sum(presence > 0),
    proportion_det = detection_days / efort_days,
    latitude = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  )

# check data
season <- whale_data %>% group_by(Season, site) %>% summarise(days = n())

# Date range
whale_data %>% summarise(min(Date), max(Date))

# Seasonal Site Distribution and Species Analysis Plots
# Required: baleen_PA_season and whale_data from your existing code

# 1. SEASONAL SUMMARY STATISTICS ----

# Number of sites per season
sites_per_season <- whale_data %>%
  group_by(Season) %>%
  summarise(
    n_sites = n_distinct(site),
    lat_min = min(latitude),
    lat_max = max(latitude),
    lon_min = min(longitude),
    lon_max = max(longitude),
    lat_range = lat_max - lat_min,
    lon_range = lon_max - lon_min
  )

# Number of years per season
years_per_season <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(Season) %>%
  summarise(
    n_years = n_distinct(year),
    year_min = min(year),
    year_max = max(year)
  )

# Combine seasonal summaries
seasonal_summary <- left_join(sites_per_season, years_per_season, by = "Season")

print("Seasonal Summary:")
print(seasonal_summary)

# 2. PLOT: Number of Sites Sampled per Year by Season ----
sites_year_season <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(Season, year) %>%
  summarise(n_sites = n_distinct(site), .groups = "drop")

p1 <- ggplot(
  sites_year_season,
  aes(
    x = factor(year), y = n_sites,
    fill = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall"))
  )
) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_manual(
    values = c(
      "Winter" = "#2E86AB", "Spring" = "#06A77D",
      "Summer" = "#F18F01", "Fall" = "#A23B72"
    )
  ) +
  labs(
    title = "",
    x = "",
    y = "Number of Sites",
    fill = "Season"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 3. PLOT: Geographic Range per Season ----
range_data <- seasonal_summary %>%
  tidyr::pivot_longer(
    cols = c(lat_range, lon_range),
    names_to = "degrees",
    values_to = "range"
  ) %>%
  mutate(degrees = ifelse(degrees == "lat_range", "Latitude", "Longitude"))

p2 <- ggplot(range_data, aes(
  x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
  y = range, fill = degrees
)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_fill_manual(values = c("Latitude" = "#06A77D", "Longitude" = "#D4A373")) +
  labs(
    title = "Distribution of Sites by Season",
    x = "",
    y = "Range *",
    fill = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.minor = element_blank()
  )

# 4. SEASONAL SUMMARY BY SPECIES ----

species_season_summary <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(Season, species) %>%
  summarise(
    n_sites = n_distinct(site),
    n_years = n_distinct(year),
    detection_days = sum(presence > 0), # Days with actual detections
    effort_days = n(), # Total days sampled
    .groups = "drop"
  )

print("\nSpecies by Season Summary:")
print(species_season_summary)

# 5. PLOT: Sites with Detections per Season by Species ----
# Check if there's actual variability in site counts by species
site_species_check <- species_season_summary %>%
  select(Season, species, n_sites) %>%
  tidyr::pivot_wider(names_from = species, values_from = n_sites)

print("\nSites per season by species (check for variability):")
print(site_species_check)

# Better plot: Show proportion of sites with detections for each species
species_detection_prop <- baleen_PA_season %>%
  group_by(Season, species) %>%
  summarise(
    sites_with_detections = sum(detection_days > 0),
    total_sites = n(),
    prop_sites_detected = sites_with_detections / total_sites,
    .groups = "drop"
  )

p3 <- ggplot(
  species_detection_prop,
  aes(
    x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
    y = prop_sites_detected, fill = species
  )
) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1)) +
  labs(
    title = "% Sites with Detections ",
    x = "",
    y = "% Sites",
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# 6. PLOT: Detection Days by Season and Species ----
p4 <- ggplot(
  species_season_summary,
  aes(
    x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
    y = detection_days, fill = species
  )
) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(
    title = "Detection Days",
    x = "",
    y = "Detection Days",
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# 7. PLOT: Years Sampled by Season and Species ----
p5 <- ggplot(
  species_season_summary,
  aes(
    x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")),
    y = n_years, fill = species
  )
) +
  geom_col(position = "dodge", alpha = 0.8) +
  labs(
    title = "Sampling Years",
    x = "",
    y = "Number of Years",
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

# 8. COMBINED LAYOUT ----
# Overall seasonal patterns
top <- (p1 | p2)
bottom <- (p3 | p4) &
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom")

combined_seasonal <- top / bottom
combined_seasonal <- combined_seasonal +
  plot_annotation(
    title = "Seasonal Summary of Baleen Whale Detections",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_seasonal)
ggsave("output/Effort/PAM_Effort/combined_seasonal.png", combined_seasonal, h = 8, w = 11, units = "in")

# 9. Site Distribution Map by Season & Year ----
# FIXED VERSION WITH PROPER STUDY AREA AND OSW OVERLAY

# Load study area if not already loaded
if (!exists("study_area") && file.exists("shapefiles/study_area.shp")) {
  study_area <- st_read("shapefiles/study_area.shp", quiet = TRUE) %>%
    st_transform(UTM20)
}

# Load OSW wind areas if not already loaded
if (!exists("osw_wind") && file.exists("shapefiles/WEA/Designated_WEAs_25_07_29.shp")) {
  osw_wind <- st_read("shapefiles/WEA/Designated_WEAs_25_07_29.shp", quiet = TRUE) %>%
    st_transform(UTM20)
}

# Load land if not already loaded
if (!exists("land") && file.exists("shapefiles/Canada/NE_10mLand.shp")) {
  land <- st_read("shapefiles/Canada/NE_10mLand.shp", quiet = TRUE) %>%
    st_transform(UTM20)
} else if (exists("land") && !is.null(land)) {
  # Ensure land is in UTM20
  if (st_crs(land) != st_crs(UTM20)) {
    land <- st_transform(land, UTM20)
  }
}

# Prepare site locations and transform to UTM20
site_locations <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(site, Season) %>%
  summarise(
    latitude = first(latitude),
    longitude = first(longitude),
    n_species = n_distinct(species),
    n_years = n_distinct(year),
    .groups = "drop"
  ) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(UTM20)

# Get bounding box from study area if available, otherwise from data
if (exists("study_area") && !is.null(study_area)) {
  bbox <- st_bbox(study_area)
  xlims <- c(bbox["xmin"], bbox["xmax"])
  ylims <- c(bbox["ymin"], bbox["ymax"])
} else {
  # Fallback to data extent with buffer
  bbox <- st_bbox(site_locations)
  buffer <- 50000 # 50km buffer in meters (UTM)
  xlims <- c(bbox["xmin"] - buffer, bbox["xmax"] + buffer)
  ylims <- c(bbox["ymin"] - buffer, bbox["ymax"] + buffer)
}

p_map <- ggplot() +
  # Add land first
  {
    if (exists("land") && !is.null(land)) {
      geom_sf(data = land, fill = "grey60", color = NA)
    }
  } +
  # Add points
  geom_sf(
    data = site_locations,
    aes(color = Season, alpha = n_years),
    size = 3
  ) +
  # Add OSW areas
  {
    if (exists("osw_wind") && !is.null(osw_wind)) {
      geom_sf(
        data = osw_wind, fill = NA, color = "#FF6B35",
        linewidth = 1, inherit.aes = FALSE
      )
    }
  } +
  # Add study area outline
  {
    if (exists("study_area") && !is.null(study_area)) {
      geom_sf(
        data = study_area, fill = NA, color = "black",
        linewidth = 0.8, linetype = "dashed", inherit.aes = FALSE
      )
    }
  } +
  scale_alpha_continuous(range = c(0.3, 1), name = "Years Sampled") +
  scale_color_manual(
    values = c(
      "Winter" = "#2E86AB", "Spring" = "#06A77D",
      "Summer" = "#F18F01", "Fall" = "#A23B72"
    )
  ) +
  coord_sf(xlim = xlims, ylim = ylims, crs = st_crs(UTM20), expand = FALSE) +
  facet_wrap(~Season) +
  labs(
    title = "Site Level Annual Sample Effort",
    x = "",
    y = ""
  ) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

print(p_map)
ggsave("output/Effort/PAM_Effort/Site_effort_map.png", p_map, width = 12, height = 10, dpi = 300)

# ============================
# PAM COVERAGE: CONCENTRIC RANGE BUFFERS (single map)--------
# ============================

# 1) Build a station layer
pam_tbl <- whale_data

has_site <- any(names(pam_tbl) %in% c("site", "station", "Station", "Site"))
site_col <- intersect(names(pam_tbl), c("site", "station", "Station", "Site"))[1]

stations_sf <- if (has_site) {
  pam_tbl %>%
    mutate(site_id = .data[[site_col]]) %>%
    group_by(site_id) %>%
    summarise(
      longitude = first(longitude),
      latitude = first(latitude),
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(UTM20)
} else {
  pam_tbl %>%
    mutate(
      site_id = paste0(round(longitude, 4), "_", round(latitude, 4))
    ) %>%
    group_by(site_id) %>%
    summarise(
      longitude = first(longitude),
      latitude = first(latitude),
      .groups = "drop"
    ) %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(UTM20)
}

# Store total stations before clipping
total_stations <- nrow(stations_sf)

# Clip stations to study area BEFORE buffering
if (exists("study_area") && !is.null(study_area)) {
  stations_sf <- st_intersection(stations_sf, study_area)
}

study_area_stations <- nrow(stations_sf)

# SUMMARY STATISTICS FOR STUDY AREA STATIONS ----
cat("\n=== PAM STATION SUMMARY ===\n")
cat("Total PAM stations (all data):", total_stations, "\n")
cat("Stations within study area:", study_area_stations, "\n")
cat("Stations excluded:", total_stations - study_area_stations, "\n")
cat("Percentage in study area:", round(100 * study_area_stations / total_stations, 1), "%\n\n")

# Get date range and species from whale_data
date_summary <- whale_data %>%
  summarise(
    min_date = min(UTC, na.rm = TRUE),
    max_date = max(UTC, na.rm = TRUE),
    total_days = n_distinct(Date),
    n_species = n_distinct(species)
  )

species_list <- whale_data %>%
  distinct(species) %>%
  pull(species) %>%
  sort()

cat("Date range:", format(date_summary$min_date, "%Y-%m-%d"), "to", 
    format(date_summary$max_date, "%Y-%m-%d"), "\n")
cat("Total unique detection days:", date_summary$total_days, "\n")
cat("Species detected:", paste(species_list, collapse = ", "), "\n")

# Seasonal coverage summary
seasonal_station_summary <- whale_data %>%
  filter(site %in% stations_sf$site_id) %>%
  group_by(Season) %>%
  summarise(
    n_stations = n_distinct(site),
    n_detection_days = n_distinct(Date),
    n_years = n_distinct(year(UTC)),
    .groups = "drop"
  ) %>%
  arrange(match(Season, c("Winter", "Spring", "Summer", "Fall")))

cat("\n=== SEASONAL COVERAGE (Study Area Stations Only) ===\n")
print(seasonal_station_summary)

# Species-specific summary
species_station_summary <- whale_data %>%
  filter(site %in% stations_sf$site_id) %>%
  group_by(species) %>%
  summarise(
    n_stations = n_distinct(site),
    n_detection_days = sum(presence > 0),
    pct_stations = round(100 * n_distinct(site) / study_area_stations, 1),
    .groups = "drop"
  )

cat("\n=== SPECIES DETECTION SUMMARY (Study Area Stations Only) ===\n")
print(species_station_summary)
cat("\n")

# 2) Define range scenarios (km)-----
ranges_km <- c(50, 25, 10, 5) # plot big-to-small so circles look concentric
ranges_m <- ranges_km * 1000

# 3) Create buffers (after clipping stations)
buffers_sf <- purrr::map2_dfr(
  ranges_km, ranges_m,
  ~ {
    st_buffer(stations_sf, dist = .y) %>%
      mutate(range_km = .x)
  }
)

# 4) Style: blue gradient by range (darker near station)
# Need to match the order we'll use in the legend
range_fill <- c(
  "50" = scales::alpha("#6BAED6", 0.15),
  "25" = scales::alpha("#2171B5", 0.20),
  "10" = scales::alpha("#08519C", 0.30),
  "5" = scales::alpha("#08306B", 0.35)
)

buffers_sf$range_km <- factor(buffers_sf$range_km, levels = c(50, 25, 10, 5))

# 5) Get proper extent for the map -----
# Use study area bounds if available, otherwise use buffer extent
if (exists("study_area") && !is.null(study_area)) {
  map_bbox <- st_bbox(study_area)
  map_xlims <- c(map_bbox["xmin"], map_bbox["xmax"])
  map_ylims <- c(map_bbox["ymin"], map_bbox["ymax"])
  
  # Crop land to bbox (not study area polygon since study area is marine)
  if (exists("land") && !is.null(land)) {
    # Add buffer to bbox to ensure land edges are included
    bbox_buffered <- st_as_sfc(map_bbox) %>% 
      st_buffer(10000) %>%  # 10km buffer
      st_bbox()
    land_cropped <- suppressWarnings(st_crop(land, bbox_buffered))
  }
} else {
  # Use the outermost buffer extent
  map_bbox <- st_bbox(buffers_sf)
  map_xlims <- c(map_bbox["xmin"], map_bbox["xmax"])
  map_ylims <- c(map_bbox["ymin"], map_bbox["ymax"])
  
  if (exists("land") && !is.null(land)) {
    bbox_buffered <- st_as_sfc(map_bbox) %>% 
      st_buffer(10000) %>%
      st_bbox()
    land_cropped <- suppressWarnings(st_crop(land, bbox_buffered))
  }
}

# 6) Plot---
p_pam_range <- ggplot() +
  # Land first
  {
    if (exists("land_cropped") && !is.null(land_cropped)) {
      geom_sf(data = land_cropped, fill = "grey75", color = "grey50", 
              linewidth = 0.3, inherit.aes = FALSE)
    }
  } +
  # Range buffers
  geom_sf(
    data = buffers_sf,
    aes(fill = as.character(range_km)),
    color = NA,
    inherit.aes = FALSE
  ) +
  # Stations on top
  geom_sf(data = stations_sf, shape = 21, fill = "white", color = "darkblue", 
          size = 2, stroke = 0.8, inherit.aes = FALSE) +
  # OWA outline
  {
    if (exists("osw_wind") && !is.null(osw_wind)) {
      geom_sf(data = osw_wind, fill = NA, color = "#FF6B35", 
              linewidth = 1, inherit.aes = FALSE)
    }
  } +
  # Study area outline
  {
    if (exists("study_area") && !is.null(study_area)) {
      geom_sf(data = study_area, fill = NA, color = "black", 
              linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE)
    }
  } +
  scale_fill_manual(
    values = range_fill,
    name = "Range from \nRecorder (Km)",
    breaks = c("10", "25", "50"),
    labels = c("10", "25", "50")
  ) +
  guides(fill = guide_legend(
    override.aes = list(alpha = 0.5), 
    nrow = 1,
    label.position = "bottom",
    title.position = "top"
  )) +
  coord_sf(xlim = map_xlims, ylim = map_ylims, crs = st_crs(UTM20), expand = FALSE) +
  ggspatial::annotation_scale(location = "br", width_hint = 0.25) +
  labs(
    title = "Potential PAM station coverage scenarios",
    subtitle = "Actual detectability varies with species, call type, noise, and propagation."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = c(0.98, 0.058),
    legend.key.spacing.x = unit(.02, "cm"),
    legend.key.justification = "center",
    legend.justification = c(1, 0),
    legend.background = element_blank(),
    legend.box.background = element_rect(fill = "white", colour = NA),
    legend.key.size = unit(0.58, "cm"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    plot.subtitle = element_text(size = 12, color = "grey30", margin = margin(b = 10))
  )

ggsave(
  filename = "output/Effort/PAM_Effort/PAM_station_range_scenarios.png",
  plot = p_pam_range,
  width = 10, height = 8, dpi = 300
)

print(p_pam_range)