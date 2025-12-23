#look at season in baleen presence days

# Load necessary libraries

options(scipen = 999)
#libraries------
pacman::p_load(terra,spatstat,sf,dplyr,patchwork, ggplot2, lubridate, grafify)

filepath = "input/2025/baleen_presence_laura_2025.csv"

#projections------
UTM20 <- "EPSG:32620" # CODE FOR UTM Zone 20

# Read baleen whale PAM season data------
baleen_DOY <- read.csv(filepath)
baleen_DOY%>%group_by(species)%>%summarise(count = n())

whale_data <-baleen_DOY%>%mutate(UTC = as.POSIXct(rec_date, format = "%Y-%m-%d"))%>%
  mutate(Date = format(as_date(rec_date), "%Y-%m-%d"))%>%
           mutate(month = month(UTC),  
                  Season = case_when(
               month %in%  c(12, 1, 2) ~ "Winter",
               month %in%  c(3, 4, 5) ~ "Spring",
               month %in%  c(6, 7, 8) ~ "Summer",
               month %in%  c(9,10,11) ~ "Fall",
               TRUE ~ NA_character_
             )
           )%>%group_by(site, Season) %>% 
  mutate(
    season_days = n_distinct(Date)
  ) %>%
  ungroup()

baleen_PA_season = whale_data%>%
  group_by(site, Season, species) %>%
  summarise(
    efort_days = first(season_days),
    detection_days = sum(presence > 0),
    proportion_det = detection_days / efort_days,
    latitude  = first(latitude),
    longitude = first(longitude),
    .groups = "drop"
  )


#check data
season = whale_data%>%group_by(Season, site, )%>%summarise(days = n())

# Date range
whale_data%>%summarise(min(Date), max(Date))

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

p1 <- ggplot(sites_year_season, 
             aes(x = factor(year), y = n_sites, 
                 fill = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")))) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_manual(
    values = c("Winter" = "#2E86AB", "Spring" = "#06A77D", 
               "Summer" = "#F18F01", "Fall" = "#A23B72")
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

p2 <- ggplot(range_data, aes(x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")), 
                             y = range, fill = degrees)) +
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
    detection_days = sum(presence > 0),  # Days with actual detections
    effort_days = n(),  # Total days sampled
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

p3 <- ggplot(species_detection_prop, 
             aes(x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")), 
                 y = prop_sites_detected, fill = species)) +
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
p4 <- ggplot(species_season_summary, 
             aes(x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")), 
                 y = detection_days, fill = species)) +
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
p5 <- ggplot(species_season_summary, 
             aes(x = factor(Season, levels = c("Winter", "Spring", "Summer", "Fall")), 
                 y = n_years, fill = species)) +
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
top  <- (p1 | p2)
bottom <- (p3 | p4)&
  plot_layout(guides = "collect") +
  theme(legend.position = "bottom")

combined_seasonal <- top / bottom
combined_seasonal <- combined_seasonal + 
  plot_annotation(
    title = "Seasonal Summary of Baleen Whale Detections",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_seasonal)
ggsave( "output/FIGS/Effort_plots/combined_seasonal.png", combined_seasonal, h = 8, w = 11, units = "in")

# # Species-specific patterns
# species_plot <- (p3 / p4 / p5) & 
#   plot_layout(guides = "collect") +
#   theme(legend.position = "bottom")
# 
# species_plot <- species_plot + 
#   plot_annotation(
#     title = "Species-Specific Seasonal Patterns",
#     theme = theme(plot.title = element_text(size = 16, face = "bold"))
#   )
# 
# 
# print(species_plot)

# 9. Site Distribution Map by Season & Year ----
site_locations <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(site, Season) %>%
  summarise(
    latitude = first(latitude),
    longitude = first(longitude),
    n_species = n_distinct(species),
    n_years = n_distinct(year),
    .groups = "drop"
  )

# Get coastline data
land <- rnaturalearth::ne_coastline(scale = "medium", returnclass = "sf")

# Create bounding box for the map based on your data
bbox <- st_bbox(c(xmin = min(site_locations$longitude) - 1, 
                  xmax = max(site_locations$longitude) + 1,
                  ymin = min(site_locations$latitude) - 1, 
                  ymax = max(site_locations$latitude) + 1),
                crs = st_crs(4326))

p_map <- ggplot() +
  geom_sf(data = land, fill = "gray90", color = "gray50", size = 0.3) +
  geom_point(data = site_locations, 
             aes(x = longitude, y = latitude, color = Season, alpha = n_years),
             size = 3) +
  scale_alpha_continuous(range = c(0.3, 1), name = "Years Sampled") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  facet_wrap(~Season) +
  labs(
    title = "Site Level Annual Sample Effort",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")+
  scale_x_continuous(breaks = scales::pretty_breaks(4))


print(p_map)
ggsave( "output/FIGS/Effort_plots/Site_effort map.png", p_map)
