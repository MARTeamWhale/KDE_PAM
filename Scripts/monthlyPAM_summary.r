# Monthly Whale Detection Analysis Plots
# Required: whale_data and baleen_PA_season from SeasonPam.r

# Load necessary libraries

options(scipen = 999)
#libraries------
pacman::p_load(terra,spatstat,sf,dplyr,patchwork, ggplot2, lubridate)
filepath = "input/2025/baleen_presence_laura_2025.csv"


# 1. MONTHLY SUMMARY STATISTICS ----

# Number of sites per month
sites_per_month <- whale_data %>%
  group_by(month) %>%
  summarise(
    n_sites = n_distinct(site),
    lat_min = min(latitude),
    lat_max = max(latitude),
    lon_min = min(longitude),
    lon_max = max(longitude),
    lat_range = lat_max - lat_min,
    lon_range = lon_max - lon_min
  )

# Number of years per month
years_per_month <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(month) %>%
  summarise(
    n_years = n_distinct(year),
    year_min = min(year),
    year_max = max(year)
  )

# Combine monthly summaries
monthly_summary <- left_join(sites_per_month, years_per_month, by = "month")

print("Monthly Summary:")
print(monthly_summary)

# 2. PLOT: Number of Sites Sampled per Year by Month ----
sites_year_month <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(month, year) %>%
  summarise(n_sites = n_distinct(site), .groups = "drop")

# Create month labels
month_labels <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                  "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

p1 <- ggplot(sites_year_month, 
             aes(x = factor(year), y = n_sites, fill = factor(month))) +
  geom_col(position = "stack", alpha = 0.8) +
  scale_fill_viridis_d(option = "turbo", labels = month_labels) +
  labs(
    title = "Site (n) per Year by Month",
    x = "",
    y = "Number of Sites",
    fill = "Month"
  ) +
  theme_grafify(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )+
  scale_y_continuous(breaks = scales::pretty_breaks(6))

p1 

ggsave( "output/FIGS/Effort_plots/annual_siteffort_bymonth.png", p1)


# 3. PLOT: Sites per Month (overall) ----
p2 <- ggplot(monthly_summary, 
             aes(x = factor(month), y = n_sites)) +
  geom_col(fill = "#2E86AB", alpha = 0.8) +
  geom_text(aes(label = n_sites), vjust = -0.5, size = 3) +
  scale_x_discrete(labels = month_labels) +
  labs(
    title = "Total Sites Sampled by Month",
    x = "Month",
    y = "Number of Sites"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

# 4. PLOT: Years Sampled per Month ----
p3 <- ggplot(monthly_summary, 
             aes(x = factor(month), y = n_years)) +
  geom_col(fill = "#A23B72", alpha = 0.8) +
  geom_text(aes(label = n_years), vjust = -0.5, size = 3) +
  scale_x_discrete(labels = month_labels) +
  labs(
    title = "Years of Sampling by Month",
    x = "Month",
    y = "Number of Years"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank()
  )

# 5. MONTHLY SUMMARY BY SPECIES ----

species_month_summary <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(month, species) %>%
  summarise(
    n_sites = n_distinct(site),
    n_years = n_distinct(year),
    detection_days = sum(presence > 0),  # Days with actual detections
    effort_days = n(),  # Total days sampled
    .groups = "drop"
  )

print("\nSpecies by Month Summary:")
print(species_month_summary)

# Check if there's variability in site counts by species
site_species_check_month <- species_month_summary %>%
  select(month, species, n_sites) %>%
  tidyr::pivot_wider(names_from = species, values_from = n_sites)

print("\nSites per month by species (check for variability):")
print(site_species_check_month)

# 6. PLOT: Detection Days by Month and Species ----
p4 <- ggplot(species_month_summary, 
             aes(x = factor(month), y = detection_days, fill = species)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_x_discrete(labels = month_labels) +
  labs(
    title = "Detection Days by Month and Species",
    x = "Month",
    y = "Detection Days",
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 7. PLOT: Proportion of Effort with Detections by Month ----
species_month_prop <- species_month_summary %>%
  mutate(detection_rate = detection_days / effort_days)

p5 <- ggplot(species_month_prop, 
             aes(x = factor(month), y = detection_rate, color = species, group = species)) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2, alpha = 0.8) +
  scale_x_discrete(labels = month_labels) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Detection Rate by Month and Species",
    x = "Month",
    y = "% of Days with Detections",
    color = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave( "output/FIGS/Effort_plots/monthly_detections_SS.png", p5)


# 8. PLOT: Effort Days per Month and Species (to assess representation) ----
p6 <- ggplot(species_month_summary, 
             aes(x = factor(month), y = effort_days, fill = species)) +
  geom_col(position = "dodge", alpha = 0.8) +
  scale_x_discrete(labels = month_labels) +
  geom_hline(yintercept = 100, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(
    title = "Sampling Effort by Month and Species",
    subtitle = "Dashed line = 100 days threshold",
    x = "Month",
    y = "Total Effort Days",
    fill = "Species"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# 9. COMBINED LAYOUT ----
top <- p1 / (p2 | p3)
bottom <- (p4 | p5) & plot_layout(guides = "collect") + 
  theme(legend.position = "bottom")

combined_monthly <- top / bottom
combined_monthly <- combined_monthly + 
  plot_annotation(
    title = "Monthly Summary of Baleen Whale Detections",
    theme = theme(plot.title = element_text(size = 16, face = "bold"))
  )

print(combined_monthly)

# Effort assessment plot
print(p6)

# 10. MONTHLY SITE DISTRIBUTION MAP ----
site_locations_month <- whale_data %>%
  mutate(year = year(UTC)) %>%
  group_by(site, month) %>%
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
bbox <- st_bbox(c(xmin = min(site_locations_month$longitude) - 1, 
                  xmax = max(site_locations_month$longitude) + 1,
                  ymin = min(site_locations_month$latitude) - 1, 
                  ymax = max(site_locations_month$latitude) + 1),
                crs = st_crs(4326))

# Add month labels to the data
site_locations_month <- site_locations_month %>%
  mutate(month_label = factor(month, levels = 1:12, labels = month_labels))

p_map_month <- ggplot() +
  geom_sf(data = land, fill = "gray90", color = "gray50", linewidth = 0.3) +
  geom_point(data = site_locations_month, 
             aes(x = longitude, y = latitude, alpha = n_years),
             color = "#2E86AB", size = 2) +
  scale_alpha_continuous(range = c(0.3, 1), name = "Years Sampled") +
  coord_sf(xlim = c(bbox["xmin"], bbox["xmax"]), 
           ylim = c(bbox["ymin"], bbox["ymax"])) +
  facet_wrap(~month_label, ncol = 4) +
  labs(
    title = "Site-Level Annual Sample Effort by Month",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom") +
  scale_x_continuous(breaks = scales::pretty_breaks(3))

print(p_map_month)

# Save plots
# ggsave("output/FIGS/Effort_plots/combined_monthly.png", combined_monthly, h = 11, w = 11, units = "in")
# ggsave("output/FIGS/Effort_plots/monthly_effort_assessment.png", p6, h = 6, w = 10, units = "in")
# ggsave("output/FIGS/Effort_plots/Site_effort_map_monthly.png", p_map_month, h = 10, w = 12, units = "in")