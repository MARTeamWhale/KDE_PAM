# Seasonal Analysis of Cetacean Sightings
# Creates faceted maps showing all cetacean records by season
# Run after WSDB_effort.R has completed successfully

# ---------------------------
# SETUP
# ---------------------------
message("\n=== Seasonal Coverage Analysis ===\n")

# Check that main script has been run
if (!exists("dt")) {
  stop("dt object not found. Please run WSDB_effort.R first or source it.")
}

if (!exists("cell_m")) {
  stop("cell_m object not found. Please run WSDB_effort.R first or source it.")
}

# Load required packages
required_pkgs <- c("dplyr", "ggplot2", "sf", "lubridate")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---------------------------
# ADD SEASONAL CLASSIFICATION
# ---------------------------
message("Classifying records by season...")

# Define seasons (Northern Hemisphere)
# Winter: Dec, Jan, Feb
# Spring: Mar, Apr, May
# Summer: Jun, Jul, Aug
# Fall: Sep, Oct, Nov

if (use_datatable && "data.table" %in% class(dt)) {
  dt[, month := lubridate::month(date_utc)]
  dt[, season := data.table::fcase(
    month %in% c(12, 1, 2), "Winter",
    month %in% c(3, 4, 5), "Spring",
    month %in% c(6, 7, 8), "Summer",
    month %in% c(9, 10, 11), "Fall",
    default = NA_character_
  )]
  dt[, season := factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))]
  
  # Get summary counts
  seasonal_summary <- dt[, .(n_records = .N), by = season]
  
} else {
  dt <- dt %>%
    mutate(
      month = lubridate::month(date_utc),
      season = case_when(
        month %in% c(12, 1, 2) ~ "Winter",
        month %in% c(3, 4, 5) ~ "Spring",
        month %in% c(6, 7, 8) ~ "Summer",
        month %in% c(9, 10, 11) ~ "Fall",
        TRUE ~ NA_character_
      ),
      season = factor(season, levels = c("Spring", "Summer", "Fall", "Winter"))
    )
  
  # Get summary counts
  seasonal_summary <- dt %>%
    group_by(season) %>%
    summarise(n_records = n(), .groups = "drop")
}

message("\nSeasonal distribution of records:")
print(seasonal_summary)

# Calculate percentages
if (use_datatable && "data.table" %in% class(seasonal_summary)) {
  total_records <- sum(seasonal_summary$n_records, na.rm = TRUE)
  seasonal_summary[, pct := round(n_records / total_records * 100, 1)]
} else {
  total_records <- sum(seasonal_summary$n_records, na.rm = TRUE)
  seasonal_summary <- seasonal_summary %>%
    mutate(pct = round(n_records / total_records * 100, 1))
}

message("\nPercentage by season:")
print(seasonal_summary)

# ---------------------------
# CALCULATE SEASONAL GRID SUMMARIES
# ---------------------------
message("\nCalculating seasonal grid summaries...")

if (use_datatable && "data.table" %in% class(dt)) {
  cell_seasonal <- dt[, .(
    n_all = sum(w, na.rm = TRUE)
  ), by = .(cell_id, gx, gy, x0, y0, season)]
  
  cell_seasonal <- as_tibble(cell_seasonal)
} else {
  cell_seasonal <- dt %>%
    group_by(cell_id, gx, gy, x0, y0, season) %>%
    summarise(n_all = sum(w, na.rm = TRUE), .groups = "drop")
}

# Remove NA seasons
cell_seasonal <- cell_seasonal %>%
  filter(!is.na(season))

message(sprintf("  Created seasonal summaries for %d unique cell-season combinations", 
                nrow(cell_seasonal)))

# ---------------------------
# CREATE SEASONAL FACETED MAP
# ---------------------------
message("\nCreating seasonal faceted map...")

# Apply log transform to match main coverage map
cell_seasonal <- cell_seasonal %>%
  mutate(n_all_display = log10(n_all + 1))

# Create breaks at powers of 10 for intuitive reading
log_breaks <- seq(0, ceiling(max(cell_seasonal$n_all_display, na.rm = TRUE)), by = 1)
actual_values <- 10^log_breaks
actual_values[1] <- 0  # First break is 0, not 1

# Create the plot
p_seasonal <- ggplot(cell_seasonal) +
  geom_tile(aes(x = x0, y = y0, fill = n_all_display)) +
  scale_fill_viridis_c(
    option = "viridis", 
    na.value = "grey90", 
    name = "Count",
    breaks = log_breaks,
    labels = actual_values,
    begin = 0.15, 
    end = 0.95
  ) +
  facet_wrap(~season, ncol = 2) +
  labs(
    title = "Seasonal Distribution of All Cetacean Sightings",
    subtitle = paste0(year_min, "–", year_max, " | ", grid_km, " km grid | log scale"),
    caption = paste0("Seasons: Spring (Mar-May), Summer (Jun-Aug), Fall (Sep-Nov), Winter (Dec-Feb)\n",
                     "Record counts: Spring=", seasonal_summary$n_records[seasonal_summary$season=="Spring"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Spring"], "%), ",
                     "Summer=", seasonal_summary$n_records[seasonal_summary$season=="Summer"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Summer"], "%), ",
                     "Fall=", seasonal_summary$n_records[seasonal_summary$season=="Fall"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Fall"], "%), ",
                     "Winter=", seasonal_summary$n_records[seasonal_summary$season=="Winter"],
                     " (", seasonal_summary$pct[seasonal_summary$season=="Winter"], "%)")
  ) +
  coord_sf(xlim = xlims, ylim = ylims, crs = 32620, expand = FALSE, clip = "on") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 10),
    plot.caption = element_text(size = 8, hjust = 0),
    strip.text = element_text(size = 11, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "black")
  )

# Add spatial overlays to each facet
if (!is.null(land)) {
  p_seasonal <- p_seasonal + 
    geom_sf(data = land, fill = "grey60", color = NA, inherit.aes = FALSE)
}

if (!is.null(osw_wind)) {
  p_seasonal <- p_seasonal + 
    geom_sf(data = osw_wind, fill = NA, color = owa_color, 
            alpha = owa_alpha, linewidth = owa_linewidth, inherit.aes = FALSE)
}

if (!is.null(study_area)) {
  p_seasonal <- p_seasonal + 
    geom_sf(data = study_area, fill = NA, color = "black", 
            linewidth = 0.6, linetype = "dashed", inherit.aes = FALSE)
}

# Save the plot
ggsave(
  file.path(out_dir, "00_seasonal_coverage.png"),
  p_seasonal,
  width = 12,
  height = 10,
  dpi = map_dpi
)

message("  ✓ Saved seasonal coverage map")

# ---------------------------
# CALCULATE SEASONAL STATISTICS FOR OSW AREAS
# ---------------------------
if (!is.null(osw_wind)) {
  message("\nCalculating seasonal statistics for OSW areas...")
  
  # Create spatial object from seasonal grid
  cells_seasonal_sf <- cell_seasonal %>%
    mutate(
      geometry = purrr::map2(x0, y0, function(x, y) {
        sf::st_polygon(list(matrix(
          c(x, y,
            x + cell_m, y,
            x + cell_m, y + cell_m,
            x, y + cell_m,
            x, y),
          ncol = 2, byrow = TRUE
        )))
      })
    ) %>%
    sf::st_as_sf(crs = 32620)
  
  # Find cells that intersect OSW polygons
  osw_cells_seasonal <- sf::st_intersection(cells_seasonal_sf, osw_wind)
  
  if (nrow(osw_cells_seasonal) > 0) {
    # Summary by season
    osw_seasonal_summary <- osw_cells_seasonal %>%
      sf::st_drop_geometry() %>%
      group_by(season) %>%
      summarise(
        n_cells_with_data = sum(n_all > 0),
        total_records = sum(n_all, na.rm = TRUE),
        mean_records = round(mean(n_all[n_all > 0], na.rm = TRUE), 1),
        median_records = round(median(n_all[n_all > 0], na.rm = TRUE), 1),
        .groups = "drop"
      )
    
    message("\nOSW Areas - Seasonal Summary:")
    print(osw_seasonal_summary, width = Inf)
    
    # Save to CSV
    readr::write_csv(
      osw_seasonal_summary,
      file.path(out_dir, "confidence_summaries", "osw_seasonal_summary.csv")
    )
    message("  ✓ Saved osw_seasonal_summary.csv")
    
  } else {
    message("  Warning: No cells intersect OSW polygons")
  }
}

# ---------------------------
# CREATE SEASONAL COMPARISON BAR CHART
# ---------------------------
message("\nCreating seasonal comparison chart...")

# Calculate cells with data by season
seasonal_cells <- cell_seasonal %>%
  group_by(season) %>%
  summarise(
    cells_with_data = sum(n_all > 0),
    total_records = sum(n_all),
    .groups = "drop"
  )

# Create bar chart
p_seasonal_bars <- ggplot(seasonal_cells, aes(x = season, y = total_records, fill = season)) +
  geom_col() +
  geom_text(aes(label = paste0(format(total_records, big.mark = ","), " records\n",
                               cells_with_data, " cells")),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c(
    "Spring" = "#78c679",
    "Summer" = "#fd8d3c", 
    "Fall" = "#d95f0e",
    "Winter" = "#6baed6"
  )) +
  labs(
    title = "Seasonal Distribution of Cetacean Sightings",
    subtitle = paste0(year_min, "–", year_max),
    x = "Season",
    y = "Total Records",
    caption = "Number of grid cells with data shown above each bar"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.15)))

ggsave(
  file.path(out_dir, "00_seasonal_comparison.png"),
  p_seasonal_bars,
  width = 8,
  height = 6,
  dpi = 300
)

message("  ✓ Saved seasonal comparison chart")

# ---------------------------
# SUMMARY STATISTICS
# ---------------------------
message("\n=== SEASONAL ANALYSIS COMPLETE ===")
message("\nOverall seasonal distribution:")
print(seasonal_summary)

message("\nKey findings:")
message(sprintf("  - Most sampled season: %s (%.1f%%)", 
                seasonal_summary$season[which.max(seasonal_summary$pct)],
                max(seasonal_summary$pct, na.rm = TRUE)))
message(sprintf("  - Least sampled season: %s (%.1f%%)", 
                seasonal_summary$season[which.min(seasonal_summary$pct)],
                min(seasonal_summary$pct, na.rm = TRUE)))

# Calculate ratio
max_season_records <- max(seasonal_summary$n_records, na.rm = TRUE)
min_season_records <- min(seasonal_summary$n_records, na.rm = TRUE)
ratio <- round(max_season_records / min_season_records, 1)

message(sprintf("  - Sampling ratio (max/min): %.1fx more records in peak season", ratio))

message("\nFiles created:")
message("  - 00_seasonal_coverage.png (faceted map)")
message("  - 00_seasonal_comparison.png (bar chart)")
if (!is.null(osw_wind)) {
  message("  - osw_seasonal_summary.csv")
}