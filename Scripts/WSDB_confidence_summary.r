# Confidence Summary Script
# This script sources the main WSDB effort mapping script and then calculates
# confidence summaries for OSW areas and the study region by species group
#
# Run after WSDB_effort.R has completed

# ---------------------------
# SETUP
# ---------------------------
message("\n=== Confidence Summary Analysis ===\n")

# Source the main script (assumes it's in the same directory)
# Note: Make sure WSDB_effort.R has completed successfully
if (!exists("cell_all")) {
  stop("cell_all object not found. Please run WSDB_effort.R first or source it.")
}

if (!exists("targets")) {
  stop("targets object not found. Please run WSDB_effort.R first or source it.")
}

# Load required packages
required_pkgs <- c("dplyr", "sf", "tidyr", "readr")
invisible(lapply(required_pkgs, library, character.only = TRUE))

# ---------------------------
# PREPARE GRID CELLS AS SF OBJECT
# ---------------------------
message("Converting grid cells to spatial features...")

# Create sf object from cell_all
cells_sf <- cell_all %>%
  mutate(
    # Create polygon geometry for each cell
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

message(sprintf("  Created %d cell polygons", nrow(cells_sf)))

# ---------------------------
# FUNCTION: CALCULATE CONFIDENCE METRICS
# ---------------------------
calculate_confidence_metrics <- function(cells_data, region_name = "Region") {
  
  # Count cells by confidence category
  conf_counts <- cells_data %>%
    sf::st_drop_geometry() %>%
    count(confidence, .drop = FALSE) %>%
    mutate(
      percent = round(n / sum(n) * 100, 1),
      region = region_name
    )
  
  # Calculate summary statistics
  summary_stats <- cells_data %>%
    sf::st_drop_geometry() %>%
    summarise(
      region = region_name,
      total_cells = n(),
      cells_with_data = sum(n_all > 0),
      total_records = sum(n_all, na.rm = TRUE),
      mean_records_per_cell = round(mean(n_all[n_all > 0], na.rm = TRUE), 1),
      median_records_per_cell = round(median(n_all[n_all > 0], na.rm = TRUE), 1),
      mean_years_per_cell = round(mean(n_years[n_years > 0], na.rm = TRUE), 1),
      median_years_per_cell = round(median(n_years[n_years > 0], na.rm = TRUE), 1),
      pct_high_confidence = round(sum(confidence == "High", na.rm = TRUE) / n() * 100, 1),
      pct_medium_high_confidence = round(sum(confidence %in% c("High", "Medium-High"), na.rm = TRUE) / n() * 100, 1),
      pct_low_confidence = round(sum(confidence %in% c("Low", "Low-Medium"), na.rm = TRUE) / n() * 100, 1)
    )
  
  return(list(
    counts = conf_counts,
    summary = summary_stats
  ))
}

# ---------------------------
# FUNCTION: CALCULATE TARGET-SPECIFIC METRICS
# ---------------------------
calculate_target_metrics <- function(cells_data, target_name, target_regex, dt_data, region_name = "Region") {
  
  message(sprintf("  Processing %s...", target_name))
  
  # Get cell IDs in this region
  region_cell_ids <- cells_data %>%
    sf::st_drop_geometry() %>%
    pull(cell_id)
  
  # Filter dt to cells in this region and calculate target metrics
  if (use_datatable && "data.table" %in% class(dt_data)) {
    target_summary <- dt_data[cell_id %in% region_cell_ids][, .(
      n_all = sum(w, na.rm = TRUE),
      n_target = sum(w * stringr::str_detect(name_blob, target_regex), na.rm = TRUE)
    ), by = .(cell_id)]
    
    target_years <- dt_data[
      cell_id %in% region_cell_ids & stringr::str_detect(name_blob, target_regex),
      .(n_target_years = data.table::uniqueN(year, na.rm = TRUE)),
      by = .(cell_id)
    ]
    
    target_summary <- as_tibble(target_summary)
    target_years <- as_tibble(target_years)
    
  } else {
    target_summary <- dt_data %>%
      filter(cell_id %in% region_cell_ids) %>%
      mutate(is_target = stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(
        n_all = sum(w, na.rm = TRUE),
        n_target = sum(w[is_target], na.rm = TRUE),
        .groups = "drop"
      )
    
    target_years <- dt_data %>%
      filter(cell_id %in% region_cell_ids,
             stringr::str_detect(name_blob, target_regex)) %>%
      group_by(cell_id) %>%
      summarise(n_target_years = n_distinct(year, na.rm = TRUE), .groups = "drop")
  }
  
  # Merge and calculate evidence categories
  target_summary <- target_summary %>%
    left_join(target_years, by = "cell_id") %>%
    mutate(n_target_years = tidyr::replace_na(n_target_years, 0)) %>%
    mutate(
      evidence = case_when(
        n_target == 0 ~ "None",
        n_target >= 5 & n_target_years >= 3 ~ "High",
        n_target >= 2 & n_target_years >= 2 ~ "Medium",
        TRUE ~ "Low"
      ),
      evidence = factor(evidence, levels = c("High", "Medium", "Low", "None"))
    )
  
  # Merge with confidence data
  target_summary <- cells_data %>%
    sf::st_drop_geometry() %>%
    select(cell_id, confidence) %>%
    right_join(target_summary, by = "cell_id") %>%
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
  
  # Count cells by target confidence
  conf_counts <- target_summary %>%
    count(target_conf, .drop = FALSE) %>%
    mutate(
      percent = round(n / sum(n) * 100, 1),
      target = target_name,
      region = region_name
    )
  
  # Summary statistics
  summary_stats <- target_summary %>%
    summarise(
      target = target_name,
      region = region_name,
      total_cells = n(),
      cells_with_target = sum(n_target > 0),
      total_target_records = sum(n_target, na.rm = TRUE),
      mean_target_per_cell = round(mean(n_target[n_target > 0], na.rm = TRUE), 1),
      median_target_per_cell = round(median(n_target[n_target > 0], na.rm = TRUE), 1),
      pct_high_confidence = round(sum(target_conf == "High", na.rm = TRUE) / n() * 100, 1),
      pct_moderate_plus = round(sum(target_conf %in% c("High", "Moderate"), na.rm = TRUE) / n() * 100, 1),
      pct_low_confidence = round(sum(target_conf %in% c("Low", "Low (coverage)"), na.rm = TRUE) / n() * 100, 1),
      pct_no_evidence = round(sum(target_conf == "No evidence", na.rm = TRUE) / n() * 100, 1)
    )
  
  return(list(
    counts = conf_counts,
    summary = summary_stats
  ))
}

# ---------------------------
# STUDY REGION SUMMARY
# ---------------------------
message("\nCalculating study region confidence metrics...")

study_metrics <- calculate_confidence_metrics(cells_sf, "Study Region")

message("  Study region summary:")
message(sprintf("    Total cells: %d", study_metrics$summary$total_cells))
message(sprintf("    Cells with data: %d (%.1f%%)", 
                study_metrics$summary$cells_with_data,
                study_metrics$summary$cells_with_data / study_metrics$summary$total_cells * 100))
message(sprintf("    High confidence: %.1f%%", study_metrics$summary$pct_high_confidence))
message(sprintf("    Medium-High+ confidence: %.1f%%", study_metrics$summary$pct_medium_high_confidence))

# Target-specific metrics for study region
study_target_metrics <- list()
for (nm in names(targets)) {
  study_target_metrics[[nm]] <- calculate_target_metrics(
    cells_sf, nm, target_patterns[[nm]], dt, "Study Region"
  )
}

# ---------------------------
# OSW AREAS SUMMARY
# ---------------------------
if (!is.null(osw_wind)) {
  message("\nCalculating OSW area confidence metrics...")
  
  # Find cells that intersect OSW polygons
  osw_cells <- sf::st_intersection(cells_sf, osw_wind)
  
  if (nrow(osw_cells) > 0) {
    osw_metrics <- calculate_confidence_metrics(osw_cells, "OSW Areas")
    
    message("  OSW areas summary:")
    message(sprintf("    Total cells intersecting OSW: %d", osw_metrics$summary$total_cells))
    message(sprintf("    Cells with data: %d (%.1f%%)", 
                    osw_metrics$summary$cells_with_data,
                    osw_metrics$summary$cells_with_data / osw_metrics$summary$total_cells * 100))
    message(sprintf("    High confidence: %.1f%%", osw_metrics$summary$pct_high_confidence))
    message(sprintf("    Medium-High+ confidence: %.1f%%", osw_metrics$summary$pct_medium_high_confidence))
    
    # Target-specific metrics for OSW areas
    osw_target_metrics <- list()
    for (nm in names(targets)) {
      osw_target_metrics[[nm]] <- calculate_target_metrics(
        osw_cells, nm, target_patterns[[nm]], dt, "OSW Areas"
      )
    }
    
  } else {
    message("  Warning: No cells intersect OSW polygons")
    osw_metrics <- NULL
    osw_target_metrics <- NULL
  }
} else {
  message("\nSkipping OSW summary (no OSW polygons loaded)")
  osw_metrics <- NULL
  osw_target_metrics <- NULL
}

# ---------------------------
# COMPILE RESULTS
# ---------------------------
message("\nCompiling results...")

# Overall confidence summaries
overall_conf_summary <- bind_rows(
  study_metrics$summary,
  if (!is.null(osw_metrics)) osw_metrics$summary else NULL
)

overall_conf_counts <- bind_rows(
  study_metrics$counts,
  if (!is.null(osw_metrics)) osw_metrics$counts else NULL
)

# Target-specific summaries
target_conf_summary <- bind_rows(
  lapply(names(study_target_metrics), function(nm) {
    study_target_metrics[[nm]]$summary
  }),
  if (!is.null(osw_target_metrics)) {
    lapply(names(osw_target_metrics), function(nm) {
      osw_target_metrics[[nm]]$summary
    })
  } else NULL
)

target_conf_counts <- bind_rows(
  lapply(names(study_target_metrics), function(nm) {
    study_target_metrics[[nm]]$counts
  }),
  if (!is.null(osw_target_metrics)) {
    lapply(names(osw_target_metrics), function(nm) {
      osw_target_metrics[[nm]]$counts
    })
  } else NULL
)

# ---------------------------
# SAVE RESULTS
# ---------------------------
message("\nSaving results...")

# Create output directory if it doesn't exist
summary_dir <- file.path(out_dir, "confidence_summaries")
dir.create(summary_dir, showWarnings = FALSE, recursive = TRUE)

# Save CSV files
readr::write_csv(
  overall_conf_summary,
  file.path(summary_dir, "overall_confidence_summary.csv")
)

readr::write_csv(
  overall_conf_counts,
  file.path(summary_dir, "overall_confidence_counts.csv")
)

readr::write_csv(
  target_conf_summary,
  file.path(summary_dir, "target_confidence_summary.csv")
)

readr::write_csv(
  target_conf_counts,
  file.path(summary_dir, "target_confidence_counts.csv")
)

message("  ✓ Saved overall_confidence_summary.csv")
message("  ✓ Saved overall_confidence_counts.csv")
message("  ✓ Saved target_confidence_summary.csv")
message("  ✓ Saved target_confidence_counts.csv")

# ---------------------------
# CREATE SUMMARY VISUALIZATIONS
# ---------------------------
message("\nCreating summary visualizations...")

# Plot 1: Overall confidence comparison (Study Region vs OSW)
if (!is.null(osw_metrics)) {
  p_comparison <- ggplot(overall_conf_counts, 
                         aes(x = confidence, y = percent, fill = region)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = c("Study Region" = "#2c7fb8", "OSW Areas" = "#ff7f00")) +
    labs(
      title = "Data Confidence Comparison: Study Region vs OSW Areas",
      x = "Confidence Category",
      y = "Percentage of Grid Cells (%)",
      fill = "Region"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top"
    )
  
  ggsave(
    file.path(summary_dir, "confidence_comparison.png"),
    p_comparison,
    width = 10,
    height = 6,
    dpi = 300
  )
  message("  ✓ Saved confidence_comparison.png")
}

# Plot 2: Target-specific confidence heatmap (OSW areas only)
# Note: "No Evidence" = zero records for that species/group in the cell
# Even if the cell has good overall coverage, if there are no records of the 
# target species, we can't say anything about that species' presence/absence
if (!is.null(osw_target_metrics)) {
  p_target_heatmap <- target_conf_summary %>%
    filter(region == "OSW Areas") %>%
    select(target, region, pct_high_confidence, pct_moderate_plus, 
           pct_low_confidence, pct_no_evidence) %>%
    pivot_longer(cols = starts_with("pct_"), 
                 names_to = "metric", 
                 values_to = "value") %>%
    mutate(
      metric = case_when(
        metric == "pct_high_confidence" ~ "High",
        metric == "pct_moderate_plus" ~ "Moderate+",
        metric == "pct_low_confidence" ~ "Low",
        metric == "pct_no_evidence" ~ "No Evidence"
      ),
      # Keep the order with No Evidence on the right
      metric = factor(metric, levels = c("High", "Moderate+", "Low", "No Evidence")),
      # Create color category: blue shades for data, grey for no evidence
      color_cat = ifelse(metric == "No Evidence", "no_evidence", "has_data"),
      # For blue scale: higher confidence = darker blue
      # Set to NA if value is 0 so it becomes white
      blue_value = case_when(
        value == 0 ~ NA_real_,
        metric == "High" ~ value * 1.0,
        metric == "Moderate+" ~ value * 0.66,
        metric == "Low" ~ value * 0.33,
        TRUE ~ NA_real_
      ),
      # For grey scale: higher % of no evidence = darker grey
      # But if value is 0, we want it to be white (NA)
      grey_value = ifelse(metric == "No Evidence" & value > 0, value, NA_real_)
    ) %>%
    ggplot(aes(x = metric, y = target)) +
    # Blue tiles for data categories (will be white when blue_value is NA)
    geom_tile(data = . %>% filter(color_cat == "has_data"),
              aes(fill = blue_value), color = "white") +
    # Grey tiles for no evidence (only when > 0%)
    geom_tile(data = . %>% filter(color_cat == "no_evidence" & !is.na(grey_value)),
              aes(alpha = grey_value), fill = "grey60", color = "white") +
    # White tiles for 0% cells (both data categories and no evidence)
    geom_tile(data = . %>% filter(value == 0),
              fill = "white", color = "grey80") +
    geom_text(aes(label = sprintf("%.1f%%", value)), size = 3.5) +
    scale_fill_gradient(
      low = "#deebf7", high = "#08519c",
      name = "% with Data",
      limits = c(0, 100),
      na.value = "white"
    ) +
    scale_alpha_continuous(
      range = c(0.2, 0.8),
      limits = c(0, 100),
      na.value = 0,
      guide = "none"  # Remove legend for no evidence
    ) +
    labs(
      title = "Sightings Data Confidence: OSW Areas",
      subtitle = paste0("Percentage of grid cells in each data confidence category\n",
                        "Blue shades = data available (darker = higher confidence); ",
                        "Grey shades = no records for that species/group"),
      x = "Confidence Category",
      y = "Species Group"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 10),
      plot.subtitle = element_text(size = 9, color = "grey30"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(summary_dir, "target_confidence_heatmap_OSW.png"),
    p_target_heatmap,
    width = 11,
    height = 6,
    dpi = 300
  )
  message("  ✓ Saved target_confidence_heatmap_OSW.png")
} else {
  message("  ⚠ Skipping target heatmap (no OSW data)")
}

# ---------------------------
# PRINT SUMMARY TABLES
# ---------------------------
message("\n=== SUMMARY TABLES ===\n")

message("OVERALL CONFIDENCE SUMMARY:")
print(overall_conf_summary, width = Inf)

message("\n\nTARGET-SPECIFIC CONFIDENCE SUMMARY:")
print(target_conf_summary, width = Inf)

message("\n=== ANALYSIS COMPLETE ===")
message("All results saved to: ", summary_dir)
message("\nFiles created:")
message("  - overall_confidence_summary.csv")
message("  - overall_confidence_counts.csv")
message("  - target_confidence_summary.csv")
message("  - target_confidence_counts.csv")
message("  - confidence_comparison.png")
message("  - target_confidence_heatmap.png")