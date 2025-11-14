# creat KDE maps similar to Feyrer et al. 2024...
######
#PLOT---------
######
# depends on ("Scripts/load_basemap_shapes.R")

species_names_map <- c(
  Bb = "Sei whale",
  Mn = "Humpback whale",
  Bp = "Fin whale",
  Bm = "Blue whale",
  Ba = "Minke whale",
  Eg = "Right Whale"
)

plotKDEMaps <- function(shapefile_dir, 
                        land, contour_data, 
                       baleen_stns,
                        output_dir, predpal = "mako", size = c(11, 8.5)) {
  
  # List all shapefiles in the directory
  shapefiles <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)
  
  # Loop through each shapefile and plot
  for (shapefile in shapefiles) {
    # Read the shapefile
    kde_sf_data <- st_read(shapefile, quiet = T)
    
    # Extract species name from the filename
    species_code <- str_extract(basename(shapefile), "(?<=KDE_pam_)[^\\.]+")
    
    # Choose the appropriate stations data based on species
      current_detects_data <-baleen_stns
    
    # Extract coordinates for labeling
    coords <- st_coordinates(current_detects_data)
    
    # Get proper species name using the mapping
    species_name <- species_names_map[species_code]
    
    # If species name is not found in the map, use the code as a fallback
    if (is.null(species_name)) {
      species_name <- species_code
    }
    
    
    # Extract the bounding box
    bbox <- st_bbox(current_detects_data)
    
    # Extract individual limits
    xmin <- bbox["xmin"]
    xmax <- bbox["xmax"]
    ymin <- bbox["ymin"]
    ymax <- bbox["ymax"]
    
    # Use these limits for xlims and ylims
    xlims <- c(xmin-100000, xmax+100000)
    ylims <- c(ymin- 75000, ymax +75000)
    
    # Get unique quantiles and sort them properly
    unique_quantiles <- unique(kde_sf_data$Quantil)
    cat("Unique quantiles before sorting:", unique_quantiles, "\n")
    
    # Remove NA values if any
    unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
    
    # Extract numeric values for sorting
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    cat("Quantile numbers:", quantile_nums, "\n")
    
    #sort by numeric values
    sort_order <- order(quantile_nums)
    unique_quantiles <- unique_quantiles[sort_order]
    cat("Unique quantiles after sorting:", unique_quantiles, "\n")
    
    # Check if we have valid quantiles
    if (length(unique_quantiles) == 0) {
      cat("Warning: No valid quantiles found. Skipping", basename(shapefile), "\n")
      next
    }
    
    # Convert Quantile to ordered factor
    kde_sf_data$Quantil <- factor(kde_sf_data$Quantil, 
                                   levels = unique_quantiles, 
                                   ordered = TRUE)
    
    # Define palette - match the number of colors to number of quantiles
    n_quantiles <- length(unique_quantiles)
    pal <- rev(get(predpal)(n_quantiles))
    
    # scales::show_col(pal)
    
    # Add truncated site column (first 3 characters)----
    current_detects_data <- current_detects_data %>%
      mutate(site_short = case_when(
        site == "SABW" ~ "SBW",
        site == "ROBP" ~ "RBP",
        site == "ROBV" ~ "RBV",
        site == "GDSE" ~ "GUL",
        str_detect(site, "^MG") ~ "MGU",
        str_detect(site, "^GL") ~ "GUL",
        TRUE ~ substr(site, 1, 3)   # Default: truncate to 3 chars
      ))
    
    # Create unique stations (one per truncated site name)
    unique_stations <- current_detects_data %>%
      group_by(site_short) %>%
      slice(1) %>%
      ungroup()
    
    # Extract coordinates for the unique stations
    unique_coords <- st_coordinates(unique_stations)
    
    
    # Plot gg------
    gg_map <- ggplot() +
      theme_bw() +
      geom_sf(data = contour_data %>% 
                dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), 
              col = "grey50", linewidth = 0.2) +
      geom_sf(data = kde_sf_data, aes(fill = Quantil), col = NA, alpha = .9, na.rm = TRUE) +
      geom_sf(data = current_detects_data, col = "black", shape = 24, fill = "yellow", 
              size = 2, alpha = .5) +
      # Add station labels
      geom_text(data = unique_stations, 
                aes(x = unique_coords[,1], y = unique_coords[,2] + 5000, label = site_short), 
                color = "white", size = 3, vjust = 0, check_overlap = T) +
      geom_sf(data = land, color = NA, fill = "grey50") +
      coord_sf(lims_method = "orthogonal", xlim = xlims, ylim = ylims, 
               crs = UTM20, expand = TRUE) +
      # Add species name as title 
      labs(title = paste("Acoustic Presence Map for", species_name)) +
      ylab("") + 
      xlab("") +
      theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(0.1, 0.9), 
            legend.background = element_rect(fill = NA), 
            legend.key = element_rect(fill = NA), 
            plot.margin = margin(0, 0, 0, 0), 
            plot.title = element_text(hjust = 0.5), 
            plot.subtitle = element_text(hjust = 0.5), 
            axis.title = element_blank()) +
      scale_fill_manual(values = pal, 
                        na.value = NA, 
                        na.translate = FALSE,
                        name = "Quantile",
                        drop = FALSE) +
      annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6, 
                       bar_cols = c("grey40", "white"))
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Output file name
    output_file <- file.path(output_dir, paste0(gsub("\\.shp$", "", basename(shapefile)), ".png"))
    
    # Save the plot
    ggsave(output_file, gg_map, width = size[1], height = size[2], units = "in")
    
    cat("Saved:", output_file, "\n")
  }
}

# Run the function
plotKDEMaps(shapefile_dir = "output/shapes/", 
            land = landUTM, 
            contour_data = cont_UTM, 
            baleen_stns = baleen_stns,
            output_dir = "output/FIGS/")