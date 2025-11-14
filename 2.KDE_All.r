#
# Unified script for KDE analysis of whale sightings and PAM data
# Updated 2025
#
# ...............
# DATASET SELECTION - CONFIGURE HERE------
# ...............

RUN_SIGHTINGS <- TRUE   # Set to TRUE to process sightings data
RUN_PAM <- TRUE         # Set to TRUE to process PAM data

# ..................

# Load necessary libraries
options(scipen = 999)

pacman::p_load(sf, tidyverse, readxl, here, leaflet, 
               scales, terra, ggrepel, viridis, ggspatial, spatstat, 
               dplyr, patchwork)

# depends on load_basemap_shapes.R

# Directory setup
raster_dir <- "output/tif/"
plot_dir <- "output/FIGS/"

# Projections
UTM20 <- crs("+init=epsg:32620") # CODE FOR UTM Zone 20

# ................................
# UNIFIED KDE FUNCTION------
# ................................

performKDE <- function(data_sf, 
                       species_col,
                       species_list, 
                       weight_col = NULL,
                       buffer_percent = 0.2, 
                       sigma_val = 10000,
                       window = NULL,
                       output_prefix,
                       threshold_quantile = NULL,
                       output_dir = "output/tif/",
                       plot_output_dir = "output/FIGS/") {
  
  # Ensure output directories exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
  
  shapefile_crs <- st_crs(data_sf)$wkt
  
  # Initialize storage
  raster_list <- list()
  sigma_store <- list() 
  global_min <- Inf
  global_max <- -Inf
  
  # Sanitize filename helper
  sanitize_filename <- function(name) {
    name <- trimws(name)
    name <- gsub("[/\\:*?\"<>|]", "_", name)
    name <- gsub(" ", "", name)
    return(name)
  }
  
  #.............................
  # FIRST PASS: Create KDE rasters-----
  #.............................
  
  
  for(species in species_list) {
    current_sf <- data_sf[data_sf[[species_col]] == species, ]
    
    if(nrow(current_sf) > 0) {
      raster_kd <- NULL
      sigma_used <- NULL
      
      tryCatch({
        species_coords <- st_coordinates(current_sf)
        
        # Create or use provided window
        if(is.null(window)) {
          x_range <- range(species_coords[, "X"])
          y_range <- range(species_coords[, "Y"])
          x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
          y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
          expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
          expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
          species_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
        } else {
          species_window <- window
        }
        
        # Create ppp object
        points_ppp <- ppp(x = species_coords[,"X"], 
                          y = species_coords[,"Y"], 
                          window = species_window)
        
        # Handle duplicates with jitter
        points_ppp <- rjitter(points_ppp, retry = TRUE, nsim = 1, drop = TRUE)
        
        # Apply weights if provided
        if(!is.null(weight_col)) {
          weights <- current_sf[[weight_col]] * 10000
          marks(points_ppp) <- data.frame(weights)
          weight_arg <- marks(points_ppp)
        } else {
          weight_arg <- NULL
        }
        
        sigma_used <- sigma_val
        sigma_store[[species]] <- sigma_used
        
        # Kernel Density Estimation
        kd_result <- density.ppp(points_ppp, 
                                 sigma = sigma_used, 
                                 positive = TRUE,
                                 kernel = "gaussian", 
                                 weights = weight_arg,
                                 dimyx = 500, 
                                 diggle = TRUE)
        
        # Convert to SpatRaster
        raster_kd <- rast(kd_result)
        crs(raster_kd) <- shapefile_crs
        
        # Normalize by species max
        max_val <- max(values(raster_kd), na.rm = TRUE)
        if (!is.na(max_val) && max_val > 0) {
          raster_kd <- raster_kd / max_val
          
          # Apply threshold if specified (for sightings)
          if(!is.null(threshold_quantile)) {
            kde_vals <- values(raster_kd, mat = FALSE, na.rm = TRUE)
            threshold <- quantile(kde_vals, threshold_quantile, na.rm = TRUE)
            
            cat("Species:", species, "- Using threshold:", threshold, 
                "(removes", sum(kde_vals < threshold)/length(kde_vals)*100, "% of pixels)\n")
            
            raster_kd[raster_kd < threshold] <- NA
          }
        }
        
        # Update global min/max
        local_min <- min(values(raster_kd), na.rm = TRUE)
        local_max <- max(values(raster_kd), na.rm = TRUE)
        global_min <- min(global_min, local_min)
        global_max <- max(global_max, local_max)
        
        # Save raster
        raster_filename <- paste0(output_dir, "KDE_", output_prefix, "_", 
                                  sanitize_filename(species), ".tif")
        writeRaster(raster_kd, filename = raster_filename, overwrite = TRUE)
        
        # Store raster and sigma
        raster_list[[species]] <- list(raster = raster_kd, sigma = sigma_used)
        
      }, error = function(e) {
        cat("Error processing species", species, ":", e$message, "\n")
      })
    }
  }
  
  #.............................
    # SECOND PASS: Create plots with consistent scale------
  #.............................
  
  
  plot_list <- list()
  
  for(species in species_list) {
    raster_data <- raster_list[[species]]
    if(!is.null(raster_data)) {
      raster_kd <- raster_data$raster
      sigma_used <- raster_data$sigma
      
      # Convert raster to data frame for ggplot
      kde_df <- as.data.frame(raster_kd, xy = TRUE) %>%
        filter(!is.na(lyr.1))
      
      # Create ggplot using global range for colour scale
      p <- ggplot(kde_df, aes(x = x, y = y, fill = lyr.1)) +
        geom_tile() +
        scale_fill_viridis_c(limits = c(global_min, global_max), 
                             label = scales::scientific) +
        theme_minimal() +
        theme(axis.title = element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(), 
              axis.line = element_blank()) +
        labs(fill = "", title = species) +
        annotate("text", x = Inf, y = Inf, 
                 label = paste("Bandwidth:", round(sigma_used, 0)),
                 hjust = 1.5, vjust = 2, size = 3, color = "orange") +
        coord_equal()
      
      plot_list[[species]] <- p
    } else {
      cat("No data found for", species, "\n")
    }
  }
  
  # Combine all plots using patchwork
  plot_combined <- wrap_plots(plot_list, guides = "collect")
  
  # Save combined plot
  combined_filename <- paste0(plot_output_dir, "combined_KDE_", output_prefix, "_plots.png")
  ggsave(filename = combined_filename, 
         plot = plot_combined, 
         width = 16, height = 12, dpi = 300)
  
  # Save individual plots
  for(species in names(plot_list)) {
    if(!is.null(plot_list[[species]])) {
      individual_filename <- paste0(plot_output_dir, "KDE_", output_prefix, "_", 
                                    sanitize_filename(species), ".png")
      ggsave(filename = individual_filename, 
             plot = plot_list[[species]], 
             width = 8, height = 6, dpi = 300)
    }
  }
  
  return(list(
    plot = plot_combined,
    bandwidths = sigma_store,
    global_range = c(global_min, global_max)
  ))
}


#.............................
    # CREATE COMMON WINDOW FROM PAM STATIONS-------
#.............................

# Define common window based on PAM stations
# This ensures both datasets use the same spatial extent for comparison
common_window <- NULL
baleen_sf <- NULL  # Initialize

if(RUN_PAM || RUN_SIGHTINGS) {
  cat("\n=== CREATING COMMON WINDOW FROM PAM STATIONS ===\n\n")
  
  # Read baleen whale PAM data
  filepath <- 'input/2025/baleen_presence_days_laura_2025.csv'
  baleen_PA <- read.csv(filepath)
  baleen_sf <- baleen_PA %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(crs = UTM20)
  
  # Create window from PAM station extents
  pam_coords <- st_coordinates(baleen_sf)
  x_range <- range(pam_coords[, "X"])
  y_range <- range(pam_coords[, "Y"])
  
  buffer_percent <- 0.2
  x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
  y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
  
  expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
  expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
  
  common_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
  
  cat("Common window created from PAM station extents\n")
  cat("X range:", expanded_x_range, "\n")
  cat("Y range:", expanded_y_range, "\n")
}


#.............................
# PROCESS SIGHTINGS DATA----
#.............................


if(RUN_SIGHTINGS) {
  
  cat("\n=== PROCESSING SIGHTINGS DATA ===\n\n")
  
  # Read sightings records
  folder_path <- "input/2025/sights/"
  csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
  
  df_combined <- csv_files %>% 
    map_dfr(read_csv, show_col_types = FALSE, col_types = cols(.default = "c"))
  
  whale_data <- df_combined %>%
    mutate(UTC = as.POSIXct(WS_DATETIME, format = "%Y-%m-%d %H:%M.%S", tz = "UTC")) %>%
    mutate(Date = format(as_date(UTC), "%Y-%m-%d")) %>%
    mutate(Month = format(as_date(UTC), "%m"))
  
  # Date range
  whale_data %>% summarise(min(Date), max(Date))
  
  # Add species codes & Sci names
  species <- "speciesCodes.csv" 
  SP_data <- read_csv(here("input", species), col_types = cols(.default = "c"))
  WS_data <- left_join(whale_data, SP_data)
  
  # Create point shapefile from records
  bbox <- Bound_boxB
  
  WS_data_sf <- st_as_sf(WS_data, coords = c("LONGITUDE","LATITUDE"), 
                         crs = sf::st_crs(4326)) %>%
    st_intersection(bbox) %>%
    st_transform(UTM20)
  
  # Ensure same CRS
  WS_data_sf <- st_make_valid(WS_data_sf)
  landUTM <- st_make_valid(landUTM)
  WS_data_sf <- st_transform(WS_data_sf, st_crs(landUTM))
  
  # Define vector of cetacean species codes
  cetacean_codes <- c(
    921, 922, 923, 925, 931, 932, 933, 934, 935, 936, 937, 938,
    7019, 7020, 7021, 7022, 7023, 7024, 7025, 7026, 7027, 7028,
    7029, 7031, 7033, 7034, 7035, 7037, 7038, 902
  )
  
  # Filter: cetaceans -> not on land -> species with >50 records
  WS_data_sf <- WS_data_sf %>%
    filter(SPECIES_CD %in% cetacean_codes) %>%
    filter(lengths(st_intersects(., landUTM)) == 0) %>%
    add_count(SPECIES_CD) %>%
    filter(n > 50) %>%
    select(-n)
  
  # Extract unique species names
  unique_species_sightings <- unique(WS_data_sf$COMMONNAME)
  
  cat("Processing", length(unique_species_sightings), "species from sightings data\n")
  
  # Run KDE for sightings using common window
  res_sightings <- performKDE(
    data_sf = WS_data_sf, 
    species_col = "COMMONNAME",
    species_list = unique_species_sightings, 
    weight_col = NULL,
    window = common_window,
    output_prefix = "sights",
    threshold_quantile = NULL,
    sigma_val = 10000
  )
  
  # Display results
  print(res_sightings$plot)
  # cat("\nSightings bandwidths:\n")
  # print(res_sightings$bandwidths)
}


#.............................
# PROCESS PAM DATA---------
#.............................


if(RUN_PAM) {
  
  cat("\n=== PROCESSING PAM DATA ===\n\n")
  
  # Baleen species
  # Bb = Sei, Bm = Blue, Bp = Fin, Mn = humpback, Ba = Minke, Eg = NARW
  baleen_PA %>% group_by(species) %>% summarise(count = n())
  
  species_list_pam <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")
  
  cat("Processing", length(species_list_pam), "species from PAM data\n")
  
  # Run KDE for PAM using common window
  res_pam <- performKDE(
    data_sf = baleen_sf, 
    species_col = "species",
    species_list = species_list_pam, 
    weight_col = "proportion_det",
    window = common_window,
    output_prefix = "pam",
    threshold_quantile = 0.95,  
    sigma_val = 10000
  )
  
  # Display results
  print(res_pam$plot)
  # cat("\nPAM bandwidths:\n")
  # print(res_pam$bandwidths)
}


#.............................
# SUMMARY------
#.............................

cat("\n=== ANALYSIS COMPLETE ===\n")
if(RUN_SIGHTINGS) cat("✓ Sightings data processed\n")
if(RUN_PAM) cat("✓ PAM data processed\n")
cat("\nOutputs saved to:\n")
cat("  Rasters:", raster_dir, "\n")
cat("  Plots:", plot_dir, "\n")