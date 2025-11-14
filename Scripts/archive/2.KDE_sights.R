#
#Script to read in sightings data on whales from wsdb and export KDE surfaces of presence
#Updated 2025

# Load necessary libraries
options(scipen = 999)

    #library
    pacman::p_load(sf, tidyverse, readxl, here, leaflet, 
                   scales, terra, ggrepel, viridis, ggspatial, spatstat, dplyr,patchwork)
    
# depends on load_basemap_shapes.R
    
    raster_dir = "output/tif/"
    plot_dir = "output/FIGS/"
    
#projections------
UTM20 <- crs("+init=epsg:32620") # CODE FOR UTM Zone 20

    
# Read sightings records---------
    # Set the path to folder containing CSV files
    folder_path <- "input/2025/sights/"
    
    # Get a list of CSV files in the folder
    csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    
    # Read each file and combine them
    df_combined <- csv_files %>% 
      map_dfr(read_csv, show_col_types = FALSE, col_types = cols(.default ="c"))
    
    whale_data <-df_combined%>%mutate(UTC = as.POSIXct(WS_DATETIME, format = "%Y-%m-%d %H:%M.%S", tz="UTC"))%>%
      mutate(Date = format(as_date(UTC), "%Y-%m-%d"))%>%mutate(Month = format(as_date(UTC), "%m"))
    
    # Date range
    whale_data%>%summarise(min(Date), max(Date))
    
    #add species codes & Sci names
    species = "speciesCodes.csv" 
    
    # check species, add common names and scientific names based on codes-----
    SP_data <- read_csv(here("input", species), col_types = cols(.default ="c"))
    WS_data = left_join(whale_data, SP_data)
    
    # bound boxes for study area -----
    bbox = Bound_boxB

    # Create point shapefile from records
    WS_data_sf = st_as_sf(WS_data, coords= c("LONGITUDE","LATITUDE"), crs = sf::st_crs(4326) )%>%
      st_intersection(bbox)%>%st_transform(UTM20)
    
    # Ensure same CRS
    WS_data_sf <- st_make_valid(WS_data_sf)
    landUTM    <- st_make_valid(landUTM)
    WS_data_sf <- st_transform(WS_data_sf, st_crs(landUTM))
    
   
    
    # Define vector of cetacean species codes (excluding NS categories and seals)
    cetacean_codes <- c(
      921, 922, 923, 925, 931, 932, 933, 934, 935, 936, 937, 938,
      7019, 7020, 7021, 7022, 7023, 7024, 7025, 7026, 7027, 7028,
      7029, 7031, 7033, 7034, 7035, 7037, 7038, 902
    )
    
    # 3) Filter: cetaceans -> not on land -> species with >50 records
    WS_data_sf <- WS_data_sf %>%
      filter(SPECIES_CD %in% cetacean_codes) %>%
      # keep points that DO NOT intersect land
      filter(lengths(st_intersects(., landUTM)) == 0) %>%
      add_count(SPECIES_CD) %>%
      filter(n > 50) %>%
      select(-n)
    
    # Extract unique species names
    unique_species <- unique(WS_data_sf$COMMONNAME)
    
    # Determine the number of unique species
    num_species <- length(unique_species)
    
    # Use all PAM stations to define a common window-----
    pam_coords <- st_coordinates(baleen_sf)
    
    x_range <- range(pam_coords[, "X"])
    y_range <- range(pam_coords[, "Y"])
    
    buffer_percent <- 0.2  # same as you use in KDE
    
    x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
    y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
    
    expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
    expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
    
    pam_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)

performKDE <- function(data_sf, unique_species, buffer_percent, window, sigma_func, output_dir = raster_dir) {
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  shapefile_crs <- st_crs(data_sf)$wkt
  
  # Initialize a list to store plots
  raster_list <- list()
  sigma_store <- list() 
  
  global_min <- Inf
  global_max <- -Inf
  
  # Add sanitize function
  sanitize_filename <- function(name) {
    name <- trimws(name)
    name <- gsub("[/\\:*?\"<>|]", "_", name)
    name <- gsub(" ", "", name)
    return(name)
  }
  
  # Loop through each species 1rst time
  for(species in unique_species) {
    current_sf <- data_sf[data_sf$COMMONNAME == species, ]
    
    if(nrow(current_sf) > 0) {
      # Initialize to NULL in case of errors
      raster_kd <- NULL
      sigma_val <- NULL
      
      tryCatch({
      species_coords <- st_coordinates(current_sf)
      
      # Create ppp object within the precomputed window
      points_ppp <- ppp(x = species_coords[,"X"], y = species_coords[,"Y"], window = window)
      
      #check for duplicates
      jitter_bur <- rjitter(points_ppp, retry=TRUE, nsim=1, drop=TRUE)

      # Apply the bandwidth function to the point pattern object
      sigma_val <- 10000
      sigma_store[[species]] <- sigma_val
      
      # Kernel Density Estimation
      kd_result <- density.ppp(jitter_bur, 
                               sigma = sigma_val, 
                               positive = T,
                               kernel = "gaussian", 
                               dimyx = 500, 
                               diggle = T)
      
      
      # Convert to SpatRaster
      raster_kd <- rast(kd_result)
      # Set the CRS for the raster
      crs(raster_kd) <- shapefile_crs
      
      # Normalize the raster values by the maximum value within the grid for the species
      max_val <- max(values(raster_kd), na.rm = TRUE)
      if (!is.na(max_val) && max_val > 0) {
        raster_kd <- raster_kd / max_val
        # Threshold: remove %tails of values-----
        # This keeps only the top 5% most probable areas to remove noise
        kde_vals <- values(raster_kd, mat = FALSE, na.rm = TRUE)
        threshold <- quantile(kde_vals, 0.95, na.rm = TRUE)  # 95th percentile
        
        cat("Species:", species, "- Using threshold:", threshold, 
            "(removes", sum(kde_vals < threshold)/length(kde_vals)*100, "% of pixels)\n")
        
        raster_kd[raster_kd < threshold] <- NA
      }

 
      
      # Extract min/max from the raster
      local_min <- min(values(raster_kd), na.rm = TRUE)
      local_max <- max(values(raster_kd), na.rm = TRUE)
      
      # Update global min/max
      global_min <- min(global_min, local_min)
      global_max <- max(global_max, local_max)
      
      
      # Save the raster to a file - WITH SANITIZED FILENAME
      raster_filename <- paste0(output_dir,  "KDE_sights_", sanitize_filename(species),".tif")
      writeRaster(raster_kd, filename = raster_filename, overwrite = T)
      
      # Store the raster AND sigma in a list - THIS WAS MISSING
      raster_list[[species]] <- list(raster = raster_kd, sigma = sigma_val)
      
      }, error = function(e) {
        cat("Error processing species", species, ":", e$message, "\n")
      })
    }
  }


  # Second pass: Create plots using stored rasters and global min/max
  plot_list <- list()

  for(species in unique_species) {
    raster_data <- raster_list[[species]]  # Get the list containing both raster and sigma
    if(!is.null(raster_data)) {
      raster_kd <- raster_data$raster      # Extract the raster
      sigma_val <- raster_data$sigma       # Extract the sigma value
      
      # Convert raster to data frame for ggplot
      kde_df <- as.data.frame(raster_kd, xy = TRUE) %>%
        filter(!is.na(lyr.1))

      # Create ggplot using global range for colour scale
      p <- ggplot(kde_df, aes(x = x, y = y, fill = lyr.1)) +
        geom_tile() +
        scale_fill_viridis_c(limits = c(global_min, global_max), label = scales::scientific) +
        theme_minimal() +
        theme(axis.title= element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(), axis.line = element_blank()) +
        labs(fill = "", title =  species) +
        annotate("text", x = Inf, y = Inf, label = paste("Bandwidth:", round(sigma_val, 0)),
                 hjust = 1.5, vjust = 2, size = 3, color = "orange")+
        coord_equal()


      # Add the plot to the list
      plot_list[[species]] <- p

    } else {
      cat("No data found for", species, "\n")
    }
  }
  # Combine all plots using patchwork
  plot_combined <- wrap_plots(plot_list, guides = "collect")

  # Save the combined plot
  combined_filename <- paste0(plot_dir, "combined_KDE_plots.png")
  ggsave(filename = combined_filename, 
         plot = plot_combined, 
         width = 16, height = 12, dpi = 300)
  
  # Also save individual plots
  for(species in names(plot_list)) {
    if(!is.null(plot_list[[species]])) {
      individual_filename <- paste0(plot_dir,  "KDE_", sanitize_filename(species),".png")
      ggsave(filename = individual_filename, 
             plot = plot_list[[species]], 
             width = 8, height = 6, dpi = 300)
    }
  }
  
  list(plot = plot_combined,
       bandwidths = sigma_store)
}


#run KDE function-----
res = performKDE(data_sf = WS_data_sf, unique_species = unique_species, window = pam_window,
           sigma_func = bw.diggle)

# Plot
res$plot

# Bandwidths per species
res$bandwidths


