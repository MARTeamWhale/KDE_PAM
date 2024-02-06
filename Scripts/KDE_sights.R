#
#Script to read in sightings data on whales from wsdb and export KDE surfaces of presence
#February, 2024

# Load necessary libraries
# install.packages("spatstat")
options(scipen = 999)

    #library
    pacman::p_load(sf, tidyverse, readxl, here, leaflet, scales, terra, ggrepel, viridis, ggspatial, spatstat, dplyr,patchwork)
    

#projections------
UTM20 <- crs("+init=epsg:32620") # CODE FOR UTM Zone 20

    
# Read sightings records---------
    # Set the path to folder containing CSV files
    folder_path <- "input/sights/"
    
    # Get a list of CSV files in the folder
    csv_files <- list.files(folder_path, pattern = "\\.csv$", full.names = TRUE)
    
    # Read each file and combine them
    df_combined <- csv_files %>% 
      map_dfr(read_csv, show_col_types = FALSE, col_types = cols(.default ="c"))
    
    whale_data <-df_combined%>%mutate(UTC = as.POSIXct(WS_DATETIME, format = "%Y-%m-%d %H:%M:%S", tz="UTC"))%>%
      mutate(Date = format(as_date(UTC), "%Y-%m-%d"))%>%mutate(Month = format(as_date(UTC), "%m"))
    
    #filter out pre-GPS years
    Pre_GPS = whale_data%>%filter(Date <= "1980-01-01")
    whale_data =  whale_data%>%filter(!Date %in% Pre_GPS$Date)
    
    #add species codes & Sci names
    species = "speciesCodes.csv" 
    
    # check species, add common names and scientific names based on codes-----
    SP_data <- read_csv(here("input", species), col_types = cols(.default ="c"))
    WS_data = left_join(whale_data, SP_data)
    
    # bound boxes for study area -----
    bbox = st_bbox( c(xmin = -67.5,ymin = 40, xmax = -55, ymax =47.5 ), crs = st_crs(4326))%>% st_as_sfc()%>% #turns the bounding box into a sfc object, that just describes a specific geometry
      st_sf()
    # # Extract the bounding box to be the same as PAM stations
    # # bbox <- st_bbox(baleen_stns)
    # # Expand the bounding box
    # # Convert 10 km to meters (as the CRS is in UTM)
    # distance_m <- 10000 
    # expanded_bbox <- bbox + c(-distance_m, -distance_m, distance_m, distance_m)
    # bbox =  expanded_bbox%>% st_as_sfc()%>% #turns the bounding box into a sfc object, that just describes a specific geometry
    #     st_sf()%>%st_transform(4326)
    # ext(bbox)

    # Create point shapefile from records
    WS_data_sf = st_as_sf(WS_data, coords= c("LONGITUDE","LATITUDE"), crs = sf::st_crs(4326) )%>%
      st_intersection(bbox)%>%st_transform(UTM20)
    
    # plot(st_geometry(WS_data_sf%>%filter(SPECIES_CD_1 == "922")))
    
    
    #filter out rare species of beaked whales
    WS_data_sf = WS_data_sf%>%dplyr::filter(SPECIES_CD_1 != "925" , SPECIES_CD_1 != "924", SPECIES_CD_1 != "7041",
                                              SPECIES_CD_1 !="7039")
    
    # Extract unique species names
    unique_species <- unique(WS_data_sf$COMMONNAME)
    
    # Determine the number of unique species
    num_species <- length(unique_species)

performKDE <- function(data_sf, unique_species, buffer_percent, sigma_func, output_dir = "output/tif/sights/") {

  shapefile_crs <- st_crs(data_sf)$wkt
  
  # Initialize a list to store plots
  raster_list <- list()
  global_min <- Inf
  global_max <- -Inf
  
  # Calculate coordinates for the entire dataset once
  coords <- st_coordinates(data_sf)
  
  # Determine the overall expanded extent
  x_range <- range(coords[, "X"])
  y_range <- range(coords[, "Y"])
  x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
  y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
  expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
  expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
  
  # Create a window object based on the expanded extent
  window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
  
  
  # Loop through each species 1rst time
  for(species in unique_species) {
    current_sf <- data_sf[data_sf$COMMONNAME == species, ]
    
    if(nrow(current_sf) > 0) {
      # coords <- st_coordinates(current_sf)
      # 
      # # Expand the extent
      # x_range <- range(coords[, "X"])
      # y_range <- range(coords[, "Y"])
      # x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
      # y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
      # expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
      # expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
      # 
      # 
      # #make data into pp object
      # window <- owin(xrange = x_range, yrange = y_range)
      # points_ppp <- ppp(x = coords[,"X"], y = coords[,"Y"], window = window)
      
      species_coords <- st_coordinates(current_sf)
      
      # Create ppp object within the precomputed window
      points_ppp <- ppp(x = species_coords[,"X"], y = species_coords[,"Y"], window = window)
      
      # Apply the bandwidth function to the point pattern object
      sigma_val <- sigma_func(points_ppp)
      
      # Kernel Density Estimation
      kd_result <- density.ppp(points_ppp, sigma = sigma_val, positive = T,
                               kernel = "gaussian", dimyx = 500, diggle = T)
      
      # Convert to SpatRaster
      raster_kd <- rast(kd_result)
      # Set the CRS for the raster
      crs(raster_kd) <- shapefile_crs
      
      # Normalize the raster values by the maximum value within the grid for the species
      max_val <- max(values(raster_kd), na.rm = TRUE)
      if (!is.na(max_val) && max_val > 0) {
        raster_kd <- raster_kd / max_val
      }
      # Store the raster in a list
      raster_list[[species]] <- raster_kd
      
      # Extract min/max from the raster
      local_min <- min(values(raster_kd), na.rm = TRUE)
      local_max <- max(values(raster_kd), na.rm = TRUE)
      
      # Update global min/max
      global_min <- min(global_min, local_min)
      global_max <- max(global_max, local_max)
      
      
      # Save the raster to a file
      raster_filename <- paste0(output_dir, species, "_sights_KDE.tif")
      writeRaster(raster_kd, filename = raster_filename, overwrite = T)
    }
  }


  # Second pass: Create plots using stored rasters and global min/max
  plot_list <- list()

  for(species in unique_species) {
    raster_kd <- raster_list[[species]]
    if(!is.null(raster_kd)) {
      # Convert raster to data frame for ggplot
      kde_df <- as.data.frame(raster_kd, xy = TRUE) %>%
        # mutate(sqrtlyr.1 = sqrt(lyr.1))%>%
        filter(!is.na(lyr.1))

      # Create ggplot using global range for colour scale
      p <- ggplot(kde_df, aes(x = x, y = y, fill = lyr.1)) +
        geom_tile() +
        scale_fill_viridis_c(limits = c(global_min, global_max), label = scales::scientific) +
        theme_minimal() +
        theme(axis.title= element_blank(),
              axis.ticks = element_blank(),
              axis.text = element_blank(), axis.line = element_blank()) +
        # legend.key.size = unit(0.5, "lines"),  # Adjust legend key size
        # legend.text = element_text(size = 8)  # Adjust legend text size
        labs(fill = "", title =  species) +
        annotate("text", x = Inf, y = Inf, label = paste("Bandwidth:", round(sigma_val, 0)),
                 hjust = 1.5, vjust = 2, size = 3, color = "white")+
        coord_equal()


      # Add the plot to the list
      plot_list[[species]] <- p

    } else {
      cat("No data found for", species, "\n")
    }
  }
  # Combine all plots using patchwork
  plot_combined <- wrap_plots(plot_list, guides = "collect")
  return(plot_combined)

}


      

#run KDE function
performKDE(data_sf = WS_data_sf, unique_species = unique_species, buffer_percent = .2, 
           sigma_func = bw.ppl)


