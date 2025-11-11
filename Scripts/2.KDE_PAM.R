#
#Script to read in PAM data on detections and export KDE surfaces of presence
#January 11, 2024
#update 2025

# Load necessary libraries

options(scipen = 999)

pacman::p_load(terra,spatstat,sf,dplyr,patchwork, ggplot2)

filepath = 'input/2025/baleen_presence_days_laura_2025.csv'


#projections------
UTM20 <- "EPSG:32620" # CODE FOR UTM Zone 20

## For each species code in species_list code:
      # Extracts that species’ points and computes a buffered bounding box.
      # Builds a spatstat point pattern within that box.
      # Pulls weights from weight_col and multiplies them by 10,000.
      # Chooses a kernel bandwidth with a user-supplied function (bw.diggle).
      # Runs a Gaussian kernel density estimate on a 1000×1000 grid with Diggle edge
      # correction.
      # Converts the result to a SpatRaster, stamps the CRS, normalizes values by that
      # species’ own max, writes a GeoTIFF, and stores it in memory.
      # After looping, it plots one panel per species with a common color scale and
      # combines them with patchwork.

performKDE <- function(data_sf, species_list, weight_col, buffer_percent, sigma_func, output_dir = "output/tif/") {

  shapefile_crs <- st_crs(data_sf)$wkt

  # Initialize a list to store plots
  raster_list <- list()
  sigma_store <- list() 
  global_min <- Inf
  global_max <- -Inf
  
      # Loop through each species 1rst time
  for(species in species_list) {
    current_sf <- data_sf[data_sf$species == species, ]

    if(nrow(current_sf) > 0) {
      coords <- st_coordinates(current_sf)
      
      # Expand the extent
      x_range <- range(coords[, "X"])
      y_range <- range(coords[, "Y"])
      x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
      y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
      expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
      expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
      
      
      #make data into pp object
       window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
       points_ppp <- ppp(x = coords[,"X"], y = coords[,"Y"], window = window)
      
      # Assign weights
      weights <- current_sf[[weight_col]]*10000
      marks(points_ppp) <- data.frame(weights)
      
      # Apply the bandwidth function to the point pattern object and store it
      sigma_val <- sigma_func(points_ppp)
      sigma_store[[species]] <- sigma_val
      
      
      # Kernel Density Estimation
      kd_result <- density.ppp(points_ppp, sigma = sigma_val, positive = T,
                               kernel = "gaussian", weights = marks(points_ppp), dimyx = 1000, diggle = T)
      
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
      raster_filename <- paste0(output_dir, species, "_kernel_density.tif")
      writeRaster(raster_kd, filename = raster_filename, overwrite = T)
    }
  }
  
  
  # Second pass: Create plots using stored rasters and global min/max
  plot_list <- list()
  
  for(species in species_list) {
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
        labs(fill = "", title = paste0("Relative density of ", species)) +
        annotate("text", x = Inf, y = Inf, label = paste("Bandwidth:", round(sigma_val, 0)),
                 hjust = 1.1, vjust = 2, size = 3, color = "white")+
        coord_equal() 
      
      
      # Add the plot to the list
      plot_list[[species]] <- p
     
    } else {
      cat("No data found for", species, "\n")
    }
  }
  # Combine all plots using patchwork
  plot_combined <- wrap_plots(plot_list, ncol = 2, guides = "collect")
  return(plot_combined)
  
}

# sigma_func = bw.scott outputs 2 d values so need to change code annotation to:
# annotate("text", x = Inf, y = Inf, label = paste("Bandwidth:",round(sigma_val[1], 0), "/", round(sigma_val[2], 0)),

# Reading baleen whale data------
baleen_PA <- read.csv(filepath)

# # Normalize the proportions within each group
# baleen_PA <- baleen_PA %>%
#   group_by(species) %>%
#   mutate(sumProp = sum(proportion_det),
#     NormalizedProportion = proportion_det / sumProp
#   ) %>%
#   ungroup() # Always a good practice to ungroup after you're done


baleen_sf = baleen_PA%>%st_as_sf(coords = c("longitude", "latitude"), crs =4326 )%>%st_transform(crs = UTM20)


# List of Baleen species-----
#Bb = Sei, Bm = Blue, Bp = Fin, Mn = humpback, Ba = Minke, Eg = NARW
baleen_PA%>%group_by(species)%>%summarise(count = n())

species_list <- c("Bb", "Bm", "Bp", "Mn", "Eg","Ba")

#run KDE function
performKDE(data_sf = baleen_sf, species_list = species_list, 
           weight_col = "proportion_det", sigma_func = bw.diggle,  buffer_percent = .2)


#       QA TESTS--------

raster_kd <- rast("output/tif/Bb_kernel_density.tif")
# 
# kde_df <- as.data.frame(raster_kd, xy = TRUE)%>%dplyr::filter(!is.na(lyr.1))
# str(kde_df)
# summary(kde_df)

