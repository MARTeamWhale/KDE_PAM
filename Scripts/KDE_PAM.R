#
#Script to read in PAM data on detections and export KDE surfaces of presence
#January 11, 2024

# Load necessary libraries
install.packages("spatstat")
library(terra)
library(spatstat)
library(sf)
library(dplyr)


performKDE <- function(data_sf, species_list, weight_col, buffer_percent, sigma_val, output_dir = "output/") {

    # Loop through each species
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
      window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
      points_ppp <- ppp(x = coords[,"X"], y = coords[,"Y"], window = window)
      
      # Assign weights
      weights <- current_sf[[weight_col]]
      marks(points_ppp) <- data.frame(weights)
      
      # Kernel Density Estimation
      kd_result <- density.ppp(points_ppp, sigma = sigma_val, weights = marks(points_ppp))
      
      # Convert to SpatRaster
      raster_kd <- rast(kd_result)
      
      # Save the raster to a file
      raster_filename <- paste0(output_dir, species, "_kernel_density.tif")
      writeRaster(raster_kd, filename = raster_filename, overwrite = T)
      
      # Optionally plot the KDE
      plot(raster_kd, main = paste0("KDE for ", species))
    } else {
      cat("No data found for", species, "\n")
    }
  }
}

# Reading baleen whale data------
baleen_PA <- read.csv('input/baleen_presence_days_laura.csv')

baleen_sf = baleen_PA%>%st_as_sf(coords = c("longitude", "latitude"), crs =4326 )%>%st_transform(crs = 2961)

# List of Baleen species-----
#Bb = Sei, Bm = Blue, Bp = Fin, Mn = humpback
species_list <- c("Bb", "Bm", "Bp", "Mn")

#run KDE function
performKDE(data_sf = baleen_sf, species_list = species_list, 
           weight_col = "proportion_det", buffer_percent = .2, sigma_val = 50000)


#read in NBW & SBW data for comparison

bw_sf = read_sf("input/DOY.shp" )%>%st_transform(crs = 2961)

# List of Beaked species-----
species_list <- c("Ha", "Mb")

#run KDE function
performKDE(data_sf = bw_sf, species_list = species_list, 
           weight_col = "pr_DOY_", buffer_percent = .2, sigma_val = 50000)