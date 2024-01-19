#
#Script to read in PAM data on detections and export KDE surfaces of presence
#January 11, 2024

# Load necessary libraries
# install.packages("spatstat")
options(scipen = 999)

library(terra)
library(spatstat)
library(sf)
library(dplyr)
library(patchwork)

#projections------
UTM20 <- crs("+init=epsg:32620") # CODE FOR UTM Zone 20

#uses spatstat


performKDE <- function(data_sf, species_list, weight_col, buffer_percent, sigma_func, output_dir = "output/tif/") {
  # Preparing the plot window
  # par(mfrow = c(ceiling(length(species_list) / 2), 2))  # Adjust the layout as needed
  # Extract the CRS from the shapefile
  shapefile_crs <- st_crs(data_sf)$wkt

  # Initialize a list to store plots
  raster_list <- list()
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
      
      # Apply the bandwidth function to the point pattern object
      sigma_val <- sigma_func(points_ppp)
      
      # Kernel Density Estimation
      kd_result <- density.ppp(points_ppp, sigma = sigma_val, positive = T,
                               kernel = "gaussian", weights = marks(points_ppp), dimyx = 1000, diggle = T)
      
      # Convert to SpatRaster
      raster_kd <- rast(kd_result)
      # Set the CRS for the raster
      crs(raster_kd) <- shapefile_crs
      
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
        mutate(sqrtlyr.1 = sqrt(lyr.1))%>% filter(!is.na(lyr.1))
      
      # Create ggplot using global range for colour scale
      p <- ggplot(kde_df, aes(x = x, y = y, fill = sqrtlyr.1)) + 
        geom_tile() +
        scale_fill_viridis_c(trans = "sqrt", limits = sqrt(c(global_min, global_max)), label = scales::scientific) + 
        theme_minimal() + 
        theme(axis.title= element_blank(), 
                                axis.ticks = element_blank(),
                                 axis.text = element_blank(), axis.line = element_blank()) +
                                # legend.key.size = unit(0.5, "lines"),  # Adjust legend key size
                                # legend.text = element_text(size = 8)  # Adjust legend text size
        labs(fill = "", title = paste0("KDE for ", species)) +
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
baleen_PA <- read.csv('input/baleen_presence_days_laura.csv')


baleen_sf = baleen_PA%>%st_as_sf(coords = c("longitude", "latitude"), crs =4326 )%>%st_transform(crs = UTM20)


# List of Baleen species-----
#Bb = Sei, Bm = Blue, Bp = Fin, Mn = humpback
species_list <- c("Bb", "Bm", "Bp", "Mn")

#run KDE function
performKDE(data_sf = baleen_sf, species_list = species_list, 
           weight_col = "proportion_det", sigma_func = bw.diggle,  buffer_percent = .2)


#read in NBW & SBW data for comparison

bw_sf = read_sf("input/DOY.shp" )%>%st_transform(crs = UTM20)

# List of Beaked species-----
species_list <- c("Ha", "Mb")

#run KDE function
performKDE(data_sf = bw_sf, species_list = species_list, sigma_func = bw.diggle,
           weight_col = "pr_DOY_", buffer_percent = .2)




#       #TESTS--------

# raster_kd <- rast("output/tif/Bb_kernel_density.tif")
# plot(raster_kd)
# 
# kde_df <- as.data.frame(raster_kd, xy = TRUE)%>%dplyr::filter(!is.na(lyr.1))
# str(kde_df)
# summary(kde_df)

