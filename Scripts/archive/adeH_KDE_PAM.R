#trying another kernel function does not work as well due to being in SP

pacman::p_load(adehabitatHR, sf, ggplot2, patchwork, sp)

filepath = 'input/2025/baleen_presence_days_laura_2025.csv'

adehab_KDE <- function(data_sf, species_list, weight_col, buffer_percent, output_dir = "output/tif/") {
  # Initialize a list to store plots
  plot_list <- list()
  
  # Loop through each species
  for(species in species_list) {
    current_sf <- data_sf[data_sf$species == species, ]
    
    if(nrow(current_sf) > 0) {
      # Extract coordinates and data
      

      # Convert sf object to SpatialPointsDataFrame with only species and weight columns
      coordinates <- st_coordinates(current_sf)
      data_for_sp <- data.frame(weight = current_sf[[weight_col]],
                                stringsAsFactors = FALSE)
      sp_data <- SpatialPointsDataFrame(coordinates, data_for_sp)
      
      # Estimate the utilization distribution (UD) using kernelUD
      ud <- kernelUD(sp_data, h = "href")
      
      # Get the raster of the UD
      raster_ud <- rast(ud)
      
      # Optionally, save the raster to a file
      raster_filename <- paste0(output_dir, species, "_kernel_density.tif")
      writeRaster(raster_ud, filename = raster_filename, overwrite = TRUE)
      
      # Convert raster to data frame for ggplot
      kde_df <- as.data.frame(raster_ud, xy = TRUE)
      
      # Create ggplot
      p <- ggplot(kde_df, aes(x = x, y = y, fill = layer)) + 
        geom_tile() +
        scale_fill_viridis_c() + 
        theme_minimal() + 
        theme(axis.title = element_blank(), 
              axis.ticks = element_blank(),
              axis.text = element_blank(),
              legend.key.size = unit(0.5, "lines"),
              legend.text = element_text(size = 8)) +
        labs(fill = "Density", title = paste0("KDE for ", species)) +
        coord_equal()
      
      # Add the plot to the list
      plot_list[[species]] <- p
    } else {
      cat("No data found for", species, "\n")
    }
  }
  
  # Combine all plots using patchwork
  plot_combined <- wrap_plots(plot_list, ncol = 2, panel_spacing = unit(2, 'mm'))
  return(plot_combined)
}

# Reading baleen whale data------
baleen_PA <- read.csv(filepath)


baleen_sf = baleen_PA%>%st_as_sf(coords = c(6,7), crs =4326 )%>%st_transform(crs = UTM20)


# List of Baleen species-----
#Bb = Sei, Bm = Blue, Bp = Fin, Mn = humpback
species_list <- c("Bb", "Bm", "Bp", "Mn")

#run KDE function
adehab_KDE(data_sf = baleen_sf, species_list = species_list, 
           weight_col = "proportion_det", buffer_percent = .2)



#test
# Create SpatialPointsDataFrame

coords <- coordinates(baleen_PA[,6:7])
sp_data <- SpatialPointsDataFrame(coords = coords, baleen_PA, proj4string = CRS(UTM20))
ud <- kernelUD(sp_data, id = species, h = "href")


# kde estimation in adeH

