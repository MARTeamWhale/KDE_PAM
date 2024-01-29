# create shapefiles of kde density rasters by binning into quantiles similar to ArcGIS...

options(scipen = 999)
library(stars)
library(mapsf)
library(stringr)
library(terra)
library(dplyr)
library(viridis)
library(classInt)

#Create shapefiles of Quantile bins for RASTER KDE --------



processKDE <- function(tif_file, output_dir = "output/shapes/") {

  # Load and plot KDE raster
  kde_raster <- rast(tif_file)
  plot(kde_raster)
  
  # Extract the CRS from the raster
  raster_crs <- crs(kde_raster)  
  
 
  kde_raster[kde_raster <= 0] <- NA
  
  # get all data values
  kde_values <- values(kde_raster)
  
  # Calculate quantile breaks (excluding NA values)
  probs = c(0.5, 0.75, 0.8, 0.85,0.9, 0.95, 1)
  quantiles <- quantile(kde_values, probs = probs, type = 6, na.rm = TRUE)
  
  # Cut the values into quantile classes
  kde_binned_values <- cut(kde_values, breaks = quantiles, include.lowest = TRUE)
  
  # Create a new raster with these binned values
  kde_quant_rast <- kde_raster
  values(kde_quant_rast) <- kde_binned_values
  
  plot(kde_quant_rast)
  
  # Convert raster to polygon and set the CRS from the raster
  kde_sf <- st_make_valid(st_as_sf(terra::as.polygons(kde_quant_rast), merge = TRUE))
  kde_sf <- st_set_crs(kde_sf, raster_crs)
  
  plot(st_geometry(kde_sf))  # Plotting for visualization
  
  
  # # Clean and process the data
  names(kde_sf)[1]<-"quants"
  
  kde_sf <- kde_sf %>%
    mutate(quants = str_replace_all(quants, "[\\(\\]\\[]", "" ))%>%
    mutate(quants = as.numeric(gsub(",.*", "", quants))) %>%
    mutate(quants.1 = as.numeric(factor(quants))) %>%
    mutate(Quantile = case_when(
      quants.1 == 1  ~ paste(probs[1]*100, "%"),
      quants.1 == 2  ~ paste(probs[2]*100, "%"),
      quants.1 == 3 ~ paste(probs[3]*100, "%"),
                            quants.1 == 4 ~ paste(probs[4]*100, "%"),
      quants.1 == 5 ~ paste(probs[5]*100, "%"),
      quants.1 == 6 ~ paste(probs[6]*100, "%"),
      TRUE                       ~ "Other"
    ))

  # Extract species name from file path
  species_name <- str_extract(basename(tif_file), "^[^_]+")
  
  # Save the shapefile
  output_path <- paste0(output_dir, species_name, "_KDE_Quant.shp")
  write_sf(kde_sf, output_path, overwrite = TRUE)
}

        
        # # Example usage for a single file---------
        processKDE("output/tif/Bb_kernel_density.tif")
        
        
        # List of .tif files------
        tif_files <- list.files("output/tif/", pattern = "\\.tif$", full.names = TRUE)
        
        # Process each .tif file
        for (tif_file in tif_files) {
          processKDE(tif_file)
        }





#test---------
Ha = rast("output/tif/Ha_kernel_density.tif")


plot(Ha)
kde_raster = Ha
kde_raster[kde_raster <= 0] <- NA

# get all data values
kde_values <- values(kde_raster)

# Calculate quantile breaks (excluding NA values)
quantiles <- quantile(kde_values, probs = c(0.5, 0.75, 0.95, 1), type = 6, na.rm = TRUE)

# Cut the values into quantile classes
kde_binned_values <- cut(kde_values, breaks = quantiles, include.lowest = TRUE)

# Create a new raster with these binned values
kde_quant_rast <- kde_raster
values(kde_quant_rast) <- kde_binned_values

plot(kde_quant_rast)
# Convert raster to polygon
kde_sf <- st_as_sf(as.polygons(kde_quant_rast, values = T), merge = TRUE)

plot(st_geometry(kde_sf))

# # Clean and process the data
names(kde_sf)[1]<-"quants"

kde_sf <- kde_sf %>%
  mutate(quants = str_replace_all(quants, "[\\(\\]\\[]", "" ))%>%
  mutate(quants = as.numeric(gsub(",.*", "", quants))) %>%
  mutate(quants = factor(quants)) 
