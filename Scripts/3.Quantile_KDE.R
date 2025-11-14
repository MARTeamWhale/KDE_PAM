# create shapefiles of kde density rasters by binning into quantiles similar to ArcGIS...
#change input to change source maps


# create shapefiles of kde density rasters by binning into quantiles
#change input to change source maps
options(scipen = 999)
pacman::p_load(sf,stars, dplyr, stringr, terra, classInt)
input = "output/tif/"
output = "output/shapes/"
# Create shapefiles of Quantile bins for RASTER KDE --------
processKDE <- function(tif_file, output_dir = output) {
  
  # Load KDE raster
  kde_raster <- rast(tif_file)
  
  # Extract the CRS from the raster
  raster_crs <- crs(kde_raster)  
  
  # Set values <= 0 to NA
  kde_raster[kde_raster <= 0] <- NA
  
  # Define quantile probabilities
  probs <- c(.50,0.75, 0.85, 0.9, 0.95, 1)
  
  # Calculate quantile breaks (excluding NA values)
  kde_values <- values(kde_raster, mat = FALSE)
  quantiles <- quantile(kde_values, probs = probs, type = 6, na.rm = TRUE)#similar to arcgis cutoffs
  
  # Create labels for the quantile bins
  quant_labels <- paste0(probs[-length(probs)] * 100, "%")
  
  # Cut the values into quantile classes with labels
  kde_binned_values <- cut(kde_values, 
                           breaks = quantiles, 
                           labels = quant_labels,
                           include.lowest = TRUE)
  
  # Create a new raster with binned values
  kde_quant_rast <- kde_raster
  values(kde_quant_rast) <- as.numeric(kde_binned_values)
  
  # Convert raster to polygon
  kde_polys <- as.polygons(kde_quant_rast, dissolve = TRUE)
  kde_sf <- st_as_sf(kde_polys)
  kde_sf <- st_make_valid(kde_sf)
  kde_sf <- st_set_crs(kde_sf, raster_crs)
  
  # Debug: Check column names
  cat("Column names after conversion:", names(kde_sf), "\n")
  cat("First column values:", kde_sf[[1]], "\n")
  
  # Add quantile labels - the first column should be the bin number
  names(kde_sf)[1] <- "quant_level"
  
  # Map the numeric levels back to percentage labels
  kde_sf <- kde_sf %>%
    mutate(Quantile = quant_labels[quant_level])
  
  # Debug: Check if Quantile column was created
  cat("Quantile values:", unique(kde_sf$Quantile), "\n")
  
  # Extract species name from file path
  species_name <- str_extract(basename(tif_file), "(?<=KDE_)[^\\.]+")
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save the shapefile
  output_path <- file.path(output_dir, paste0( "KDE_", species_name,".shp"))
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  cat("Processed:", species_name, "\n\n")
  
  return(invisible(kde_sf))
}
# Process all .tif files------
tif_files <- list.files(input, pattern = "\\.tif$", full.names = TRUE)
# Process each .tif file with error handling
for (tif_file in tif_files) {
  tryCatch({
    processKDE(tif_file)
  }, error = function(e) {
    cat("Error processing", basename(tif_file), ":", conditionMessage(e), "\n")
  })
}

