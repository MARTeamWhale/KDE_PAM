# creat KDE maps similar to ArcGIS...
options(scipen = 999)
library(stars)
library(mapsf)
library(stringr)
library(terra)
library(dplyr)
library(viridis)

#Create shapefiles of Quantile bins for RASTER KDE --------

processKDE <- function(tif_file, output_dir = "output/shapes/") {

  # Load and plot KDE raster
  kde_raster <- rast(tif_file)
  # Extract the CRS from the raster
  # Extract the CRS from the raster and convert it to an sf-compatible format
  raster_crs <- crs(kde_raster)  
  plot(kde_raster)
  

  # Process the raster, change name of data layer to be standard across tifs
  kde_stars <- st_as_stars(kde_raster, crs = crs(raster_crs)) %>%
    rename(KDE = names(st_as_stars(kde_raster))[1]) %>%
    mutate(KDE = ifelse(KDE == 0, NA, KDE))
  
  # Set the CRS for the stars object
  st_crs(kde_stars) <- raster_crs
  
  # Calculate quantile breaks
  breaks_quant <- mf_get_breaks(x = c(min(kde_stars$KDE, na.rm = TRUE), kde_stars$KDE), breaks = "q6")
  breaks_quant <- breaks_quant[c(2,3,4,5,6,7)]
  
  # Create quantile classes with NA values
  kde_quant <- mutate(kde_stars, breaks = cut(KDE, breaks_quant)) 
  kde_quant_rast <- rast(kde_quant[2])
  plot(kde_quant_rast)
  
  # Convert raster to polygon and set the CRS from the raster
  kde_sf <- st_make_valid(st_as_sf(terra::as.polygons(kde_quant_rast), merge = TRUE))
  kde_sf <- st_set_crs(kde_sf, raster_crs)
  
  plot(st_geometry(kde_sf))
  
  # Clean and process the data
  kde_sf <- kde_sf %>%
    mutate(quants = str_replace_all(breaks, "[\\(\\]]", "" )) %>%
    mutate(quants = as.numeric(gsub(",.*", "", quants))) %>%
    mutate(quants = factor(quants)) %>%
    mutate(Quantile = c(0.05, 0.25, 0.50, 0.75, 0.95))
  
  # Extract species name from file path
  species_name <- str_extract(basename(tif_file), "^[^_]+")
  
  # Save the shapefile
  output_path <- paste0(output_dir, species_name, "_KDE_Quant.shp")
  write_sf(kde_sf, output_path, overwrite = TRUE)
}

# # Example usage for a single file
# processKDE("output/tif/Bb_kernel_density.tif")


# List of .tif files
tif_files <- list.files("output/tif/", pattern = "\\.tif$", full.names = TRUE)

# Process each .tif file
for (tif_file in tif_files) {
  processKDE(tif_file)
}

######
#PLOT---------
######
# source("Scripts/load_basemap_shapes.R")

species_names_map <- c(
  Ha = "Northern bottlenose whale",
  Mb = "Sowerby's beaked whale",
  Bb = "Sei whale",
  Ma = "Humpback whale",
  Bp = "Fin whale",
  Bm = "Blue whale"
 )

plotKDEMaps <- function(shapefile_dir, 
                       land, contour_data, 
                       beaked_stns, baleen_stns,
                       output_dir, predpal = "mako", size = c(11, 8.5)) {

  # List all shapefiles in the directory
  shapefiles <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)
  
  
  
  
  # Loop through each shapefile and plot
  for (shapefile in shapefiles) {
    # Read the shapefile
    kde_sf_data <- st_read(shapefile)
    
    # Extract species name from the filename
    species_code <- str_extract(basename(shapefile), "^[^_]+")
    
    # Choose the appropriate stations data based on species
    if(species_code %in% c("Ha", "Mb")) {
      current_detects_data <- beaked_stns
    } else {
      current_detects_data<- baleen_stns
    }
    
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
    xlims <- c(xmin-10000, xmax+10000)
    ylims <- c(ymin- 10000, ymax +10000)
    
    
    

  # Define palette
  pal <- rev(get(predpal)(7))
  pal <- c(pal[c(2,3,4,5,6)])
  
  # Plot gg------
  gg_map <- ggplot() +
    theme_bw() +
    # geom_sf(data = contour_data %>% dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), col = "grey50", linewidth = 0.2) +
    geom_sf(data = kde_sf_data, aes(fill = quants), col = NA, alpha = .9, na.rm = TRUE) +
    geom_sf(data = current_detects_data, col = "black", shape = 24, fill = "yellow", size = 2, alpha = .5) +
    # Add station labels
    geom_text(data = current_detects_data, aes(x = coords[,1], y = coords[,2]+5000, label = site), 
              color = "white", size = 3, vjust = 0, position = position_dodge(0.9)) +  # Adjust as needed
    geom_sf(data = land, color = NA, fill = "grey50") +
    coord_sf(lims_method = "orthogonal", xlim = xlims, ylim = ylims, crs = UTM20, expand = T) +
    # Add species name as title 
    labs(title = paste("KDE Map for", species_name)) +
    
     ylab("") + 
    xlab("") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.1, 0.9), 
          legend.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), plot.margin = margin(0, 0, 0, 0), 
          plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title = element_blank()) +
    scale_fill_manual(values = pal, na.value = NA, na.translate = FALSE, labels = c("0.05","0.25", "0.50","0.75","0.95"), 
                      name = "", guide_legend(title.position = "bottom")) +
    annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6, bar_cols = c("grey40", "white"))
  
  # Output file name
  output_file <- paste0(output_dir, basename(shapefile), ".png")
  
  # Save the plot
  ggsave(output_file, gg_map, width = size[1], height = size[2], units = "in")
  }
}

plotKDEMaps(shapefile_dir = "output/shapes/", land = landUTM, contour_data = cont_UTM, 
            beaked_stns = beaked_stns,
            baleen_stns = baleen_stns,
            output_dir = "output/FIGS/")


kde_sf_data <- st_read("output/shapes/Bb_KDE_Quant.shp")
