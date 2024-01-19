# creat KDE maps similar to ArcGIS version used in NBW habitat...
######
#PLOT---------
######
# source("Scripts/load_basemap_shapes.R")

species_names_map <- c(
  Ha = "Northern bottlenose whale",
  Mb = "Sowerby's beaked whale",
  Bb = "Sei whale",
  Mn = "Humpback whale",
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
    xlims <- c(xmin-100000, xmax+100000)
    ylims <- c(ymin- 75000, ymax +75000)
    
    
    
    
    # Define palette
    pal <- rev(get(predpal)(7))
    pal <- c(pal[c(2,3,5)])
    
    # scales::show_col(pal)
    
    # Plot gg------
    gg_map <- ggplot() +
      theme_bw() +
      geom_sf(data = contour_data %>% dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), col = "grey50", linewidth = 0.2) +
      geom_sf(data = kde_sf_data, aes(fill = Quantile), col = NA, alpha = .9, na.rm = TRUE) +
      geom_sf(data = current_detects_data, col = "black", shape = 24, fill = "yellow", size = 2, alpha = .5) +
      # Add station labels
      geom_text(data = current_detects_data, aes(x = coords[,1], y = coords[,2]+5000, label = site), 
                color = "white", size = 3, vjust = 0, position = position_dodge(0.9)) +  # Adjust as needed
      geom_sf(data = land, color = NA, fill = "grey50") +
      coord_sf(lims_method = "orthogonal", xlim = xlims, ylim = ylims, crs = UTM20, expand = T) +
      # Add species name as title 
      labs(title = paste("Acoustic Presence Map for", species_name)) +
      
      ylab("") + 
      xlab("") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.1, 0.9), 
            legend.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), plot.margin = margin(0, 0, 0, 0), 
            plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title = element_blank()) +
      scale_fill_manual(values = pal, na.value = NA, na.translate = FALSE, labels = c( "0.50","0.75","0.95"), 
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
