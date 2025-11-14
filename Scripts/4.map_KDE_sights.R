# creat KDE maps for sightings...
######
#PLOT---------
######
 # source("Scripts/1.load_basemap_shapes.R")


plotKDEMaps <- function(shapefile_dir, 
                        land, contour_data, 
                        output_dir, predpal = "mako", size = c(11, 8.5)) {
  
  # List all shapefiles in the directory
  shapefiles <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)
  
  # Loop through each shapefile and plot
  for (shapefile in shapefiles) {
    # Read the shapefile
    kde_sf_data <- st_read(shapefile)
    
    # Extract species name from the filename
    species <- str_extract(basename(shapefile), "^(.*?)_KDE_Quant")
    
    # Extract the bounding box (same as KDE_sights)
    bbox <- st_bbox( c(xmin = -67.5,ymin = 40, xmax = -55, ymax =47.5 ), crs = st_crs(4326))%>% st_as_sfc()%>% #turns the bounding box into a sfc object, that just describes a specific geometry
      st_sf()%>%st_transform(UTM20)
    
    bbox = st_bbox(bbox)
    
    # Extract individual limits
    xmin <- bbox["xmin"]
    xmax <- bbox["xmax"]
    ymin <- bbox["ymin"]
    ymax <- bbox["ymax"]
    
    # Use these limits for xlims and ylims
    xlims <- c(xmin, xmax)
    ylims <- c(ymin, ymax)
    
    
    # Define palette
    pal <- rev(get(predpal)(6))
    # pal <- c(pal[c(2,3,5)])
    
    # scales::show_col(pal)
    
    #define labels
    unique_quantiles <- sort(unique(kde_sf_data$Quantil))
    
    # Plot gg------
    gg_map <- ggplot() +
      theme_bw() +
      geom_sf(data = contour_data %>% dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), col = "grey50", linewidth = 0.2) +
      geom_sf(data = kde_sf_data, aes(fill = Quantil), col = NA, alpha = .9, na.rm = TRUE) +
      
      geom_sf(data = land, color = NA, fill = "grey50") +
      coord_sf(lims_method = "orthogonal", xlim = xlims, ylim = ylims, crs = UTM20, expand = T) +
      # Add species name as title 
      labs(title = paste("Relative sighting density of", species)) +
      
      ylab("") + 
      xlab("") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = c(0.1, 0.9), 
            legend.background = element_rect(fill = NA), legend.key = element_rect(fill = NA), plot.margin = margin(0, 0, 0, 0), 
            plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5), axis.title = element_blank()) +
      scale_fill_manual(values = pal, na.value = NA, na.translate = FALSE, labels = unique_quantiles, 
                        name = "", guide_legend(title.position = "bottom")) +
      annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6, bar_cols = c("grey40", "white"))
    
    # Output file name
    output_file <- paste0(output_dir, basename(shapefile), ".png")
    
    # Save the plot
    ggsave(output_file, gg_map, width = size[1], height = size[2], units = "in")
  }
}

plotKDEMaps(shapefile_dir = "output/shapes/sights/", land = landUTM, contour_data = cont_UTM, 
            output_dir = "output/FIGS/sights/")
