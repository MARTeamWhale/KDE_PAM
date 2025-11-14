# Unified KDE mapping for both PAM and sighting data
######
# PLOT---------
######
# depends on ("Scripts/load_basemap_shapes.R")

# Species name mappings for PAM data (abbreviations)
species_names_pam <- c(
  Bb = "Sei whale",
  Mn = "Humpback whale",
  Bp = "Fin whale",
  Bm = "Blue whale",
  Ba = "Minke whale",
  Eg = "Right Whale"
)

# Create point shapefile from records
# INPUTS NEED TO BE CSV files, file names are specified in quotations here:-----
input_file <-"WSDB_DFO_Sep2025.csv"
species = "speciesCodes.csv" 

whale_data =  read_csv(here("input/2025/sights/", input_file), col_types = cols(.default ="c"))

# check species, add common names and scientific names based on codes-----
SP_data <- read_csv(here("input", species), col_types = cols(.default ="c"))
WS_data = left_join(whale_data, SP_data)
# bound boxes for study area -----
bbox = Bound_boxB

WS_data_sf = st_as_sf(WS_data, coords= c("LONGITUDE","LATITUDE"), crs = sf::st_crs(4326) )%>%
  st_intersection(bbox)%>%st_transform(UTM20)
# Ensure same CRS
WS_data_sf <- st_make_valid(WS_data_sf)
landUTM    <- st_make_valid(landUTM)
WS_data_sf <- st_transform(WS_data_sf, st_crs(landUTM))



# Define vector of cetacean species codes (excluding NS categories and seals)
cetacean_codes <- c(
  921, 922, 923, 925, 931, 932, 933, 934, 935, 936, 937, 938,
  7019, 7020, 7021, 7022, 7023, 7024, 7025, 7026, 7027, 7028,
  7029, 7031, 7033, 7034, 7035, 7037, 7038, 902
)

# 3) Filter: cetaceans -> not on land -> species with >50 records
WS_data_sf <- WS_data_sf %>%
  filter(SPECIES_CD %in% cetacean_codes) %>%
  # keep points that DO NOT intersect land
  filter(lengths(st_intersects(., landUTM)) == 0) %>%
  add_count(SPECIES_CD) %>%
  filter(n > 50) %>%
  select(-n)

plotKDEMaps <- function(shapefile_dir, 
                        land, 
                        baleen_stns = NULL,
                        sighting_data = NULL,
                        output_dir, 
                        predpal = "mako", 
                        size = c(11, 8.5),
                        exclude_quantiles = "50%") {  # NEW PARAMETER
  
  # List all shapefiles in the directory
  shapefiles <- list.files(shapefile_dir, pattern = "\\.shp$", full.names = TRUE)
  
  # Loop through each shapefile and plot
  for (shapefile in shapefiles) {
    # Read the shapefile
    kde_sf_data <- st_read(shapefile, quiet = TRUE)
    
    # Detect data type from filename
    filename <- basename(shapefile)
    is_pam <- str_detect(filename, "KDE_pam_")
    is_sighting <- str_detect(filename, "KDE_sights_")
    
    # Set up bounding box (same for both data types)
    bbox <- st_bbox(c(xmin = -67.5, ymin = 40, xmax = -55, ymax = 47.5), 
                    crs = st_crs(4326)) %>% 
      st_as_sfc() %>%
      st_sf() %>% 
      st_transform(UTM20)
    
    bbox <- st_bbox(bbox)
    
    xmin <- bbox["xmin"]
    xmax <- bbox["xmax"]
    ymin <- bbox["ymin"]
    ymax <- bbox["ymax"]
    
    xlims <- c(xmin, xmax)
    ylims <- c(ymin, ymax)
    
    # Extract species name and set title based on data type
    if (is_pam) {
      # PAM data processing
      species_code <- str_extract(filename, "(?<=KDE_pam_)[^\\.]+")
      species_name <- species_names_pam[species_code]
      
      # If species name is not found in the map, use the code as a fallback
      if (is.null(species_name) || is.na(species_name)) {
        species_name <- species_code
      }
      
      # Check if baleen_stns is provided
      if (is.null(baleen_stns)) {
        stop("baleen_stns must be provided for PAM data")
      }
      
      map_title <- paste("Acoustic Presence Map for", species_name)
      
    } else if (is_sighting) {
      # Sighting data processing
      raw_name <- str_extract(filename, "KDE_sights_(.+?)(?=\\.shp)", group = 1)
      
      # Split by hyphen to get animal type and species
      name_parts <- str_split(raw_name, "-", n = 2)[[1]]
      
      if (length(name_parts) == 2) {
        animal_type <- str_to_title(name_parts[1])  # e.g., "Whales", "Dolphins"
        species <- str_to_title(str_replace_all(name_parts[2], "[_-]", " "))  # e.g., "Atlantic Bottlenose"
        
        # Singularize animal type if needed (simple approach)
        animal_type <- str_replace(animal_type, "s$", "")
        
        species_name <- paste(species, animal_type)
      } else {
        # Fallback if format is different
        species_name <- str_to_title(str_replace_all(raw_name, "[_-]", " "))
      }
      
      map_title <- paste("Relative sighting density of", species_name)
      
    } else {
      cat("Warning: Cannot determine data type for", filename, ". Skipping.\n")
      next
    }
    
    # FILTER OUT EXCLUDED QUANTILES BEFORE PROCESSING
    kde_sf_data <- kde_sf_data %>%
      filter(!Quantil %in% exclude_quantiles)
    
    # Get unique quantiles and sort them properly
    unique_quantiles <- unique(kde_sf_data$Quantil)
    unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
    
    # Extract numeric values for sorting
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    
    # Sort by numeric values
    sort_order <- order(quantile_nums)
    unique_quantiles <- unique_quantiles[sort_order]
    
    # Check if we have valid quantiles
    if (length(unique_quantiles) == 0) {
      cat("Warning: No valid quantiles found after filtering. Skipping", filename, "\n")
      next
    }
    
    # Convert Quantile to ordered factor
    kde_sf_data$Quantil <- factor(kde_sf_data$Quantil, 
                                  levels = unique_quantiles, 
                                  ordered = TRUE)
    
    # Define palette - match the number of colors to number of quantiles
    n_quantiles <- length(unique_quantiles)
    pal <- rev(get(predpal)(n_quantiles))
    
    # Start building the plot
    gg_map <- ggplot() +
      theme_bw() +
      geom_sf(data = cont_UTM %>% 
                dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), 
              col = "grey50", linewidth = 0.2) +
      geom_sf(data = kde_sf_data, aes(fill = Quantil), col = NA, alpha = .9, na.rm = TRUE)
    
    # Add sighting points for sighting data
    if (is_sighting && !is.null(sighting_data)) {
      # Extract the raw name from filename for better matching
      raw_species <- str_extract(filename, "KDE_sights_(.+?)(?=\\.shp)", group = 1)
      
      cat("\nProcessing sighting map:", filename)
      cat("\nRaw species string from filename:", raw_species, "\n")
      
      # Normalize both the filename and data for matching
      # Remove all spaces, hyphens, and apostrophes, convert to uppercase
      normalize_name <- function(x) {
        str_to_upper(str_replace_all(x, "[ -']", ""))
      }
      
      normalized_filename <- normalize_name(raw_species)
      
      # Find matching species in the data
      species_sightings <- sighting_data %>%
        filter(normalize_name(COMMONNAME) == normalized_filename)
      
      cat("Number of matching sightings found:", nrow(species_sightings), "\n")
      
      if (nrow(species_sightings) > 0) {
        gg_map <- gg_map +
          geom_sf(data = species_sightings, col = "white", shape = 21, 
                  fill = "yellow", size = 1.5, alpha = 0.6, stroke = 0.5)
        cat("Added", nrow(species_sightings), "sighting points to map\n")
      } else {
        cat("WARNING: No matching sightings found for", species_name, "\n")
      }
    }
    
    # Add station points and labels ONLY for PAM data
    if (is_pam && !is.null(baleen_stns)) {
      
      cat("\n=== Processing PAM Station Labels ===\n")
      
      # Add truncated site column
      current_detects_data <- baleen_stns %>%
        mutate(site_short = case_when(
          site == "SABW" ~ "SBW",
          site == "ROBP" ~ "RBP",
          site == "ROBV" ~ "RBV",
          site == "GDSE" ~ "GUL",
          str_detect(site, "^MG") ~ "MGU",
          str_detect(site, "^GL") ~ "GUL",
          TRUE ~ substr(site, 1, 3)
        ))
      
      cat("Total detects:", nrow(current_detects_data), "\n")
      
      # Create unique stations (one per truncated site name)
      unique_stations <- current_detects_data %>%
        group_by(site_short) %>%
        slice(1) %>%
        ungroup()
      
      cat("Unique stations:", nrow(unique_stations), "\n")
      cat("Station names:", paste(unique_stations$site_short, collapse=", "), "\n")
      
      # Extract coordinates
      station_coords <- st_coordinates(unique_stations)
      cat("Coordinate matrix dimensions:", dim(station_coords), "\n")
      cat("Coordinate column names:", colnames(station_coords), "\n")
      
      # Create label dataframe using numeric indexing
      label_df <- data.frame(
        x = as.numeric(station_coords[, 1]),
        y = as.numeric(station_coords[, 2]) + 5000,
        label = as.character(unique_stations$site_short),
        stringsAsFactors = FALSE
      )
      
      cat("Sample label positions:\n")
      print(head(label_df, 3))
      
      # Add station points to map
      gg_map <- gg_map +
        geom_sf(data = current_detects_data, col = "black", shape = 24, 
                fill = "yellow", size = 2, alpha = 0.5)+
        geom_sf_text(data = unique_stations, 
                  aes(label = site_short), 
                  color = "white", size = 3, vjust = 0, 
                  check_overlap = FALSE)
      
      cat("Added", nrow(label_df), "text labels\n")
      cat("=== Done ===\n\n")
    }
    
    # Complete the plot
    gg_map <- gg_map +
      geom_sf(data = land, color = NA, fill = "grey50") +
      coord_sf(lims_method = "orthogonal", xlim = xlims, ylim = ylims, 
               crs = UTM20, expand = TRUE) +
      labs(title = map_title) +
      ylab("") + 
      xlab("")+
     
       # Continue theme settings
            theme(panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(), 
            legend.position = c(0.1, 0.9), 
            legend.background = element_rect(fill = NA), 
            legend.key = element_rect(fill = NA), 
            plot.margin = margin(5, 5, 5, 5),  # Increased margins
            plot.title = element_text(hjust = 0.5), 
            plot.subtitle = element_text(hjust = 0.5), 
            axis.title = element_blank()) +
      scale_fill_manual(values = pal, 
                        na.value = NA, 
                        na.translate = FALSE,
                        name = if(is_pam) "Quantile" else "",
                        labels = unique_quantiles,
                        drop = FALSE) +
      annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6, 
                       bar_cols = c("grey40", "white"))
    
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Output file name
    output_file <- file.path(output_dir, paste0(gsub("\\.shp$", "", basename(shapefile)), ".png"))
    
    # Save the plot
    ggsave(output_file, gg_map, width = size[1], height = size[2], units = "in")
    
    cat("Saved:", output_file, "\n")
  }
}

# Example usage:
# Run for both PAM and sighting data from the same directory
# Exclude the 50% quantile from plots
plotKDEMaps(shapefile_dir = "output/shapes/", 
            land = landUTM, 
            baleen_stns = baleen_stns,  # Only used for PAM data
            sighting_data = WS_data_sf, # Sighting points for sighting KDE maps
            output_dir = "output/FIGS/",
            exclude_quantiles = "50%")  # NEW: Exclude 50% quantile

# If you want to exclude multiple quantiles, pass a vector:
# exclude_quantiles = c("50%", "25%")