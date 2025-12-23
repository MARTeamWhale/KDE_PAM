Share#
# Unified script for KDE analysis of whale sightings and PAM data
# Updated 2025
#
# ...............
# DATASET SELECTION - CONFIGURE HERE------
# ...............

RUN_SIGHTINGS <- TRUE   # Set to TRUE to process sightings data
RUN_PAM <- TRUE         # Set to TRUE to process PAM data

# ..........pastedok i need to update this script to run a KDE on the aerial sightings data which has slightly different structure than WSDB and includes a weight based on cluster size. can you help?3:38 PMI'd be happy to help you update the script for aerial sightings data! To create the best solution, I need to understand a few things about your aerial data:

File location and format: Where is the aerial sightings data file located and what format is it in (CSV, Excel, etc.)?
  Column structure: What are the key column names, particularly:
  
  Latitude/longitude columns
Species identifier column
Cluster size column (for weighting)
Any date/time columns


Species handling:
  
  Does it use the same species codes/naming as WSDB, or does it have its own system?
  Should it use the same species list filtering approach (>50 records)?
  
  
  Integration approach: Would you like:
  
  A separate section (like the existing PAM and WSDB sections) controlled by a RUN_AERIAL <- TRUE flag?
  The aerial data to be combined with the WSDB sightings data?
  Both options available?
  
  
  Spatial filtering: Should the aerial data:
  
  Use the same common window from PAM stations?
  Be filtered to remove land intersections like WSDB?
  Use the same bounding box?
  
  
  
  Once you provide these details, I can modify the script to handle your aerial sightings data with appropriate cluster size weighting!deduplicated_sightings.csvcsvsightings csv will be in the "input/2025/aerial/" folder. file is attached for the data var infos. uses dif sp_code than WSDB, yes should use same filtering process and only keep records from 2015 onwards. for #4. I would like a control with tru, and an option to run independently as not sure it makes sense to run together with WSDB sights given the lack of cluster size info in that dataset. For #5 yes to all questions3:42 PMPerfect! I'll update your script to add an aerial sightings section. Here's the modified version:Updated KDE Analysis Script with Aerial DataCode #
# Unified script for KDE analysis of whale sightings and PAM data
# Updated 2025
#
# ...............
# DATASET SELECTION - CONFIGURE HERE------
# ...............

RUN_SIGHTINGS <- TRUE   # Set to TRUE to process WSDB sightings data
RUN_AERIAL <- TRI've updated your script with the following key changes for aerial sightings:
Key Updates:

New control flag: RUN_AERIAL <- TRUE at the top to independently control aerial data processing
Aerial sightings section that:

Reads CSV files from input/2025/aerial/ folder
Filters to 2015 onwards based on Date_time_loc
Uses Long_sight and Lat_sight for coordinates (the actual sighting location, not platform location)
Uses Sp_name_clean as the species identifier column
Applies the same spatial filtering (bounding box, remove land intersections, >50 records)
Weights by Clustsize (cluster size) - this is the key difference from WSDB
Removes records with missing or zero cluster sizes


Cluster size summary: Prints statistics showing mean, median, and max cluster sizes by species before running KDE
Same spatial extent: Uses the common window from PAM stations, ensuring all three datasets are comparable
Separate outputs: Creates files with aerial prefix (e.g., KDE_aerial_Humpbackwhale.tif)

The aerial section runs independently, so you can compare weighted aerial results against unweighted WSDB results and PAM detection proportions. Let me know if you need any adjustments!this code cat("Date range in aerial data:", min(aerial_data$Date, na.rm = TRUE), 
      "to", max(aerial_data$Date, na.rm = TRUE), "\n") produces "Date range in aerial data: 16659 to 19926 "4:57 PMThe dates are being printed as numeric values instead of actual dates. This happens when the date conversion isn't working properly. Let's fix this:Updated KDE Analysis Script with Aerial DataCode ∙ Version 2   # Process dates and filter to 2015 onwards
  aerial_data <- aerial_combined %>%
    mutate(Date = as_date(Date_time_loc)) %>%
    mutate(Year = year(Date)) %>%
    filter(Year >= 2015)
  
  cat("Date range in aerial data:", 
      format(min(aerialI've added format() to convert the dates to a readable format. However, if you're still seeing numeric values, it might mean the Date_time_loc column isn't in a format that as_date() can parse directly.
If the issue persists, you might need to specify the date format explicitly. Can you check what the Date_time_loc column looks like in your CSV? If it's in a specific format (like "YYYY-MM-DD HH:MM:SS" or "DD/MM/YYYY"), we can adjust the code to parse it correctly using:
rmutate(Date = as_date(Date_time_loc, format = "%Y-%m-%d %H:%M:%S"))
Let me know what format the dates are in and I can update it accordingly!date_time_utc should be used and it looks fine4:57 PMUpdated KDE Analysis Script with Aerial DataCode ∙ Version 3   # Process dates and filter to 2015 onwards
  aerial_data <- aerial_combined %>%
    mutate(Date = as_date(Date_time_utc)) %>%
    mutate(Year = year(Date)) %>%
    filter(Year >= 2015)
  
  cat("Date range in aerial data:", 
      format(min(aerialPerfect! I've updated the code to use Date_time_utc instead of Date_time_loc. This should now properly parse the dates and display them correctly.can we omit all unknown "" species? filter them out before processing4:58 PMUpdated KDE Analysis Script with Aerial DataCode ∙ Version 4   # Filter: not on land -> remove unknown species -> valid cluster size -> species with >50 records
aerial_sf <- aerial_sf %>%
  filter(lengths(st_intersects(., landUTM)) == 0) %>%
  filter(!is.na(Sp_name_clean) & Sp_name_clean != "" & Sp_name_Done! I've updated the filtering to remove any records where Sp_name_clean is NA, empty string (""), or "Unknown" before the other filters are applied. This ensures only identified species are processed.no it will need to be a filter to detect the word "Unknown" as there are mutltiple types of unknown...Unknown Dolphin/Porpoise / Dauphin/Marsouin inconnu, Common Dolphin / Dauphin commun, Harbour Porpoise / Marsouin commun, Unknown Whale/Cetacean / Baleine/Cetace inconnu, Risso's Dolphin / Dauphin de Risso, Fin Whale / Rorqual commun, Humpback Whale / Rorqual a bosse, Pilot Whale / Globicephale, Minke Whale / Petit rorqual, Atlantic White-sided Dolphin / Dauphin a flanc blanc, Grey Seal / Phoque gris, Unknown Seal / Phoque inconnu 
         >5:00 PMUpdated KDE Analysis Script with Aerial DataCode ∙ Version 5   # Filter: not on land -> remove unknown species -> valid cluster size -> species with >50 records
         aerial_sf <- aerial_sf %>%
           filter(lengths(st_intersects(., landUTM)) == 0) %>%
           filter(!is.na(Sp_name_clean) & Sp_name_clean != "") %>%
           fPerfect! I've updated the filter to use grepl() to detect any species name containing "Unknown" or "inconnu" (case-insensitive). This will catch all variations like "Unknown Dolphin/Porpoise", "Unknown Whale/Cetacean", "Phoque inconnu", etc., and exclude them from the analysis.#
# Unified script for KDE analysis of whale sightings and PAM data
# Updated 2025
#
# ...............
# DATASET SELECTION - CONFIGURE HERE------
# ...............

RUN_SIGHTINGS <- TRUE   # Set to TRUE to process WSDB sightings data
RUN_AERIAL <- TRUE      # Set to TRUE to process aerial sightings pastedok for this version of the script i am getting this error? === PROCESSING PAM DATA ===
Processing 6 species from PAM data
Error processing species Bb : number of rows of data frame != number of points 
Error processing species Bm : number of rows of data frame != number of points 
Error processing species Bp : number of rows of data frame != number of points 
Error processing species Mn : number of rows of data frame != number of points 
Error processing species Eg : number of rows of data frame != number of points 
Error processing species Ba : number of rows of data frame != number of points 
No data found for Bb 
No data found for Bm 
No data found for Bp 
No data found for Mn 
No data found for Eg 
No data found for Ba 
Error in unit(rep(0, TABLE_COLS * dims[2]), "null") : 
  'x' and 'units' must have length > 0
In addition: There were 23 warnings (use warnings() to see them)
> species_list_pam
[1] "Bb" "Bm" "Bp" "Mn" "Eg" "Ba"
> baleen_sf
Simple feature collection with 240 features and 5 fields
Geometry type: POINT
Dimension:     XY
Bounding box:  xmin: 134988 ymin: 4571754 xmax: 1060413 ymax: 5278400
Projected CRS: WGS 84 / UTM zone 20N
First 10 features:
   site species effort_days detection_days proportion_det                 geometry
1   CBN      Ba         426              5           0.01 POINT (697593.7 5226464)
2   CBN      Bb         426              0           0.00 POINT (697593.7 5226464)
3   CBN      Bm         426            139           0.33 POINT (697593.7 5226464)
4   CBN      Bp         426            186           0.44 POINT (697593.7 5226464)
5   CBN      Eg         426              5           0.01 POINT (697593.7 5226464)
6   CBN      Mn         426             11           0.03 POINT (697593.7 5226464)
7   CCU      Ba         346              2           0.01 POINT (236314.6 4581828)
8   CCU      Bb         346            100           0.29 POINT (236314.6 4581828)
9   CCU      Bm         346             32           0.09 POINT (236314.6 4581828)
10  CCU      Bp         346            221           0.64 POINT (236314.6 4581828)Show more2:23 PMThe error occurs because when you filter baleen_sf by species, some points have the same coordinates (multiple species detected at the same PAM station), and after rjitter() is applied, the number of jittered points doesn't match the number of weight values in the marks.
         The issue is that rjitter() can drop or modify points, but the weights dataframe still has the original number of rows. Let's fix this by adjusting how weights are applied after jittering:Updated KDE Analysis Script with Aerial DataCode ∙ Version 6         # Apply weights first if provided
        if(!is.null(weight_col)) {
          weights <- current_sf[[weight_col]] * 10000
          # Create ppp object with marks
          points_ppp <- ppp(x = species_coords[,"X"], 
                       The fix changes the order of operations:

Before: Create points → jitter → add weights (this caused the mismatch because jitter can change the number of points)
After: Create points with weights already attached → jitter (which preserves the marks) → use the weights

This ensures the weights stay synchronized with the points throughout the jittering process. The PAM data should now process correctly!can we make it so the french names are omitted from the facet plots?2:29 PMUpdated KDE Analysis Script with Aerial DataCode ∙ Version 7   # Filter: not on land -> remove unknown species -> valid cluster size -> species with >50 records
  aerial_sf <- aerial_sf %>%
    filter(!grepl("Unknown|inconnu", Sp_name_clean, ignore.case = TRUE)) %>%
    filter(lengths(st_intersects(., landUTM)Perfect! I've added a line that removes everything after the " /" in the species names. This uses sub() to replace " / French name" with nothing, keeping only the English name.
         For example:
           
           "Humpback Whale / Rorqual a bosse" → "Humpback Whale"
         "Fin Whale / Rorqual commun" → "Fin Whale"
         "Harbour Porpoise / Marsouin commun" → "Harbour Porpoise"
         
         The facet plot titles will now show only the English names!#
           # Unified script for KDE analysis of whale sightings and PAM data
           # Updated 2025
           #
           # ...............
           # DATASET SELECTION - CONFIGURE HERE------
         # ...............
         
         RUN_SIGHTINGS <- T   # Set to TRUE to process WSDB sightings data
         RUN_AERIAL <- TRUE      # Set to TRUE to process aerial sightings datpastedok i made a few changes to the script and i dont see where you added the sub so can you identify it and update to use my variation on the script?2:47 PMYou're right - I don't see the sub() function in your current script. Let me add it back in to remove the French names:Updated KDE Analysis Script with Aerial DataCode ∙ Version 8   # Filter: not on land -> remove unknown species -> remove French names -> valid cluster size -> species with >50 records
         aerial_sf <- aerial_sf %>%
           filter(lengths(st_intersects(., landUTM)) == 0) %>%
           filter(!is.na(Sp_name_clean) & Sp_namePerfect! I've added the sub() function back in. It will remove everything from " /" onwards in the species names, so your facet plot titles will show only the English names (e.g., "Humpback Whale" instead of "Humpback Whale / Rorqual a bosse").where?2:47 PMI apologize for the confusion! Let me update the artifact with your current script version:Updated KDE Analysis Script with Aerial DataCode ∙ Version 9 #
# Unified script for KDE analysis of whale sightings and PAM data
# Updated 2025
#
# ...............
# DATASET SELECTION - CONFIGURE HERE------
# ...............

RUN_SIGHTINGS <- T   # Set to TRUE to process WSDB sightings data
RUN_AERIAL <- TRUE Now I can see the full script! The line that removes French names is on line 424:
rmutate(Sp_name_clean = sub(" /.*$", "", Sp_name_clean)) %>%  # Remove French names after "/"
It's in the aerial sightings section, right after filtering for unknown species and before filtering for valid cluster sizes. This line will strip everything from " /" onwards, leaving only the English species names for the plot titles.can we fix the wsdb sightings species names too? RN they are all caps and weird db format: Species: WHALE-SPERM, WHALE- LONG-FINNED PILOT, DOLPHINS-ATLANTIC BOTTLENOSE, WHALE-FIN, DOLPHINS-ATLANTIC WHITE-SIDED, SEAL-GREY, WHALE-NORTHERN BOTTLENOSE, WHALE-MINKE, WHALE-HUMPBACK, PORPOISE-HARBOUR, WHALE-SEI, WHALE-BLUE, DOLPHINS-COMMON, WHALE-NORTH ATLANTIC RIGHT, DOLPHINS-STRIPED, DOLPHINS-RISSO'S, WHALE-KILLER, WHALE-SOWERBY'S BEAKED 
                  it would be nice to have them match the other species codes in PAM and aerial?2:51 PMGreat idea! Let's clean up the WSDB species names to make them more readable and consistent with the other datasets. I'll convert them to title case and reformat them:Updated KDE Analysis Script with Aerial DataCode ∙ Version 10   # Filter: cetaceans -> not on land -> species with >25 records
                  WS_data_sf <- WS_data_sf %>%
                    filter(SPECIES_CD %in% cetacean_codes) %>%
                    filter(lengths(st_intersects(., landUTM)) == 0) %>%
                    mutate(COMMONNAME = str_to_title(COMMONNAME)) %Actually, that might be too complicated. Let me simplify it to just make them cleaner and more readable:Updated KDE Analysis Script with Aerial DataCode ∙ Version 11   # Filter: cetaceans -> not on land -> clean species names -> species with >25 records
                  WS_data_sf <- WS_data_sf %>%
                    filter(SPECIES_CD %in% cetacean_codes) %>%
                    filter(lengths(st_intersects(., landUTM)) == 0) %>%
                    mutate(COMMONNAME = str_This will transform the names like:
                             
                             "WHALE-SPERM" → "Sperm Whale"
                           "WHALE-HUMPBACK" → "Humpback Whale"
                           "DOLPHINS-ATLANTIC BOTTLENOSE" → "Atlantic Bottlenose Dolphin"
                           "PORPOISE-HARBOUR" → "Harbour Porpoise"
                           "SEAL-GREY" → "Grey Seal"
                           
                           This should match the format used in the aerial data (e.g., "Humpback Whale", "Harbour Porpoise")!Long Finned Pilot Whale has a space before thename, Atlantic Bottlenose Whale"           68       19.7         15          80
 3 "Atlantic White Sided Whale"  should be dolphins, Grey Whale" should be seal, "Harbour Whale"  should be porpoise, 14 "Risso's Whale"       and Striped Whale"             should be dolphin,3:00 PMLet me fix those issues with a better cleaning approach:Updated KDE Analysis Script with Aerial DataCode ∙ Version 12   # Filter: cetaceans -> not on land -> clean species names -> species with >25 records
  WS_data_sf <- WS_data_sf %>%
    filter(SPECIES_CD %in% cetacean_codes) %>%
    filter(lengths(st_intersects(., landUTM)) == 0) %>%
    mutate(COMMONNAME = str_This improved version:

Trims whitespace to fix the leading space issue
Identifies the original animal type from the prefix
Removes all prefixes (Whale, Dolphins, Porpoise, Seal)
Uses special cases to correctly identify dolphins, porpoise, and seal species
Rebuilds the names with the correct suffix

So now you'll get:

"WHALE- LONG-FINNED PILOT" → "Long Finned Pilot Whale"
"DOLPHINS-ATLANTIC BOTTLENOSE" → "Atlantic Bottlenose Dolphin"
"DOLPHINS-ATLANTIC WHITE-SIDED" → "Atlantic White Sided Dolphin"
"SEAL-GREY" → "Grey Seal"
"PORPOISE-HARBOUR" → "Harbour Porpoise"
"DOLPHINS-RISSO'S" → "Risso's Dolphin"
"DOLPHINS-STRIPED" → "Striped Dolphin"
options(scipen = 999)

pacman::p_load(sf, tidyverse, terra, spatstat, lubridate, patchwork, viridis)

#========================
# INPUT / OUTPUT SETUP
#========================

filepath      <- "input/2025/baleen_presence_laura_2025.csv"
raster_dir    <- "output/tif/seasonal/"
plot_dir      <- "outpastedcan you give me a descriptive comment of what this script does that I could paste into the header for future ref>3:05 PMHere's a descriptive header comment for your seasonal KDE script:
r#
# Seasonal KDE Analysis of Baleen Whale PAM Detections
# Updated 2025
#
# This script performs kernel density estimation (KDE) on passive acoustic monitoring (PAM)
# data for baleen whale species, stratified by season (Winter, Spring, Summer, Fall).
#
# Key features:
# - Aggregates daily PAM detections by site, species, and season
# - Calculates proportion of days with detections per season at each station
# - Generates weighted KDE surfaces using detection proportions as weights
# - Creates seasonal comparison plots (2x2 grid) for each species
# - Normalizes density values within each species for visual comparison
# - Uses common spatial window across all species/seasons for consistency
# - Outputs: seasonal raster files (.tif) and combined seasonal plots (.png)
#
# Species analyzed: Sei (Bb), Blue (Bm), Fin (Bp), Humpback (Mn), Minke (Ba), 
#                   North Atlantic Right Whale (Eg)
#
# Seasons defined: Winter (Dec-Feb), Spring (Mar-May), Summer (Jun-Aug), Fall (Sep-Nov)
#
This gives a clear overview of what the script does, its main features, and the key outputs someone would get from running it.# Monthly KDE with Quantile Binning - Matches general map style
# Uses 80%, 90%, 95% quantiles to reduce noise

options(scipen = 999)
pacman::p_load(sf, stars, dplyr, stringr, terra, classInt, ggplot2, 
               patchwork, viridis, ggspatial)

# Directories
raster_input <- "output/tif/monthly/pastedcan you do the same for this script?3:06 PMHere's a descriptive header comment for your monthly quantile KDE script:
  r#
# Monthly KDE Analysis with Quantile Binning for Baleen Whale PAM Detections
# Updated 2025
#
# This script converts monthly kernel density estimation (KDE) rasters into quantile-based
# polygon shapefiles and creates faceted monthly comparison plots for baleen whale species.
#
# Key features:
# - Reads monthly KDE rasters from performMonthlyKDE() output
# - Bins KDE values into quantile classes (80%, 90%, 95%) to reduce noise
# - Converts quantile-binned rasters to vector polygons (shapefiles)
# - Creates 12-month facet plots (4x3 grid) showing acoustic presence by quantile
# - Matches general map style with bathymetric contours, scale bars, and PAM stations
# - Allows filtering of specific quantile levels for cleaner visualization
#
# Processing steps:
# 1. Load monthly KDE rasters (one per species per month)
# 2. Calculate quantile breaks and bin pixel values
# 3. Convert to polygons and save as shapefiles
# 4. Generate faceted plots combining all months for each species
#
# Species analyzed: Sei (Bb), Blue (Bm), Fin (Bp), Humpback (Mn), Minke (Ba), 
#                   North Atlantic Right Whale (Eg)
#
# Inputs:  Monthly KDE rasters (.tif) from output/tif/monthly/
# Outputs: Quantile shapefiles (.shp) in output/shapes/monthly/
#          Monthly facet plots (.png) in output/FIGS/monthly_quantiles/

options(scipen = 999)
pacman::p_load(sf, stars, dplyr, stringr, terra, classInt, ggplot2, 
               patchwork, viridis, ggspatial)

# Directories
raster_input <- "output/tif/monthly/"
shapefile_output <- "output/shapes/monthly/"
plot_output <- "output/FIGS/monthly_quantiles/"

# Create directories if needed
if (!dir.exists(shapefile_output)) dir.create(shapefile_output, recursive = TRUE)
if (!dir.exists(plot_output)) dir.create(plot_output, recursive = TRUE)

# Month names for labeling
month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Species name mappings
species_names <- list(
  "Bb" = "Sei Whale",
  "Bm" = "Blue Whale",
  "Bp" = "Fin Whale",
  "Mn" = "Humpback Whale",
  "Ba" = "Minke Whale",
  "Eg" = "North Atlantic Right Whale"
)

#............................
# STEP 1: Process Monthly KDE Rasters to Quantile Polygons------
#............................

processMonthlyKDE_toQuantiles <- function(tif_file, output_dir = shapefile_output) {
  
  # Load KDE raster
  kde_raster <- rast(tif_file)
  
  # Extract the CRS from the raster
  raster_crs <- crs(kde_raster)  
  
  # Set values <= 0 to NA
  kde_raster[kde_raster <= 0] <- NA
  
  # Define quantile probabilities - 80%, 90%, 95%
  probs <- c(0.80, 0.90, 0.95, 1)
  
  # Calculate quantile breaks (excluding NA values)
  kde_values <- values(kde_raster, mat = FALSE)
  quantiles <- quantile(kde_values, probs = probs, type = 6, na.rm = TRUE)
  
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
  
  # Add quantile labels
  names(kde_sf)[1] <- "quant_level"
  
  # Map the numeric levels back to percentage labels
  kde_sf <- kde_sf %>%
    mutate(Quantil = quant_labels[quant_level])
  
  # Extract species and month from filename
  filename <- basename(tif_file)
  species_code <- str_extract(filename, "(?<=pam_)[A-Za-z]+")
  month_num <- as.numeric(str_extract(filename, "(?<=month)\\d+"))
  
  # Add metadata
  kde_sf <- kde_sf %>%
    mutate(species = species_code,
           month = month_num,
           month_name = month_names[month_num])
  
  # Save the shapefile
  output_path <- file.path(output_dir, 
                           paste0("KDE_pam_", species_code, 
                                  "_month", sprintf("%02d", month_num), ".shp"))
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  cat("Processed:", species_code, "- Month", month_num, "\n")
  
  return(kde_sf)
}

#............................
# STEP 2: Process All Monthly TIF Files------
#............................

cat("\n=== PROCESSING MONTHLY RASTERS TO QUANTILE SHAPEFILES ===\n\n")

# Get all monthly TIF files
monthly_tif_files <- list.files(raster_input, pattern = "KDE_pam_.*_month\\d+\\.tif$", 
                                full.names = TRUE)

cat("Found", length(monthly_tif_files), "monthly TIF files\n\n")

# Process each file with error handling
for (tif_file in monthly_tif_files) {
  tryCatch({
    processMonthlyKDE_toQuantiles(tif_file)
  }, error = function(e) {
    cat("Error processing", basename(tif_file), ":", conditionMessage(e), "\n")
  })
}

#............................
# STEP 3: Create Facet Plots by Species (matching general map style)------
#............................

cat("\n=== CREATING FACET PLOTS ===\n\n")

plotMonthlyQuantileFacets <- function(shapefile_dir, 
                                      land,
                                      contours = NULL,
                                      output_dir = plot_output,
                                      predpal = "mako",
                                      exclude_quantiles = NULL) {
  
  # Load all monthly shapefiles
  shapefiles <- list.files(shapefile_dir, pattern = "KDE_pam_.*_month\\d+\\.shp$", 
                           full.names = TRUE)
  
  if (length(shapefiles) == 0) {
    cat("No monthly shapefiles found in", shapefile_dir, "\n")
    return(NULL)
  }
  
  cat("Reading", length(shapefiles), "shapefiles...\n")
  
  # Read all shapefiles
  all_data <- lapply(shapefiles, function(f) {
    tryCatch({
      shp <- st_read(f, quiet = TRUE)
      
      # Extract metadata from filename if missing
      filename <- basename(f)
      species_code <- str_extract(filename, "(?<=pam_)[A-Za-z]+")
      month_num <- as.numeric(str_extract(filename, "(?<=month)\\d+"))
      
      if (!"species" %in% names(shp)) shp$species <- species_code
      if (!"month" %in% names(shp)) shp$month <- month_num
      if (!"month_name" %in% names(shp)) shp$month_name <- month_names[month_num]
      
      return(shp)
    }, error = function(e) {
      return(NULL)
    })
  })
  
  # Remove NULL entries and combine
  all_data <- all_data[!sapply(all_data, is.null)]
  if (length(all_data) == 0) {
    cat("No valid shapefiles could be read\n")
    return(NULL)
  }
  
  combined_data <- do.call(rbind, all_data)
  
  # Filter out excluded quantiles if specified
  if (!is.null(exclude_quantiles) && length(exclude_quantiles) > 0) {
    combined_data <- combined_data %>%
      filter(!Quantil %in% exclude_quantiles)
  }
  
  # Set up bounding box (matching general map style)
  bbox <- st_bbox(c(xmin = -68, ymin = 41, xmax = -55, ymax = 47.8), 
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
  
  # Convert month to ordered factor
  combined_data$month_name <- factor(combined_data$month_name,
                                     levels = month_names,
                                     ordered = TRUE)
  
  # Get unique species
  unique_species <- sort(unique(combined_data$species))
  
  # Create a plot for each species
  for (sp in unique_species) {
    
    species_data <- combined_data %>% filter(species == sp)
    
    # Get species display name
    species_display <- if (!is.null(species_names[[sp]])) {
      species_names[[sp]]
    } else {
      sp
    }
    
    cat("Creating facet plot for", species_display, "\n")
    
    # Get unique quantiles and sort properly
    unique_quantiles <- unique(species_data$Quantil)
    unique_quantiles <- unique_quantiles[!is.na(unique_quantiles)]
    
    # Extract numeric values for sorting
    quantile_nums <- as.numeric(gsub("%", "", unique_quantiles))
    sort_order <- order(quantile_nums)
    unique_quantiles <- unique_quantiles[sort_order]
    
    # Convert Quantile to ordered factor
    species_data$Quantil <- factor(species_data$Quantil, 
                                   levels = unique_quantiles, 
                                   ordered = TRUE)
    
    # Define palette - match the number of colors to number of quantiles
    n_quantiles <- length(unique_quantiles)
    pal <- rev(get(predpal)(n_quantiles +1))
    
    # Create the faceted plot (matching general map style)
    p <- ggplot() +
      theme_bw()
    
    # Add contours if provided
    if (!is.null(contours)) {
      p <- p + 
        geom_sf(data = contours %>% 
                  dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)), 
                col = "grey70", linewidth = 0.2)+
        geom_sf(data = baleen_stns, color = "grey20", fill = "grey20", shape = 19, alpha = .5) 
    }
    
    # Add KDE quantile polygons
    p <- p +
      geom_sf(data = species_data, aes(fill = Quantil), 
              col = NA, alpha = 0.8, na.rm = TRUE) +
      # Add land
      geom_sf(data = land, color = NA, fill = "grey20") +
      # Facet by month
      facet_wrap(~ month_name, ncol = 4) +
      # Color scale
      scale_fill_manual(values = pal,
                        na.value = NA,
                        na.translate = FALSE,
                        name = "Quantile",
                        labels = unique_quantiles,
                        drop = FALSE) +
      # Coordinate system (matching general map style)
      coord_sf(lims_method = "orthogonal", 
               xlim = xlims, 
               ylim = ylims, 
               crs = st_crs(land), 
               expand = TRUE) +
      # Labels
      labs(title = paste("Monthly Acoustic Presence -", species_display)) +
      ylab("") + 
      xlab("") +
      # Theme settings (matching general map style)
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.background = element_rect(fill = NA),
        legend.key = element_rect(fill = NA),
        plot.margin = margin(5, 5, 5, 5),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        strip.text = element_text(size = 11, face = "bold"),
        strip.background = element_rect(fill = "grey90"),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      ) +
      # Scale bar
      annotation_scale(location = "br", width_hint = 0.25, text_cex = 0.6, 
                       bar_cols = c("grey40", "white"))
    
    # Save the plot
    output_file <- file.path(output_dir, 
                             paste0("monthly_quantile_facets_", sp, ".png"))
    
    ggsave(output_file, p, 
           width = 16, height = 12, dpi = 300, bg = "white", units = "in")
    
    cat("Saved:", output_file, "\n\n")
  }
}

#............................
# STEP 4: Execute Plotting------
#............................

#proj
UTM20 <- terra::crs("+init=epsg:32620") # CODE FOR UTM Zone 20

#load land as UTM
land <- read_sf(here::here("~/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp"))%>%
  dplyr::filter(CONTINENT == "North America")%>%st_transform(UTM20)


# Plot with contours if available (set to NULL if you don't have cont_UTM)
plotMonthlyQuantileFacets(
  shapefile_dir = shapefile_output,
  land = land,
  contours = if(exists("cont_UTM")) cont_UTM else NULL,
  output_dir = plot_output,
  predpal = "mako",
  exclude_quantiles = NULL  # Show all quantiles (80%, 90%, 95%)
)

cat("\n=== MONTHLY QUANTILE ANALYSIS COMPLETE ===\n")
cat("Shapefiles saved to:", shapefile_output, "\n")
cat("Plots saved to:", plot_output, "\n")
