# HEADER --------------------------------------------
#
# Author: Laura Joan Feyrer
# Email:  ljfeyrer@dal.ca
# Most Recent Date Updated: 2026-01-15
#
# Script Name: run_all_kdes_optimized_SIMPLIFIED_v3.R
#
# Description:
## SIMPLIFIED version matching KDES_CE_Collab.R exactly
## Fixed coastline extent to match CE_Collab plotting
#
# Changes from previous version:
# - Fixed coastline cropping to study_area (using st_crop)
# - Fixed plot extent to use study_area bbox (not auto)
# - Matches CE_Collab's coastline handling exactly
#

# SET OPTIONS ---------------------------------------
cat("SETTING OPTIONS... \n\n", sep = "")
options(scipen = 999)
options(encoding = "UTF-8")

# RESOLUTION - MUST MATCH CE_COLLAB EXACTLY
RESOLUTION <- 500  # 500m cells

# ...............
# PERFORMANCE SETTINGS------
# ...............

ENABLE_CACHING <- TRUE
FAST_MODE <- FALSE
USE_PARALLEL <- T
CLEAR_CACHE <- F
CACHE_DIR <- "cache/"

# OUTPUT ORGANIZATION
USE_PARAMETER_FOLDERS <- F
PARAMETER_LABEL <- ""

# ...............
# DATASET AND PARAMETER SELECTION------
# ...............

RUN_SIGHTINGS <- T
RUN_PAM <- T

# ANALYSIS TYPE OPTIONS
RUN_INDIVIDUAL_SPECIES <- T
RUN_FAMILY_GROUPS <- T
RUN_INDIVIDUAL_PAM <- T
RUN_GROUPED_PAM <- T
RUN_BEAKED_PAM <- TRUE
RUN_GROUPED_BEAKED_PAM <- T

# BANDWIDTH OPTIONS
USE_BW_DIGGLE_INDIVIDUAL <- FALSE
USE_BW_DIGGLE_FAMILY <- FALSE

# PROXIMITY WEIGHTING OPTIONS
USE_PROXIMITY_WEIGHTING_INDIVIDUAL <- TRUE
USE_PROXIMITY_WEIGHTING_FAMILY <- FALSE

# BANDWIDTH CONFIGURATION
DEFAULT_SIGHTINGS_SIGMA <- 10000
DEFAULT_PAM_SIGMA <- 10000

SPECIES_BANDWIDTH <- list()

# PROXIMITY WEIGHTING CONFIGURATION
PROXIMITY_K_NEIGHBORS <- 8
PROXIMITY_MAX_DISTANCE <- 5000
PROXIMITY_MIN_WEIGHT <- 0.05

# VISUALIZATION CONFIGURATION
CONTOUR_QUANTILE <- 0.90

# Load necessary libraries
pacman::p_load(sf, tidyverse, readxl, here, leaflet, 
               scales, terra, ggrepel, viridis, ggspatial, spatstat, 
               dplyr, patchwork, digest, lubridate)

if(USE_PARALLEL) {
  pacman::p_load(parallel)
  n_cores <- max(1, detectCores() - 1)
  cat("Parallel processing enabled with", n_cores, "cores\n")
}

# Create parameter-based folder structure
if(USE_PARAMETER_FOLDERS) {
  if(PARAMETER_LABEL != "") {
    param_suffix <- PARAMETER_LABEL
  } else {
    param_suffix <- paste0(
      "bw", DEFAULT_PAM_SIGMA[1],
      "_q", gsub("\\.", "", as.character(CONTOUR_QUANTILE[1]))
    )
  }
  
  raster_dir <- file.path("output/tif", param_suffix,"")
  plot_dir <- file.path("output/FIGS", param_suffix,"")
  shapefile_dir <- file.path("output/shapes", param_suffix,"")
  
  cat("\nUsing parameter folder:", param_suffix, "\n")
} else {
  raster_dir <- "output/tif/"
  plot_dir <- "output/FIGS/"
  shapefile_dir <- "output/shapes/"
}

# Projections
UTM20 <- crs("+init=epsg:32620")

# Create cache directory if needed
if(ENABLE_CACHING && !dir.exists(CACHE_DIR)) {
  dir.create(CACHE_DIR, recursive = TRUE)
}

# ................................
# CACHING FUNCTIONS------
# ................................

save_to_cache <- function(obj, name, subdir = NULL) {
  if(!ENABLE_CACHING) return(invisible(NULL))
  
  cache_path <- CACHE_DIR
  if(!is.null(subdir)) {
    cache_path <- file.path(CACHE_DIR, subdir)
    if(!dir.exists(cache_path)) dir.create(cache_path, recursive = TRUE)
  }
  
  cache_file <- file.path(cache_path, paste0(name, ".rds"))
  saveRDS(obj, cache_file)
  cat("  Cached:", cache_file, "\n")
}

load_from_cache <- function(name, subdir = NULL) {
  if(!ENABLE_CACHING) return(NULL)
  
  cache_path <- CACHE_DIR
  if(!is.null(subdir)) {
    cache_path <- file.path(CACHE_DIR, subdir)
  }
  
  cache_file <- file.path(cache_path, paste0(name, ".rds"))
  if(file.exists(cache_file)) {
    cat("  Loading from cache:", cache_file, "\n")
    return(readRDS(cache_file))
  }
  return(NULL)
}

get_cache_hash <- function(...) {
  params <- list(...)
  digest::digest(params, algo = "md5")
}

# ................................
# CLEANUP------
# ................................

if(!FAST_MODE) {
  cat("\n=== CLEANING OUTPUT FOLDERS ===\n\n")
  
  prefixes_to_clean <- c()
  
  if(RUN_SIGHTINGS) {
    if(RUN_INDIVIDUAL_SPECIES) {
      prefixes_to_clean <- c(prefixes_to_clean, "sightings_species")
    }
    if(RUN_FAMILY_GROUPS) {
      prefixes_to_clean <- c(prefixes_to_clean, "sightings_odontocetes", "sightings_mysticetes")
    }
  }
  
  if(RUN_PAM) {
    if(RUN_INDIVIDUAL_PAM) {
      prefixes_to_clean <- c(prefixes_to_clean, "pam_baleen_")
    }
    if(RUN_GROUPED_PAM) {
      prefixes_to_clean <- c(prefixes_to_clean, "grouped_baleen_pam")
    }
  }
  
  if(RUN_BEAKED_PAM) {
    prefixes_to_clean <- c(prefixes_to_clean, "pam_beaked_")
    if(RUN_GROUPED_BEAKED_PAM) {
      prefixes_to_clean <- c(prefixes_to_clean, "grouped_beaked_pam")
    }
  }
  
  if (!ENABLE_CACHING || CLEAR_CACHE) {
    if (dir.exists(raster_dir[1])) {
      tif_files <- list.files(
        raster_dir[1],
        pattern = "\\.tif$",
        full.names = TRUE,
        recursive = TRUE
      )
      
      files_to_remove <- character(0)
      
      if (length(tif_files) > 0 && length(prefixes_to_clean) > 0) {
        rm_idx <- vapply(tif_files, function(f) {
          any(vapply(prefixes_to_clean, function(p) grepl(p, basename(f)), logical(1)))
        }, logical(1))
        
        files_to_remove <- tif_files[rm_idx]
      }
      
      if (length(files_to_remove) > 0) {
        file.remove(files_to_remove)
        cat("Removed", length(files_to_remove), "raster files for enabled analyses\n")
      } else {
        cat("No matching raster files to remove in", raster_dir[1], "\n")
      }
      
    } else {
      dir.create(raster_dir[1], recursive = TRUE)
    }
    
  } else {
    cat("Keeping cached raster files\n")
    if (!dir.exists(raster_dir[1])) {
      dir.create(raster_dir[1], recursive = TRUE)
    }
  }
  
  if (dir.exists(shapefile_dir[1])) {
    shp_files <- list.files(shapefile_dir[1], pattern = "\\.(shp|shx|dbf|prj|cpg)$", 
                            full.names = TRUE, recursive = TRUE)
    
    if(length(shp_files) > 0 && length(prefixes_to_clean) > 0) {
      # Use vapply to ensure logical vector output
      should_remove <- vapply(shp_files, function(f) {
        any(vapply(prefixes_to_clean, function(p) grepl(p, basename(f)), logical(1)))
      }, logical(1))
      
      files_to_remove <- shp_files[should_remove]
    } else {
      files_to_remove <- character(0)
    }
    
    if (length(files_to_remove) > 0) {
      file.remove(files_to_remove)
      cat("Removed", length(files_to_remove), "shapefile components for enabled analyses\n")
    }
  } else {
    dir.create(shapefile_dir[1], recursive = TRUE)
  }

  if (dir.exists(plot_dir[1])) {
    png_files <- list.files(plot_dir[1], pattern = "\\.png$", full.names = TRUE)
    
    plot_prefixes <- c()
    if(RUN_SIGHTINGS) {
      if(RUN_INDIVIDUAL_SPECIES) plot_prefixes <- c(plot_prefixes, "sightings_species")
      if(RUN_FAMILY_GROUPS) plot_prefixes <- c(plot_prefixes, "sightings_odontocetes", "sightings_mysticetes")
    }
    if(RUN_PAM) {
      if(RUN_INDIVIDUAL_PAM) plot_prefixes <- c(plot_prefixes, "pam_baleen_")
      if(RUN_GROUPED_PAM) plot_prefixes <- c(plot_prefixes, "grouped_baleen_pam")
    }
    if(RUN_BEAKED_PAM) {
      plot_prefixes <- c(plot_prefixes, "pam_beaked_")
      if(RUN_GROUPED_BEAKED_PAM) plot_prefixes <- c(plot_prefixes, "grouped_beaked_pam")
    }
    
    # Use vapply instead of sapply to ensure logical vector
    if(length(png_files) > 0 && length(plot_prefixes) > 0) {
      should_remove <- vapply(png_files, function(f) {
        any(vapply(plot_prefixes, function(p) grepl(p, basename(f)), logical(1)))
      }, logical(1))
      
      files_to_remove <- png_files[should_remove]
    } else {
      files_to_remove <- character(0)
    }
    
    if (length(files_to_remove) > 0) {
      file.remove(files_to_remove)
      cat("Removed", length(files_to_remove), "plot files for enabled analyses\n")
    }
  } else {
    dir.create(plot_dir[1], recursive = TRUE)
  }

if(CLEAR_CACHE && dir.exists(CACHE_DIR)) {
  cache_files <- list.files(CACHE_DIR, full.names = TRUE, recursive = TRUE)
  if(length(cache_files) > 0) {
    file.remove(cache_files)
    cat("Cleared", length(cache_files), "cached files\n")
  }
}
}
# ................................
# CALCULATE PROXIMITY WEIGHTS------
# ................................

calculate_proximity_weights <- function(data_sf, 
                                        k_neighbors = PROXIMITY_K_NEIGHBORS, 
                                        max_distance = PROXIMITY_MAX_DISTANCE,
                                        min_weight = PROXIMITY_MIN_WEIGHT) {
  
  cat("  Calculating proximity weights with:\n")
  cat("    - K neighbors:", k_neighbors, "\n")
  cat("    - Max distance:", max_distance, "m\n")
  cat("    - Min weight:", min_weight, "\n")
  
  coords <- st_coordinates(data_sf)
  dist_matrix <- as.matrix(dist(coords))
  
  proximity_scores <- apply(dist_matrix, 1, function(row) {
    sorted_dists <- sort(row)
    nearest_k <- sorted_dists[2:min(k_neighbors + 1, length(sorted_dists))]
    mean(nearest_k)
  })
  
  weights <- 1 / (1 + (proximity_scores / max_distance))
  weights <- (weights - min(weights)) / (max(weights) - min(weights))
  weights <- pmax(weights, min_weight)
  
  cat("  Weight statistics:\n")
  cat("    - Min:", round(min(weights), 3), "\n")
  cat("    - Mean:", round(mean(weights), 3), "\n")
  cat("    - Max:", round(max(weights), 3), "\n")
  cat("    - Points below 0.5 weight:", sum(weights < 0.5), 
      "(", round(100 * sum(weights < 0.5) / length(weights), 1), "%)\n")
  
  return(weights)
}

# ................................
# SIMPLIFIED CONTOUR CREATION------
# ................................

create_contour_90_baseline <- function(kde_raster, output_path) {
  kde_raster[kde_raster <= 0] <- NA
  
  kde_values <- values(kde_raster, mat = FALSE, na.rm = TRUE)
  threshold_90 <- quantile(kde_values, probs = 0.90, type = 6, na.rm = TRUE)
  
  cat("    Threshold value:", threshold_90, "\n")
  
  kde_binary <- kde_raster
  kde_binary[kde_binary < threshold_90] <- NA
  kde_binary[!is.na(kde_binary)] <- 1
  
  kde_poly <- as.polygons(kde_binary, dissolve = TRUE)
  kde_sf <- st_as_sf(kde_poly)
  kde_sf <- st_make_valid(kde_sf)
  
  cat("    Created", nrow(kde_sf), "polygon(s)\n")
  
  kde_sf <- kde_sf %>%
    mutate(contour = "0.90") %>%
    select(contour, geometry)
  
  write_sf(kde_sf, output_path, delete_dsn = TRUE)
  
  cat("    Saved to:", output_path, "\n")
  
  return(kde_sf)
}

# ................................
# UNIFIED KDE FUNCTION------
# ................................

performKDE <- function(data_sf, 
                       species_col,
                       species_list, 
                       weight_col = NULL,
                       buffer_percent = 0.05,
                       sigma_val = 10000,
                       window = NULL,
                       output_prefix,
                       threshold_quantile = NULL,
                       output_dir = raster_dir,
                       plot_output_dir = plot_dir,
                       shapefile_output_dir = shapefile_dir,
                       use_bw_diggle = FALSE,
                       data_source_label = "KDE",
                       resolution = RESOLUTION,
                       coastline_for_plot = NULL,
                       owa_for_plot = NULL,
                       study_area_for_plot = NULL) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  if (!dir.exists(plot_output_dir)) dir.create(plot_output_dir, recursive = TRUE)
  
  shapefile_crs <- st_crs(data_sf)$wkt
  
  raster_list <- list()
  sigma_store <- list() 
  global_min <- Inf
  global_max <- -Inf
  
  sanitize_filename <- function(name) {
    name <- trimws(name)
    name <- gsub("[/\\:*?\"<>|]", "_", name)
    name <- gsub(" ", "_", name)
    return(name)
  }
  
  cache_params <- list(
    species_list = species_list,
    weight_col = weight_col,
    sigma_val = sigma_val,
    use_bw_diggle = use_bw_diggle,
    output_prefix = output_prefix,
    resolution = resolution,
    buffer_percent = buffer_percent,
    proximity_settings = if(!is.null(weight_col)) {
      list(PROXIMITY_K_NEIGHBORS, PROXIMITY_MAX_DISTANCE, PROXIMITY_MIN_WEIGHT)
    } else NULL
  )
  cache_key <- get_cache_hash(cache_params)
  cache_name <- paste0(output_prefix, "_", cache_key)
  
  #.............................
  # FIRST PASS: Create or Load KDE rasters-----
  #.............................
  
  if(FAST_MODE) {
    cat("\n=== FAST MODE: Loading existing rasters ===\n")
    for(species in species_list) {
      species_clean <- sanitize_filename(species)
      bw_label <- ifelse(use_bw_diggle, "diggle_bw*", paste0("fixed_bw", round(sigma_val, 0)))
      pattern <- paste0("KDE_", output_prefix, "_", species_clean, "_", bw_label, ".tif")
      tif_files <- list.files(output_dir, pattern = glob2rx(pattern), full.names = TRUE)
      
      if(length(tif_files) > 0) {
        raster_kd <- rast(tif_files[1])
        bw_match <- regmatches(basename(tif_files[1]), regexpr("bw[0-9]+", basename(tif_files[1])))
        sigma_used <- as.numeric(gsub("bw", "", bw_match))
        bw_method <- ifelse(grepl("diggle", basename(tif_files[1])), "diggle", "fixed")
        bw_label_actual <- paste0(bw_method, "_bw", sigma_used)
        
        raster_list[[species]] <- list(
          raster = raster_kd, 
          sigma = sigma_used, 
          method = bw_method,
          bw_label = bw_label_actual,
          raster_path = tif_files[1]
        )
        
        local_min <- min(values(raster_kd), na.rm = TRUE)
        local_max <- max(values(raster_kd), na.rm = TRUE)
        global_min <- min(global_min, local_min)
        global_max <- max(global_max, local_max)
        
        cat("  Loaded raster for", species, "\n")
      } else {
        cat("  Warning: No raster found for", species, "- will need to calculate\n")
      }
    }
  } else {
    cached_metadata <- load_from_cache(cache_name, "rasters")
    
    if(!is.null(cached_metadata)) {
      cat("\n=== Using cached KDE results ===\n")
      
      raster_list <- list()
      for(species in names(cached_metadata$raster_list)) {
        raster_path <- cached_metadata$raster_list[[species]]$raster_path
        if(file.exists(raster_path)) {
          raster_kd <- rast(raster_path)
          
          raster_list[[species]] <- cached_metadata$raster_list[[species]]
          raster_list[[species]]$raster <- raster_kd
          
          cat("  Loaded", species, "from", basename(raster_path), "\n")
        } else {
          cat("  Warning: Cached raster not found for", species, "- will recalculate\n")
          cached_metadata <- NULL
          break
        }
      }
      
      if(!is.null(cached_metadata)) {
        sigma_store <- cached_metadata$sigma_store
        global_min <- cached_metadata$global_min
        global_max <- cached_metadata$global_max
      }
    }
    
    if(is.null(cached_metadata)) {
      cat("\n=== Calculating KDE rasters ===\n")
      
      for(species in species_list) {
        current_sf <- data_sf[data_sf[[species_col]] == species, ]
        
        if(nrow(current_sf) > 0) {
          cat("Processing", species, "(", nrow(current_sf), "records)...\n")
          
          raster_kd <- NULL
          sigma_used <- NULL
          bw_method <- NULL
          
          tryCatch({
            species_coords <- st_coordinates(current_sf)
            
            if(is.null(window)) {
              x_range <- range(species_coords[, "X"])
              y_range <- range(species_coords[, "Y"])
              x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
              y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
              expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
              expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
              species_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
            } else {
              species_window <- window
            }
            
            if(!is.null(weight_col)) {
              weights <- current_sf[[weight_col]] * 10000
              points_ppp <- ppp(x = species_coords[,"X"], 
                                y = species_coords[,"Y"], 
                                window = species_window,
                                marks = data.frame(weights = weights))
            } else {
              points_ppp <- ppp(x = species_coords[,"X"], 
                                y = species_coords[,"Y"], 
                                window = species_window)
            }
            
            points_ppp <- rjitter(points_ppp, retry = TRUE, nsim = 1, drop = TRUE)
            
            if(!is.null(weight_col)) {
              weight_arg <- marks(points_ppp)
            } else {
              weight_arg <- NULL
            }
            
            if (isTRUE(use_bw_diggle)) {
              sigma_used <- bw.diggle(points_ppp)
              bw_method <- "diggle"
            } else {
              if(exists("SPECIES_BANDWIDTH") && species %in% names(SPECIES_BANDWIDTH)) {
                sigma_used <- SPECIES_BANDWIDTH[[species]]
                bw_method <- "fixed"
                cat("  Using species-specific bandwidth:", round(sigma_used, 0), "m\n")
              } else {
                sigma_used <- sigma_val
                bw_method <- "fixed"
              }
            }
            
            sigma_store[[species]] <- list(sigma = sigma_used, method = bw_method)
            
            # Calculate dimensions for specified resolution
            x_extent <- diff(species_window$xrange)
            y_extent <- diff(species_window$yrange)
            dimx <- round(x_extent / resolution)
            dimy <- round(y_extent / resolution)
            
            cat("  Resolution:", resolution, "m (grid dimensions:", dimx, "x", dimy, ")\n")
            
            # MATCH CE_COLLAB: Remove adaptive parameter
            kd_result <- density.ppp(points_ppp, 
                                     sigma = sigma_used, 
                                     positive = TRUE,
                                     kernel = "gaussian", 
                                     weights = weight_arg,
                                     dimyx = c(dimy, dimx),
                                     diggle = TRUE)
            
            raster_kd <- rast(kd_result)
            crs(raster_kd) <- shapefile_crs
            
            # Normalize by species max
            max_val <- max(values(raster_kd), na.rm = TRUE)
            if (!is.na(max_val) && max_val > 0) {
              raster_kd <- raster_kd / max_val
              
              if(!is.null(threshold_quantile)) {
                kde_vals <- values(raster_kd, mat = FALSE, na.rm = TRUE)
                threshold <- quantile(kde_vals, threshold_quantile, na.rm = TRUE)
                raster_kd[raster_kd < threshold] <- NA
              }
            }
            
            # Report actual resolution
            actual_res <- res(raster_kd)
            cat("  Actual raster resolution:", round(mean(actual_res), 2), "m\n")
            
            local_min <- min(values(raster_kd), na.rm = TRUE)
            local_max <- max(values(raster_kd), na.rm = TRUE)
            global_min <- min(global_min, local_min)
            global_max <- max(global_max, local_max)
            
            bw_label <- paste0(bw_method, "_bw", round(sigma_used, 0))
            raster_filename <- paste0(output_dir, "KDE_", output_prefix, "_", 
                                      sanitize_filename(species), "_", bw_label, ".tif")
            writeRaster(raster_kd, filename = raster_filename, overwrite = TRUE)
            
            raster_list[[species]] <- list(raster = raster_kd, 
                                           sigma = sigma_used, 
                                           method = bw_method,
                                           bw_label = bw_label,
                                           raster_path = raster_filename)
            
          }, error = function(e) {
            cat("  Error processing species", species, ":", e$message, "\n")
          })
        }
      }
      
      cache_obj <- list(
        raster_list = lapply(raster_list, function(x) {
          list(
            sigma = x$sigma,
            method = x$method,
            bw_label = x$bw_label,
            raster_path = x$raster_path
          )
        }),
        sigma_store = sigma_store,
        global_min = global_min,
        global_max = global_max
      )
      save_to_cache(cache_obj, cache_name, "rasters")
    }
  }
  
  #.............................
  # GENERATE QUANTILE SHAPEFILES------
  #.............................
  
  cat("\n=== GENERATING QUANTILE SHAPEFILES (BASELINE METHOD) ===\n")
  
  if (!dir.exists(shapefile_output_dir)) {
    dir.create(shapefile_output_dir, recursive = TRUE)
  }
  
  q90_dir <- file.path(shapefile_output_dir, "q90")
  if(!dir.exists(q90_dir)) dir.create(q90_dir, recursive = TRUE)
  
  shapefile_paths <- list()
  
  for(species in names(raster_list)) {
    tryCatch({
      tif_file <- raster_list[[species]]$raster_path
      kde_raster <- rast(tif_file)
      
      base_name <- tools::file_path_sans_ext(basename(tif_file))
      species_clean <- str_extract(base_name, "(?<=KDE_)[^_]+_.*")
      
      contour_filename <- paste0(species_clean, "_contour90.shp")
      contour_path <- file.path(q90_dir, contour_filename)
      
      cat("  Creating contour for", species, "\n")
      contour_sf <- create_contour_90_baseline(kde_raster, contour_path)
      
      shapefile_paths[[species]] <- contour_path
      raster_list[[species]]$shapefile_path <- contour_path
      raster_list[[species]]$contour_sf <- contour_sf
      
    }, error = function(e) {
      cat("  Error:", species, "-", conditionMessage(e), "\n")
    })
  }
  
  cat("\n✓ Created", length(shapefile_paths), "baseline contour shapefiles\n\n")
  
  #.............................
  # CREATE PLOTS------
  #.............................
  
  contour_dir <- "q90"
  
  # Create readable species name mapping
  readable_names <- list()
  for(sp in species_list) {
    if(nchar(sp) <= 6 && !grepl(" ", sp)) {
      if(sp == "Bb") readable_names[[sp]] <- "Sei Whale"
      else if(sp == "Bm") readable_names[[sp]] <- "Blue Whale"
      else if(sp == "Bp") readable_names[[sp]] <- "Fin Whale"
      else if(sp == "Mn") readable_names[[sp]] <- "Humpback Whale"
      else if(sp == "Ba") readable_names[[sp]] <- "Minke Whale"
      else if(sp == "Eg") readable_names[[sp]] <- "North Atlantic Right Whale"
      else if(sp == "Ha") readable_names[[sp]] <- "Northern Bottlenose Whale"
      else if(sp == "Mb") readable_names[[sp]] <- "Sowerby's Beaked Whale"
      else if(sp == "Zc") readable_names[[sp]] <- "Goose Beaked Whale"
      else if(sp == "MmMe") readable_names[[sp]] <- "True's/Gervais' Beaked Whale"
      else readable_names[[sp]] <- sp
    } else {
      readable_names[[sp]] <- sp
    }
  }
  
  # CRITICAL FIX: Use study_area bbox for consistent plotting extent
  if(!is.null(study_area_for_plot)) {
    study_bbox <- st_bbox(study_area_for_plot)
    xlims <- c(study_bbox["xmin"], study_bbox["xmax"])
    ylims <- c(study_bbox["ymin"], study_bbox["ymax"])
  } else {
    xlims <- NULL
    ylims <- NULL
  }
  
  plot_list <- list()
  
  cat("=== Generating plots ===\n")
  for(species in species_list) {
    raster_data <- raster_list[[species]]
    if(!is.null(raster_data)) {
      cat("  Plotting", species, "...\n")
      
      raster_kd <- raster_data$raster
      sigma_used <- raster_data$sigma
      bw_method <- raster_data$method
      bw_label <- raster_data$bw_label
      
      species_display <- readable_names[[species]]
      
      kde_df <- as.data.frame(raster_kd, xy = TRUE) %>%
        filter(!is.na(lyr.1)) %>%
        mutate(lyr.1_sqrt = sqrt(lyr.1))
      
      contour_sf <- raster_data$contour_sf
      
      if(!is.null(contour_sf)) {
        cat("    Contour available:", nrow(contour_sf), "features\n")
      }
      
      species_points <- data_sf[data_sf[[species_col]] == species, ]
      is_pam <- "site" %in% names(species_points) || "deployment" %in% names(species_points)
      
      bw_text <- paste0(toupper(bw_method), " BW: ", round(sigma_used, 0), "m")
      contour_text <- paste0(CONTOUR_QUANTILE * 100, "% contour (baseline)")
      plot_title <- paste0(data_source_label, ": ", species_display)
      
      p <- ggplot()
      
      # Add KDE tiles first
      p <- p + geom_tile(data = kde_df, aes(x = x, y = y, fill = lyr.1_sqrt)) +
        scale_fill_viridis_c(option = "viridis", 
                             name = "sqrt(Density)",
                             limits = c(0, max(kde_df$lyr.1_sqrt)),
                             begin = 0.15,
                             end = 1.0)
      
      # Add contour
      if(!is.null(contour_sf) && nrow(contour_sf) > 0) {
        p <- p + geom_sf(data = contour_sf, fill = NA, color = "white", 
                         linewidth = .8, linetype = "solid")
      }
      
      # Add coastline (CROPPED to study area)
      if(!is.null(coastline_for_plot)) {
        p <- p + geom_sf(data = coastline_for_plot, fill = "grey90", 
                         color = "grey70", linewidth = 0.3)
      }
      
      # Add OWA
      if(!is.null(owa_for_plot)) {
        p <- p + geom_sf(data = owa_for_plot, fill = NA, color = "red", 
                         linewidth = 0.5, linetype = "dashed")
      }
      
      # Add data points
      if(nrow(species_points) > 0) {
        if(is_pam) {
          p <- p + geom_sf(data = species_points, color = "white", size = 1.5, 
                           shape = 21, fill = "black", stroke = 0.5, alpha = 0.8)
        } else {
          p <- p + geom_sf(data = species_points, color = "white", size = 0.8, 
                           shape = 21, fill = "black", stroke = 0.3, alpha = 0.6)
        }
      }
      
      # Add study area outline (top layer)
      if(!is.null(study_area_for_plot)) {
        p <- p + geom_sf(data = study_area_for_plot, fill = NA, color = "black", linewidth = 0.8)
      }
      
      # CRITICAL FIX: Use explicit xlim/ylim to match study area
      if(!is.null(xlims) && !is.null(ylims)) {
        p <- p + coord_sf(crs = st_crs(UTM20), xlim = xlims, ylim = ylims, expand = FALSE)
      } else {
        p <- p + coord_sf(crs = st_crs(UTM20))
      }
      
      p <- p + 
        theme_minimal() +
        labs(title = plot_title) +
        theme(
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          legend.position = "right",
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          panel.grid = element_line(color = "grey95")
        )
      
      # Add annotations
      if(!is.null(xlims) && !is.null(ylims)) {
        p <- p +
          annotate("text", x = xlims[2], y = ylims[2], 
                   label = bw_text,
                   hjust = 1.1, vjust = 1.5, size = 3, color = "orange") +
          annotate("text", x = xlims[2], y = ylims[2], 
                   label = contour_text,
                   hjust = 1.1, vjust = 3.0, size = 3, color = "white")
      }
      
      plot_list[[species]] <- p
    }
  }
  
  # Save plots
  if(length(plot_list) > 0) {
    
    is_grouped <- grepl("grouped|combined_mysticetes|combined_odontocetes|guild", 
                        output_prefix, ignore.case = TRUE)
    
    if (is_grouped && length(plot_list) > 1) {
      plot_combined <- wrap_plots(plot_list, guides = "collect")
      combined_filename <- paste0(plot_output_dir, "combined_KDE_", output_prefix, "_plots.png")
      ggsave(filename = combined_filename, plot = plot_combined, width = 16, height = 12, dpi = 300)
      cat("  Saved combined plot:", combined_filename, "\n")
    } else {
      plot_combined <- NULL
    }
    
    for(species in names(plot_list)) {
      if(!is.null(plot_list[[species]])) {
        bw_label <- raster_list[[species]]$bw_label
        individual_filename <- paste0(plot_output_dir, "KDE_", output_prefix, "_", 
                                      sanitize_filename(species), "_", bw_label, ".png")
        ggsave(filename = individual_filename, plot = plot_list[[species]], width = 8, height = 6, dpi = 300)
      }
    }
  } else {
    plot_combined <- NULL
  }
  
  return(list(
    plot = plot_combined,
    bandwidths = sigma_store,
    global_range = c(global_min, global_max)
  ))
}

#.............................
# LOAD STUDY AREA (WITH CACHING)-------
#.............................

cat("\n=== LOADING STUDY AREA ===\n\n")

study_area <- load_from_cache("study_area", "spatial")
common_window <- load_from_cache("common_window", "spatial")

if(is.null(study_area)) {
  study_area_path <- "shapefiles/studyArea/ESS_study_area_simple.shp"
  study_area <- st_read(study_area_path, quiet = TRUE)
  
  if(st_crs(study_area)$epsg != 32620) {
    study_area <- st_transform(study_area, UTM20)
  }
  
  save_to_cache(study_area, "study_area", "spatial")
}

if(is.null(common_window)) {
  study_bbox <- st_bbox(study_area)
  buffer_percent <- 0.05
  x_range <- c(study_bbox["xmin"], study_bbox["xmax"])
  y_range <- c(study_bbox["ymin"], study_bbox["ymax"])
  x_buffer <- (x_range[2] - x_range[1]) * buffer_percent
  y_buffer <- (y_range[2] - y_range[1]) * buffer_percent
  expanded_x_range <- c(x_range[1] - x_buffer, x_range[2] + x_buffer)
  expanded_y_range <- c(y_range[1] - y_buffer, y_range[2] + y_buffer)
  common_window <- owin(xrange = expanded_x_range, yrange = expanded_y_range)
  
  save_to_cache(common_window, "common_window", "spatial")
  
  cat("Common window created:\n")
  cat("  X range:", common_window$xrange, "\n")
  cat("  Y range:", common_window$yrange, "\n")
}

cat("Study area loaded\n")
cat("Window buffer:", buffer_percent * 100, "%\n\n")

#.............................
# LOAD BASE LAYERS FOR PLOTTING (WITH CACHING)-------
#.............................

cat("=== LOADING BASE LAYERS FOR PLOTTING ===\n\n")

# CRITICAL FIX: Crop coastline to study area (matching CE_Collab Line 289-295)
coastline <- load_from_cache("coastline_cropped", "spatial")
if(is.null(coastline)) {
  coastline_path <- "shapefiles/Canada/NE_10mLand.shp"
  if(file.exists(coastline_path)) {
    coastline_full <- st_read(coastline_path, quiet = TRUE) %>%
      st_transform(UTM20)
    
    # MATCH CE_COLLAB: Use st_crop to clip to study area bbox
    coastline <- st_crop(coastline_full, st_bbox(study_area))
    
    save_to_cache(coastline, "coastline_cropped", "spatial")
    cat("Loaded and cropped coastline to study area\n")
  } else {
    coastline <- NULL
    cat("Coastline file not found\n")
  }
} else {
  cat("Loaded coastline from cache (already cropped)\n")
}

# Load OWA and crop to study area (matching CE_Collab Line 298-304)
owa <- load_from_cache("owa_cropped", "spatial")
if(is.null(owa)) {
  owa_path <- "shapefiles/WEA/Designated_WEAs_25_07_29.shp"
  if(file.exists(owa_path)) {
    owa_full <- st_read(owa_path, quiet = TRUE) %>%
      st_transform(UTM20)
    
    # MATCH CE_COLLAB: Use st_crop to clip to study area bbox
    owa <- st_crop(owa_full, st_bbox(study_area))
    
    save_to_cache(owa, "owa_cropped", "spatial")
    cat("Loaded and cropped OWA to study area\n")
  } else {
    owa <- NULL
    cat("OWA file not found\n")
  }
} else {
  cat("Loaded OWA from cache (already cropped)\n")
}

cat("\n")

#.............................
# PROCESS COMBINED SIGHTINGS DATA (WITH CACHING)----
#.............................

if(RUN_SIGHTINGS) {
  
  cat("\n=== PROCESSING COMBINED SIGHTINGS DATA ===\n\n")
  
  cetacean_sf <- load_from_cache("cetacean_sf_2015plus", "sightings")
  species_25plus <- load_from_cache("species_25plus", "sightings")
  
  if(is.null(cetacean_sf) || is.null(species_25plus)) {
    combined_data <- read_csv("input/2025/combined_sights/combined_dedup_1km_day_clipped_classified.csv",
                              show_col_types = FALSE)
    
    cat("Total records loaded:", nrow(combined_data), "\n")
    
    combined_data <- combined_data %>%
      filter(as.Date(date_utc) >= as.Date("2015-01-01"))
    
    cat("Records after 2015+ filter:", nrow(combined_data), "\n")
    
    combined_sf <- st_as_sf(combined_data, 
                            coords = c("lon", "lat"), 
                            crs = 4326) %>%
      st_transform(UTM20) %>%
      st_make_valid()
    
    cetacean_sf <- combined_sf %>%
      filter(cetacean_family %in% c("odontocete", "mysticete"))
    
    cat("Cetacean records:", nrow(cetacean_sf), "\n")
    
    species_counts <- cetacean_sf %>%
      mutate(common_name = case_when(
        str_detect(common_name, "SEI") ~ "WHALE-SEI",
        TRUE ~ common_name
      )) %>%
      st_drop_geometry() %>%
      group_by(common_name, cetacean_family) %>%
      summarise(n_records = n(), .groups = "drop") %>%
      arrange(desc(n_records))
    
    species_25plus <- species_counts %>%
      filter(n_records >= 25) %>%
      pull(common_name)
    
    save_to_cache(cetacean_sf, "cetacean_sf_2015plus", "sightings")
    save_to_cache(species_25plus, "species_25plus", "sightings")
    
    cat("\nCached processed sightings data\n")
  } else {
    cat("Loaded sightings data from cache\n")
  }
  
  cat("\n", length(species_25plus), "species with ≥25 records\n")
  
  #.............................
  # INDIVIDUAL SPECIES KDEs----
  #.............................
  
  if(RUN_INDIVIDUAL_SPECIES && length(species_25plus) > 0) {
    
    if(USE_PROXIMITY_WEIGHTING_INDIVIDUAL) {
      cat("\n=== INDIVIDUAL SPECIES KDEs WITH PROXIMITY WEIGHTING ===\n\n")
      
      species_sf <- cetacean_sf %>%
        filter(common_name %in% species_25plus)
      
      weights_cache_key <- get_cache_hash(list(
        nrow(species_sf),
        PROXIMITY_K_NEIGHBORS,
        PROXIMITY_MAX_DISTANCE,
        PROXIMITY_MIN_WEIGHT
      ))
      weights_cache_name <- paste0("proximity_weights_", weights_cache_key)
      
      cached_weights <- load_from_cache(weights_cache_name, "weights")
      
      if(!is.null(cached_weights)) {
        species_sf$proximity_weight <- cached_weights
        cat("Loaded proximity weights from cache\n")
      } else {
        cat("Calculating proximity-based weights...\n")
        species_sf$proximity_weight <- calculate_proximity_weights(species_sf)
        save_to_cache(species_sf$proximity_weight, weights_cache_name, "weights")
      }
      
      res_species <- performKDE(
        data_sf = species_sf, 
        species_col = "common_name",
        species_list = species_25plus, 
        weight_col = "proximity_weight",
        window = common_window,
        output_prefix = "sightings_species_proximity",
        sigma_val = DEFAULT_SIGHTINGS_SIGMA,
        use_bw_diggle = USE_BW_DIGGLE_INDIVIDUAL,
        data_source_label = "Sightings (Proximity Weighted)",
        resolution = RESOLUTION,
        buffer_percent = 0.05,
        coastline_for_plot = coastline,
        owa_for_plot = owa,
        study_area_for_plot = study_area
      )
      
    } else {
      cat("\n=== INDIVIDUAL SPECIES KDEs WITHOUT PROXIMITY WEIGHTING ===\n\n")
      
      species_sf <- cetacean_sf %>%
        filter(common_name %in% species_25plus)
      
      res_species <- performKDE(
        data_sf = species_sf, 
        species_col = "common_name",
        species_list = species_25plus, 
        weight_col = NULL,
        window = common_window,
        output_prefix = "sightings_species",
        sigma_val = DEFAULT_SIGHTINGS_SIGMA,
        use_bw_diggle = USE_BW_DIGGLE_INDIVIDUAL,
        data_source_label = "Sightings (Unweighted)",
        resolution = RESOLUTION,
        buffer_percent = 0.05,
        coastline_for_plot = coastline,
        owa_for_plot = owa,
        study_area_for_plot = study_area
      )
    }
  }
  
  #.............................
  # FAMILY-LEVEL KDEs----
  #.............................
  
  if(RUN_FAMILY_GROUPS) {
    
    cat("\n=== FAMILY-LEVEL KDEs ===\n\n")
    
    odontocete_sf <- cetacean_sf %>%
      filter(cetacean_family == "odontocete") %>%
      mutate(family_group = "Odontocetes (All Toothed Whales)")
    
    if(nrow(odontocete_sf) > 0) {
      if(USE_PROXIMITY_WEIGHTING_FAMILY) {
        weights_cache_key <- paste0("odonto_weights_", get_cache_hash(list(nrow(odontocete_sf))))
        cached_weights <- load_from_cache(weights_cache_key, "weights")
        
        if(!is.null(cached_weights)) {
          odontocete_sf$proximity_weight <- cached_weights
        } else {
          odontocete_sf$proximity_weight <- calculate_proximity_weights(odontocete_sf)
          save_to_cache(odontocete_sf$proximity_weight, weights_cache_key, "weights")
        }
        
        res_odontocete <- performKDE(
          data_sf = odontocete_sf, 
          species_col = "family_group",
          species_list = "Odontocetes (All Toothed Whales)", 
          weight_col = "proximity_weight",
          window = common_window,
          output_prefix = "sightings_odontocetes_proximity",
          sigma_val = DEFAULT_SIGHTINGS_SIGMA,
          use_bw_diggle = USE_BW_DIGGLE_FAMILY,
          data_source_label = "Sightings (Proximity Weighted)",
          resolution = RESOLUTION,
          buffer_percent = 0.05,
          coastline_for_plot = coastline,
          owa_for_plot = owa,
          study_area_for_plot = study_area
        )
      } else {
        res_odontocete <- performKDE(
          data_sf = odontocete_sf, 
          species_col = "family_group",
          species_list = "Odontocetes (All Toothed Whales)", 
          weight_col = NULL,
          window = common_window,
          output_prefix = "sightings_odontocetes",
          sigma_val = DEFAULT_SIGHTINGS_SIGMA,
          use_bw_diggle = USE_BW_DIGGLE_FAMILY,
          data_source_label = "Sightings (Unweighted)",
          resolution = RESOLUTION,
          buffer_percent = 0.05,
          coastline_for_plot = coastline,
          owa_for_plot = owa,
          study_area_for_plot = study_area
        )
      }
    }
    
    mysticete_sf <- cetacean_sf %>%
      filter(cetacean_family == "mysticete") %>%
      mutate(family_group = "Mysticetes (All Baleen Whales)")
    
    if(nrow(mysticete_sf) > 0) {
      if(USE_PROXIMITY_WEIGHTING_FAMILY) {
        weights_cache_key <- paste0("mysti_weights_", get_cache_hash(list(nrow(mysticete_sf))))
        cached_weights <- load_from_cache(weights_cache_key, "weights")
        
        if(!is.null(cached_weights)) {
          mysticete_sf$proximity_weight <- cached_weights
        } else {
          mysticete_sf$proximity_weight <- calculate_proximity_weights(mysticete_sf)
          save_to_cache(mysticete_sf$proximity_weight, weights_cache_key, "weights")
        }
        
        res_mysticete <- performKDE(
          data_sf = mysticete_sf, 
          species_col = "family_group",
          species_list = "Mysticetes (All Baleen Whales)", 
          weight_col = "proximity_weight",
          window = common_window,
          output_prefix = "sightings_mysticetes_proximity",
          sigma_val = DEFAULT_SIGHTINGS_SIGMA,
          use_bw_diggle = USE_BW_DIGGLE_FAMILY,
          data_source_label = "Sightings (Proximity Weighted)",
          resolution = RESOLUTION,
          buffer_percent = 0.05,
          coastline_for_plot = coastline,
          owa_for_plot = owa,
          study_area_for_plot = study_area
        )
      } else {
        res_mysticete <- performKDE(
          data_sf = mysticete_sf, 
          species_col = "family_group",
          species_list = "Mysticetes (All Baleen Whales)", 
          weight_col = NULL,
          window = common_window,
          output_prefix = "sightings_mysticetes",
          sigma_val = DEFAULT_SIGHTINGS_SIGMA,
          use_bw_diggle = USE_BW_DIGGLE_FAMILY,
          data_source_label = "Sightings (Unweighted)",
          resolution = RESOLUTION,
          buffer_percent = 0.05,
          coastline_for_plot = coastline,
          owa_for_plot = owa,
          study_area_for_plot = study_area
        )
      }
    }
  }
}

#.............................
# PROCESS PAM DATA (WITH CACHING)---------
#.............................

if(RUN_PAM) {
  
  cat("\n=== PROCESSING PAM DATA ===\n\n")
  
  # CRITICAL: Use raw file and process like CE_Collab
  filepath_raw <- 'input/2025/baleen_presence_laura_2025.csv'
  filepath_processed <- 'input/2025/baleen_presence_days_laura_2025.csv'
  
  if(file.exists(filepath_raw)) {
    cat("Using raw PAM file (matching CE_Collab format)\n")
    
    baleen_PA_raw <- read.csv(filepath_raw)
    
    # Process like CE_Collab does (exactly matching KDES_CE_Collab.R lines 335-348)
    baleen_PA <- baleen_PA_raw %>%
      mutate(rec_date = as.Date(rec_date)) %>%
      group_by(site, latitude, longitude, species) %>%
      summarise(
        detection_days = sum(presence, na.rm = TRUE),
        effort_days = n_distinct(rec_date),
        proportion_det = detection_days / effort_days,
        .groups = "drop"
      ) %>%
      filter(effort_days > 0)
    
    cat("Processed raw PAM data like CE_Collab\n")
    cat("  Total station-species combinations:", nrow(baleen_PA), "\n")
    
  } else if(file.exists(filepath_processed)) {
    cat("Warning: Using processed PAM file (may differ from CE_Collab)\n")
    baleen_PA <- read.csv(filepath_processed)
  } else {
    stop("No PAM file found!")
  }
  
  baleen_sf_full <- baleen_PA %>%
    st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
    st_transform(crs = UTM20)
  
  baleen_sf <- st_intersection(baleen_sf_full, study_area)
  
  cat("PAM records:", nrow(baleen_sf), "\n")
  cat("Unique sites:", length(unique(baleen_sf$site)), "\n\n")
  
  species_list_pam <- c("Bb", "Bm", "Bp", "Mn", "Eg", "Ba")
  
  #.............................
  # INDIVIDUAL SPECIES PAM KDEs----
  #.............................
  
  if(RUN_INDIVIDUAL_PAM) {
    cat("=== INDIVIDUAL SPECIES PAM KDEs ===\n\n")
    
    pam_species_names <- c(
      "Bb" = "Sei_Whale",
      "Bm" = "Blue_Whale", 
      "Bp" = "Fin_Whale",
      "Mn" = "Humpback_Whale",
      "Eg" = "North_Atlantic_Right_Whale",
      "Ba" = "Minke_Whale"
    )
    
    for(species_code in species_list_pam) {
      species_data <- baleen_sf %>% filter(species == species_code)
      species_readable <- pam_species_names[species_code]
      
      if(nrow(species_data) > 0) {
        cat("Processing", species_readable, "(n =", nrow(species_data), ")\n")
        
        res_pam_individual <- performKDE(
          data_sf = species_data, 
          species_col = "species",
          species_list = species_code, 
          weight_col = "proportion_det",
          window = common_window,
          output_prefix = paste0("pam_baleen_", species_readable),
          sigma_val = DEFAULT_PAM_SIGMA,
          use_bw_diggle = FALSE,
          data_source_label = "Baleen PAM",
          resolution = RESOLUTION,
          buffer_percent = 0.05,
          coastline_for_plot = coastline,
          owa_for_plot = owa,
          study_area_for_plot = study_area
        )
      }
    }
  }
  
  #.............................
  # GROUPED PAM----
  #.............................
  
  if(RUN_GROUPED_PAM) {
    cat("\n=== GROUPED PAM KDE ===\n\n")
    
    baleen_pam_grouped <- baleen_sf %>%
      group_by(geometry) %>%
      summarise(
        detection_days = sum(detection_days, na.rm = TRUE),
        effort_days = first(effort_days),
        .groups = "drop"
      ) %>%
      mutate(
        proportion_det = detection_days / effort_days,
        baleen_group = "All Baleen Whales (PAM)"
      )
    
    res_baleen_pam <- performKDE(
      data_sf = baleen_pam_grouped, 
      species_col = "baleen_group",
      species_list = "All Baleen Whales (PAM)", 
      weight_col = "proportion_det",
      window = common_window,
      output_prefix = "grouped_baleen_pam",
      sigma_val = DEFAULT_PAM_SIGMA,
      use_bw_diggle = FALSE,
      data_source_label = "PAM",
      resolution = RESOLUTION,
      buffer_percent = 0.05,
      coastline_for_plot = coastline,
      owa_for_plot = owa,
      study_area_for_plot = study_area
    )
  }
}

#.............................
# PROCESS BEAKED WHALE PAM DATA---------
#.............................

if(RUN_BEAKED_PAM) {
  
  cat("\n=== PROCESSING BEAKED WHALE PAM DATA ===\n\n")
  
  beaked_sf <- load_from_cache("beaked_sf_clipped", "pam")
  
  if(is.null(beaked_sf)) {
    beaked_filepath <- 'input/2025/beaked_PAM/beaked_pam_results_2026-01-12.csv'
    beaked_PA <- read_csv(beaked_filepath, show_col_types = FALSE)
    
    beaked_summary <- beaked_PA %>%
      group_by(deployment, station, latitude, longitude, species) %>%
      summarise(
        detection_days = sum(presence, na.rm = TRUE),
        effort_days = n_distinct(rec_date),
        proportion_det = detection_days / effort_days,
        .groups = "drop"
      )
    
    beaked_sf_full <- beaked_summary %>%
      st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
      st_transform(crs = UTM20)
    
    beaked_sf <- st_intersection(beaked_sf_full, study_area)
    
    save_to_cache(beaked_sf, "beaked_sf_clipped", "pam")
    cat("Cached beaked PAM data\n")
  } else {
    cat("Loaded beaked PAM data from cache\n")
  }
  
  cat("Beaked PAM records:", nrow(beaked_sf), "\n\n")
  
  species_list_beaked <- unique(beaked_sf$species)
  
  cat("=== INDIVIDUAL BEAKED WHALE SPECIES PAM KDEs ===\n\n")
  
  beaked_species_names <- c(
    "Ha" = "Northern_Bottlenose_Whale",
    "Mb" = "Sowerbys_Beaked_Whale",
    "Zc" = "Goose_Beaked_Whale",
    "MmMe" = "Trues_Gervais_Beaked_Whale"
  )
  
  for(species_code in species_list_beaked) {
    species_data <- beaked_sf %>% filter(species == species_code)
    species_readable <- beaked_species_names[species_code]
    
    if(nrow(species_data) > 0) {
      cat("Processing", species_readable, "(n =", nrow(species_data), ")\n")
      
      res_beaked_individual <- performKDE(
        data_sf = species_data, 
        species_col = "species",
        species_list = species_code, 
        weight_col = "proportion_det",
        window = common_window,
        output_prefix = paste0("pam_beaked_", species_readable),
        sigma_val = DEFAULT_PAM_SIGMA,
        use_bw_diggle = FALSE,
        data_source_label = "Beaked Whale PAM",
        resolution = RESOLUTION,
        buffer_percent = 0.05,
        coastline_for_plot = coastline,
        owa_for_plot = owa,
        study_area_for_plot = study_area
      )
    }
  }
  
  if(RUN_GROUPED_BEAKED_PAM) {
    cat("\n=== GROUPED BEAKED WHALE PAM KDE ===\n\n")
    
    beaked_pam_grouped <- beaked_sf %>%
      group_by(deployment, station, geometry) %>%
      summarise(
        detection_days = sum(detection_days, na.rm = TRUE),
        effort_days = first(effort_days),
        .groups = "drop"
      ) %>%
      mutate(
        proportion_det = detection_days / effort_days,
        beaked_group = "All Beaked Whales (PAM)"
      )
    
    res_beaked_pam_grouped <- performKDE(
      data_sf = beaked_pam_grouped, 
      species_col = "beaked_group",
      species_list = "All Beaked Whales (PAM)", 
      weight_col = "proportion_det",
      window = common_window,
      output_prefix = "grouped_beaked_pam",
      sigma_val = DEFAULT_PAM_SIGMA,
      use_bw_diggle = FALSE,
      data_source_label = "Beaked Whale PAM",
      resolution = RESOLUTION,
      buffer_percent = 0.05,
      coastline_for_plot = coastline,
      owa_for_plot = owa,
      study_area_for_plot = study_area
    )
  }
}

#.............................
# SUMMARY------
#.............................

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("VERSION 3: Fixed coastline/plot extent\n")
cat("- Coastline cropped to study_area using st_crop()\n")
cat("- Plot extent fixed to study_area bbox\n")
cat("- Matching CE_Collab exactly\n")
cat("Resolution:", RESOLUTION, "m\n")
cat("Window buffer: 5%\n")
cat("Mode:", ifelse(FAST_MODE, "FAST (plots only)", "FULL"), "\n")
cat("Caching:", ifelse(ENABLE_CACHING, "ENABLED", "DISABLED"), "\n")
if(RUN_SIGHTINGS) cat("✓ Sightings processed\n")
if(RUN_PAM) cat("✓ Baleen PAM processed\n")
if(RUN_BEAKED_PAM) cat("✓ Beaked whale PAM processed\n")
cat("\nOutputs:\n")
cat("  Rasters:", raster_dir, "\n")
cat("  Plots:", plot_dir, "\n")
cat("  Shapefiles:", shapefile_dir, "\n")
if(ENABLE_CACHING) cat("  Cache:", CACHE_DIR, "\n")

cat("\n==============================================================================\n")
cat("PROCESSING COMPLETE - MATCHING CE_COLLAB EXACTLY\n")
cat("Coastline cropped | Plot extent fixed | Resolution: 500m\n")
cat("==============================================================================\n")