cat("\n=== CREATING SEASONAL FACET PLOTS ===\n\n")

plotSeasonalQuantileFacets <- function(shapefile_dir,
                                       land,
                                       contours = NULL,
                                       output_dir = plot_output,
                                       predpal = "mako",
                                       exclude_quantiles = NULL) {
  
  shapefiles <- list.files(
    shapefile_dir,
    pattern = "KDE_pam_.*_season_.*\\.shp$",
    full.names = TRUE
  )
  
  if (!length(shapefiles)) {
    cat("No seasonal shapefiles found in", shapefile_dir, "\n")
    return(NULL)
  }
  
  cat("Reading", length(shapefiles), "seasonal shapefiles...\n")
  
  all_data <- lapply(shapefiles, function(f) {
    tryCatch({
      shp <- st_read(f, quiet = TRUE)
      
      filename     <- basename(f)
      species_code <- str_extract(filename, "(?<=pam_)[A-Za-z]+")
      season_name  <- str_extract(filename, "(?<=season_)[A-Za-z]+")
      
      if (!"species" %in% names(shp)) shp$species <- species_code
      if (!"Season"  %in% names(shp)) shp$Season  <- factor(season_name,
                                                            levels = season_levels)
      shp
    }, error = function(e) {
      NULL
    })
  })
  
  all_data <- all_data[!sapply(all_data, is.null)]
  if (!length(all_data)) {
    cat("No valid seasonal shapefiles could be read\n")
    return(NULL)
  }
  
  combined_data <- do.call(rbind, all_data)
  
  # optional drop of certain quantiles
  if (!is.null(exclude_quantiles) && length(exclude_quantiles)) {
    combined_data <- combined_data %>%
      filter(!Quantil %in% exclude_quantiles)
  }
  
  if (!nrow(combined_data)) {
    cat("No polygons left after filtering quantiles\n")
    return(NULL)
  }
  
  # fixed bbox like monthly script (Scotian Shelf region)
  bbox <- st_bbox(c(xmin = -68, ymin = 41, xmax = -55, ymax = 47.8),
                  crs = st_crs(4326)) %>%
    st_as_sfc() %>%
    st_sf() %>%
    st_transform(st_crs(land))
  
  bbox <- st_bbox(bbox)
  xlims <- c(bbox["xmin"], bbox["xmax"])
  ylims <- c(bbox["ymin"], bbox["ymax"])
  
  combined_data$Season <- factor(combined_data$Season,
                                 levels = season_levels,
                                 ordered = TRUE)
  
  unique_species <- sort(unique(combined_data$species))
  
  for (sp in unique_species) {
    species_data <- combined_data %>% filter(species == sp)
    if (!nrow(species_data)) next
    
    species_display <- if (!is.null(species_names[[sp]])) {
      species_names[[sp]]
    } else {
      sp
    }
    
    cat("Creating seasonal facet plot for", species_display, "\n")
    
    # quantiles present for this species
    uq <- sort(unique(species_data$Quantil))
    uq <- uq[!is.na(uq)]
    
    # ensure ordered factor
    qnums <- as.numeric(gsub("%", "", uq))
    order_idx <- order(qnums)
    uq <- uq[order_idx]
    
    species_data$Quantil <- factor(species_data$Quantil,
                                   levels = uq,
                                   ordered = TRUE)
    
    n_q <- length(uq)
    pal <- rev(get(predpal)(n_q + 1))
    
    p <- ggplot() +
      theme_bw()
    
    if (!is.null(contours)) {
      p <- p +
        geom_sf(
          data = contours %>%
            dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)),
          col  = "grey70",
          linewidth = 0.2
        )+
        geom_sf(data = baleen_stns, color = "grey20", fill = "grey20", shape = 19, alpha = .5) 
    }
    
    p <- p +
      geom_sf(data = species_data, aes(fill = Quantil),
              col = NA, alpha = 0.8, na.rm = TRUE) +
      geom_sf(data = land, color = NA, fill = "grey20") +
      facet_wrap(~ Season, ncol = 2) +
      scale_fill_manual(
        values       = pal,
        na.value     = NA,
        na.translate = FALSE,
        name         = "Quantile",
        labels       = uq,
        drop         = FALSE
      ) +
      coord_sf(
        lims_method = "orthogonal",
        xlim = xlims,
        ylim = ylims,
        crs  = st_crs(land),
        expand = TRUE
      ) +
      labs(
        title = paste("Seasonal Acoustic Presence -", species_display),
        x = "",
        y = ""
      ) +
      theme(
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        legend.position   = "right",
        legend.background = element_rect(fill = NA),
        legend.key        = element_rect(fill = NA),
        plot.margin       = margin(5, 5, 5, 5),
        plot.title        = element_text(hjust = 0.5, size = 16, face = "bold"),
        strip.text        = element_text(size = 11, face = "bold"),
        strip.background  = element_rect(fill = "grey90"),
        axis.title        = element_blank(),
        axis.text         = element_blank(),
        axis.ticks        = element_blank()
      ) +
      annotation_scale(location = "br", width_hint = 0.25,
                       text_cex = 0.6, bar_cols = c("grey40", "white"))
    
    out_png <- file.path(output_dir,
                         paste0("seasonal_quantile_facets_", sp, ".png"))
    
    ggsave(out_png, p,
           width = 12, height = 10, dpi = 300, bg = "white", units = "in")
    
    cat("Saved:", out_png, "\n\n")
  }
}

#-------------------------------
# STEP 3: run plotting
#-------------------------------

# projection same as your monthly script
UTM20 <- terra::crs("+init=epsg:32620")

land <- read_sf(here::here("~/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp")) %>%
  dplyr::filter(CONTINENT == "North America") %>%
  st_transform(UTM20)

plotSeasonalQuantileFacets(
  shapefile_dir    = shapefile_output,
  land             = land,
  contours         = if (exists("cont_UTM")) cont_UTM else NULL,
  output_dir       = plot_output,
  predpal          = "mako",
  exclude_quantiles = NULL  # or c("80%") if you want only 90/95
)

cat("\n=== SEASONAL QUANTILE ANALYSIS COMPLETE ===\n")
cat("Shapefiles saved to:", shapefile_output, "\n")
cat("Plots saved to:",      plot_output, "\n")
