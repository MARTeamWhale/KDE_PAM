# Load packages-----
pacman::p_load(sf, dplyr, here, terra, smoothr)

# ---- 1.Read in existing shapefiles--------------------------

# paths - change these to your actual files
eez_path  <- here::here("shapefiles/EEZ/EEZ_can.shp")
  nafo_path  <- "shapefiles/NAFO/NAFO_Divisions_2021_poly_clipped.shp"
land_path   <- here::here("/Users/chirp/CODE/KDE_PAM/shapefiles/Canada/NE_10mLand.shp")
# bathy_path <- "/Users/chirp/CODE/shapefiles/Bathymetry/GEBCO_bathy/gebco_2020.tif"
WEA_path = "shapefiles/WEA/Designated_WEAs_25_07_29.shp"
nafo <- st_read(nafo_path)
eez  <- st_read(eez_path)
land <- st_read(land_path)
WEA <- st_read(WEA_path)
bathy <- bathy_cropUTM #from load basemap
OSW_studyarea = st_read("/Users/chirp/CODE/shapefiles/wind/RA_OSW_Study_Area/Study_Area_NS.shp")

st_crs(OSW_studyarea)


# Pick a common CRS (WGS84 here, use UTM20 or whatever you are mapping in)
target_crs <- 32620

nafo <- st_transform(nafo, target_crs)
eez  <- st_transform(eez, target_crs)
OSW_studyarea  <- st_transform(OSW_studyarea, target_crs)
WEA  <- st_transform(WEA, target_crs)



# --2. Filter NAFO to 3V and 3W, union to one----------------------

# Adjust DIV / DIVISION to match your attribute name
nafo_ess <- nafo |>
  filter(Division %in% c("4V", "4W"))

nafo_ess_union <- nafo_ess |>
  st_make_valid() |>
  st_union() |>
  st_as_sf() |>
  mutate(id = 1)  # simple id so it has an attribute

plot(st_geometry(nafo_ess_union))

# --- 3. Clip to Regional assessment area to stay inside Canada------------------------------------------


ess_nafo_RA <- st_intersection(nafo_ess_union, OSW_studyarea)

 plot(st_geometry(ess_nafo_RA))


# # ---- 4. Remove deep water as not economically feasible for development---------------------------------------------
# 
#  # Create a polygon mask for depth ≤ 400 m
# 
#  # Note: depth values are negative, so threshold is -400
#  shallow_mask <- bathy <= -400   # returns TRUE for deep water; we invert
#  
#  # Create shallow area raster (depth < 400 m)
#  shallow_area <- bathy >= -400   # TRUE where water is shallower than 400m
#  
#  # Convert raster mask to polygons
#  shallow_poly <- as.polygons(shallow_area, dissolve = TRUE) |> 
#    st_as_sf() |> 
#    st_make_valid()
#  
#  # Keep only TRUE pixels
#  shallow_poly <- shallow_poly |> filter(gebco_2020 == 1)
#  
#  # Clip ESS study area to <400 m --------------------------------
#  
#  ess_shallow <- st_intersection(ess_nafo_RA, shallow_poly) |> 
#    st_make_valid() |> 
#    mutate(depth_filter = "<400m")
#  
#  plot(st_geometry(ess_shallow))
#  
 
 # Clip to main polygon (largest area)-----
 
 # 2. Explode multipolygons to single polygons
 ess_parts_utm <- st_cast(ess_nafo_RA, "MULTIPOLYGON")
 ess_parts_utm <- st_cast(ess_parts_utm, "POLYGON")
 
 # 3. Compute area in km2 OUTSIDE any pipe
 #    Use sf::st_area explicitly so we don't hit some other st_area method
 areas <- sf::st_area(ess_parts_utm)
 ess_parts_utm$area_km2 <- as.numeric(areas) / 1e6
 
 # 4. Keep only the largest polygon
 ess_parts_utm <- ess_parts_utm[order(ess_parts_utm$area_km2, decreasing = TRUE), ]
 
 # First row is largest
 ess_main_utm <- ess_parts_utm[1, ]
 
 # Union just in case
 ess_main_utm <- st_union(ess_main_utm)
 ess_main_utm <- st_as_sf(ess_main_utm)
 
 plot(st_geometry(ess_main_utm))
 
 
 # 4. Remove interior holes
 #    threshold is in square metres (CRS is UTM). Set it huge so all holes are filled.
 ess_main_utm <- smoothr::fill_holes(
   ess_main_utm,
   threshold = 1e12   # 1,000,000 km^2, safely bigger than your shelf polygon
 )
 plot(st_geometry(ess_main_utm))
 
 # 5. Simplify
 ess_simple_utm <- st_simplify(
   ess_main_utm,
   dTolerance = 2000,
   preserveTopology = TRUE
 )
 plot(st_geometry(ess_simple_utm))
 
 # # 6. Optional smoothing
 # ess_smooth_utm <- smoothr::smooth(
 #   ess_simple_utm,
 #   method = "ksmooth",
 #   smoothness = 2
 # )
 # plot(st_geometry(ess_smooth_utm))
 
 # 7. Save
 ess_study_final <- ess_simple_utm
 ess_study_final$study_area <- "ess_simple_utm"
 
 st_write(
   ess_study_final,
   "shapefiles/studyArea/ESS_study_area_simple.shp",
   delete_dsn = TRUE
 )
 
 
 #make a map to visualize extent-----
 bbox <- st_bbox(ess_simple_utm)
 
 xmin <- bbox["xmin"]
 xmax <- bbox["xmax"]
 ymin <- bbox["ymin"]
 ymax <- bbox["ymax"]
 
 xlims <- c(xmin-100000, xmax)
 ylims <- c(ymin, ymax)
 
 #for the legend
 OSW_studyarea$type        <- "Offshore wind regional assessment area"
 ess_study_final$type      <- "Eastern Scotian Shelf study area"
 WEA$type                   <- "Wind Energy AOIs"
 
 gg_map <- ggplot() +
    theme_bw() +
    
    # bathy
    geom_sf(
       data = cont_UTM %>% 
          dplyr::filter(level %in% c(-200, -400, -1000, -2500, -3200)),
       col = "grey50",
       linewidth = 0.2
    ) +
    
    # land
    geom_sf(data = landUTM, color = NA, fill = "grey50") +
    
    # OSW study area
    geom_sf(
       data = OSW_studyarea,
       aes(fill = type, color = type),
       alpha = 0.5,
       linewidth = 0.6
    ) +
    
    # ESS < 400 m
    geom_sf(
       data = ess_study_final,
       aes(fill = type, color = type),
       alpha = 0.5,
       linewidth = 0.6
    ) +
    
    # Designated wind energy areas
    geom_sf(
       data = WEA,
       aes(fill = type, color = type),
       alpha = 0.5,
       linewidth = 0.6
    ) +
    
    scale_fill_manual(
       name   = "",
       values = c(
          "Offshore wind regional assessment area" = "grey90",
          "Eastern Scotian Shelf study area"              = "lightblue",
          "Wind Energy AOIs"        = "pink"
       )
    ) +
    scale_color_manual(
       name   = "",
       values = c(
          "Offshore wind regional assessment area" = "black",
          "Eastern Scotian Shelf study area"              = "blue",
          "Wind Energy AOIs"        = "red"
       )
    ) +
    
    coord_sf(
       lims_method = "orthogonal",
       xlim = xlims,
       ylim = ylims,
       crs  = UTM20,
       expand = TRUE
    ) +
    
    labs(x = "", y = "") +
    
    theme(
       legend.position      = c(0.98, 0.02),  # bottom right inside plot
       legend.justification = c("right", "bottom"),
       legend.background    = element_rect(fill = "white", color = "grey80"),
       legend.key.size      = unit(0.6, "cm"),
       legend.title = element_blank()
    )
 
 
 #visualize & save plot
 
 gg_map
 ggsave("output/FIGS/proposed_StudyAreamap.png", gg_map, dpi = 300)
 
 