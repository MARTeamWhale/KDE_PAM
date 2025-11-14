#plot WSDB species sightings

source("Scripts/1.load_basemap_shapes.R")

pacman::p_load(sf, tidyverse, readxl, here, leaflet, scales, terra, ggrepel, viridis, ggspatial)

# INPUTS NEED TO BE CSV files, file names are specified in quotations here:-----
input_file <-"WSDB_DFO_Sep2025.csv"
species = "speciesCodes.csv" 

whale_data =  read_csv(here("input/2025/sights/", input_file), col_types = cols(.default ="c"))

# check species, add common names and scientific names based on codes-----
SP_data <- read_csv(here("input", species), col_types = cols(.default ="c"))
WS_data = left_join(whale_data, SP_data)

WS_data%>%group_by(COMMONNAME)%>%summarise(sp_n = n())

# Map  sightings records---------

# bound boxes for study area -----
bbox = st_bbox( c(xmin = -67.5,ymin = 40, xmax = -55, ymax =47.5 ), crs = st_crs(4326))%>% st_as_sfc()%>% #turns the bounding box into a sfc object, that just describes a specific geometry
  st_sf()

# Create point shapefile from records
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

# Extract unique species names
unique_species <- unique(WS_data_sf$COMMONNAME)

# Determine the number of unique species
num_species <- length(unique_species)

# Generate color and shape palettes
pal <- scales::hue_pal()(num_species)
sp_pts <- seq(1, num_species + 15, length.out = num_species) # Adjust shape values

# Loop over each species and create separate plots
for (i in 1:num_species) {
  species_name <- unique_species[i]
  species_data <- WS_data_sf[WS_data_sf$SCIENTIF == species_name, ]
  
  # Create plot for the current species
  Map_Species <- ggplot() +
    geom_sf(data = cont %>% dplyr::filter(level %in% c(-200, -1000, -2500)), col = "light gray", linewidth = 0.2) +
    geom_sf(data = land, color = "grey50", fill = "grey50") +
    geom_sf(data = species_data, aes(color = SCIENTIF, fill = SCIENTIF, shape = SCIENTIF), size = 2, stroke = .5, alpha = 0.65) +
    ggtitle(bquote("Sightings of" ~ italic(.(species_name)))) +  # Add title with italicized species name
    
     theme_bw(base_size = 12) +
    xlab("") + ylab("") +
    scale_color_manual(values = pal[i], name = "") +
    scale_fill_manual(values = pal[i], name = "") +
    scale_shape_manual(values = sp_pts[i], name = "") +
    theme(axis.text.x = element_text(hjust = .75), legend.position = "none") +
    scale_y_continuous(breaks = seq(40, 50, by = 2.5)) + 
    scale_x_continuous(breaks = seq(-65, -55, by = 5)) +
    coord_sf(xlim = c(-70, -55), ylim = c(48, 39.5), expand = F)+
    annotation_scale(location = "br", width_hint = 0.25)  # Add scale bar
  
  
  # Save the plot as a PNG
  file_name <- paste0("output/Figs/sights/Sightings_", species_name, ".png")
  ggsave(file_name, Map_Species, height = 7.5, width = 10, units = "in", dpi = 300)
}

#plot facet of all species ---------

#pallete
predpal="viridis"
pal = get(predpal)(4)
show_col(pal)

sp_pts = c(22,17,19, 8)

labs = factor(WS_data_sf$SCIENTIF)
Map_Species = ggplot() +
  
  # add contours--
  geom_sf(data = cont%>%dplyr::filter(level %in% c(-200, -1000, -2500)), col = "light gray", linewidth  = 0.2) +
  
  #land
  geom_sf(data = land, color= "grey50", fill="grey50") +
  
  
  #add ASO
  geom_sf(data = WS_data_sf, aes(col = SCIENTIF, fill = SCIENTIF, shape = SCIENTIF), size = 2,stroke = .5, alpha = 0.65)+
  facet_wrap(~SCIENTIF)+
  theme_bw(base_size = 12)+
  
  xlab("") + ylab("") +
  # ggtitle(eq2)+
  scale_color_manual(values = pal, name = "")+
  scale_fill_manual(values = pal, name = "")+
  scale_shape_manual(values = sp_pts, name = "")+
  #add themes to specify legend location and other attributes of display here
  theme(axis.text.x=element_text(hjust=.75), 
        legend.position= "none",
  )+
  # # set map limits
  scale_y_continuous(breaks = seq( 40,50, by = 2.5)) + 
  scale_x_continuous(breaks = seq(-65, -55, by = 5)) +
  coord_sf(xlim = c(-70, -55), ylim = c(48, 39.5), expand = F)


Map_Species

ggsave(here::here("output/FIGS/sights/Whale_sightings_Map.png"), Map_Species, height = 7.5, width = 10, units = "in", dpi = 300 )

