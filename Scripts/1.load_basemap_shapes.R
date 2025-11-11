# load shapefiles for KDE_PAM

pacman::p_load(sp, terra, dplyr, sf, viridis, ggplot2, ggrepel, stringr, here, ggtext, readr, grid,
               pals, tidyr, fuzzyjoin, patchwork,mapsf,classInt,
               ggforce, readr, ggspatial, lubridate, stars, patchwork, scales, RColorBrewer, grafify)

#projections------
UTM20 <- terra::crs("+init=epsg:32620") # CODE FOR UTM Zone 20
# UTM21 <-  terra::crs("+init=epsg:32621") # CODE FOR UTM Zone 21
prj <-  terra::crs("+init=epsg:4326") # EPSG code for WGS84

filepath = 'input/2025/baleen_presence_days_laura_2025.csv'


# read in PAM AMAR STNS for plotting-------

baleen_stns = st_as_sf(read.csv(filepath), coords = c("longitude", "latitude"), crs =4326 )%>%
  st_transform(crs = UTM20)%>%dplyr::group_by(site)%>%dplyr::summarise()

 plot(st_geometry(baleen_stns))

# #beaked
# beaked_stns = read_sf("input/DOY.shp",crs =4326 )%>%st_transform(crs = UTM20)%>%dplyr::group_by(site)%>%dplyr::summarise()
# 
# stations = rbind(baleen_stns, beaked_stns)
# # plot(st_geometry(stations))
# crs(stations)

# bound boxes for plots / study areas -----
GEO_BOUND =  st_bbox( c(xmin = -75,ymin = 38, xmax = -40, ymax =60 ), crs = st_crs(4326))%>%
  st_as_sfc()%>% st_sf()


UTM_BOUND = st_bbox(GEO_BOUND)%>%st_as_sfc()%>% st_sf()%>%st_transform(UTM20)


# bound box SS data -----
Bound_boxB <- st_bbox( c(xmin = -68,ymin = 40, xmax = -50, ymax =48 ), crs = st_crs(4326))
Bound_boxB <- Bound_boxB %>%
  st_as_sfc()%>% #turns the bounding box into a sfc object, that just describes a specific geometry
  st_sf()
Bound_boxBUTM <- Bound_boxB%>%st_transform(UTM20)

ext(Bound_boxBUTM)

# bathy data -------
#for contours use GEBCO
r <- terra::rast("~/CODE/shapefiles/Bathymetry/GEBCO_bathy/gebco_2020.tif")
ext(r)
crs(r)
# plot(r)

#need to downsample bc too big
bathy = terra::aggregate(r, fact = 2)
# bathy_df <- as.data.frame(bathy, xy = T)%>%dplyr::rename(Depth = gebco_2020)

#now crop & project to UTM extent of GEOGRAPHIC area
bathy_UTM = bathy%>%crop(GEO_BOUND)%>%terra::project("EPSG:32620")

bathy<- as.data.frame(bathy_UTM, xy = T)%>%dplyr::rename(Depth = gebco_2020)%>%
  mutate(Depth = ifelse(Depth >=-10, NA, Depth))

#now crop to extent of Gully area
bathy_crop = bathy_UTM%>%crop(Bound_boxBUTM)

crs(bathy_crop)

bathy_cropUTM <- as.data.frame(bathy_crop, xy = T)%>%dplyr::rename(Depth = gebco_2020)%>%
  mutate(Depth = ifelse(Depth >=-10, NA, Depth))

#contours----
cont <- as.contour(r, levels= c( -200, -350, -400, -500,-1000,-2000,-2500,-3000, -3200, -4000, -5000))
cont <- st_as_sf(cont)%>%st_cast("LINESTRING")
cont_UTM = cont%>%st_transform(UTM20)%>%st_intersection(UTM_BOUND)

##MPAs & conservation zones compiled from layers available as shapefiles ----
# OECMS from https://open.canada.ca/data/en/dataset/44769543-7a23-4991-a53f-c2cf7c7a946f
# Ocean's Act MPAs from https://open.canada.ca/data/en/dataset/a1e18963-25dd-4219-a33f-1a38c4971250
# Georges Bank Oil and Gas Exclusion zone from CNSOPB https://callforbids.cnsopb.ns.ca/2016/01/data-environment/gis-information
# Northern bottlenose whale Critical Habitat designations from Species at Risk
# https://open.canada.ca/data/en/dataset/db177a8c-5d7d-49eb-8290-31e6a45d786c


###### Compile all conservation zone information
#OA MPAS
ALL_MPAS = read_sf(here::here("~/CODE/shapefiles/ProtectedAreas/DFO/OA_MPAs/DFO_MPA_MPO_ZPM.shp"))
ALL_MPAS_UTM= ALL_MPAS%>%
  st_transform(UTM20) 
#OECDS
ALL_MPAS2 = read_sf(here::here("~/CODE/shapefiles/ProtectedAreas/DFO/OECMS/DFO_OECM_MPO_AMCEZ.shp"))
ALL_MPAS2_UTM= ALL_MPAS2%>%
  st_transform(UTM20) 
#O&G exclusions
gbPA = read_sf(here::here("~/CODE/shapefiles/ProtectedAreas/OilnGas/Georges_Bank_exclusion_zone.shp"))
gbPA_UTM= gbPA%>%
  st_transform(UTM20) 
#NAFO Closures
#from http://www.atlas-horizon2020.eu/layers/geonode:a__2019_Closures_sponge_coral

Coral = read_sf(here::here("~/CODE/shapefiles/ProtectedAreas/NAFO_closures/a__2019_Closures_sponge_coral.shp"))
Coral_UTM= Coral%>%
  st_transform(UTM20) 

#FUNDIAN AOI
Fundian <- read_sf(here::here("~/CODE/shapefiles/ProtectedAreas/DFO/FundianAOI/FundianChannel_BrownsBank_AOI_poly.shp"))
Fundian_UTM=Fundian%>%
  st_transform(UTM20)

##CH - No Zone 1-------
NBW_CH <- read_sf(here::here("~/CODE/shapefiles/SAR_CH/NBW_CH/NorthernBottlenoseWhale_CH.shp"))
NBW_CH_UTM=NBW_CH%>%
  st_transform(UTM20)

#read in Gully Zones, transform to UTM Zone 20 
Gully_UTM <- read_sf(here::here("~/CODE/shapefiles/ProtectedAreas/DFO/Gully/Gully_MPA.shp"))%>%st_transform(UTM20)
GullyZ1= Gully_UTM%>%filter(NAME == "Zone 1")

# Land ------
#all countries

land <- read_sf(here::here("~/CODE/shapefiles/coastline/worldcountries/ne_50m_admin_0_countries.shp"))%>%dplyr::filter(CONTINENT == "North America")
landUTM = land%>%st_transform(UTM20) %>%st_intersection(UTM_BOUND)

#NAFO------------

NAFO <- read_sf(here::here("~/CODE/shapefiles/nafo_divisions/divisions/NAFO_Divisions_2021_poly_clipped.shp"))
NAFOUTM = NAFO%>%st_transform(UTM20) %>%st_intersection(UTM_BOUND)
# NAFO_label = sf::st_centroid(NAFO)
# # get coordinates and add back to dataframe
# NAFO_label$Long<-st_coordinates(NAFO_label)[,1] # get coordinates
# NAFO_label$Lat<-st_coordinates(NAFO_label)[,2] # get coordinates

# NAFO_label = NAFO_label%>%filter(!is.na(ZONE), ZONE !="3Pn", ZONE != "4Vn", ZONE != "4S",  ZONE != "4R")

#EEZ------------
EEZ <- read_sf("~/CODE/shapefiles/EEZ/EEZ_can.shp")
EEZ_label = st_point_on_surface(EEZ)

EEZ_UTM = EEZ%>%st_transform(UTM20)%>%st_intersection(UTM_BOUND)

#sable Island
sable <- read_sf("~/CODE/shapefiles/coastline/Sable/Sable Island Low Water Mark 1998.shp")%>% 
  st_transform(4326)
sable_label = st_point_on_surface(sable)
# get coordinates and add back to dataframe
sable_label$Long<-st_coordinates(sable_label)[,1] # get coordinates
sable_label$Lat<-st_coordinates(sable_label)[,2] # get coordinates
