# KDE Analysis Outputs for Collaborative Project

This directory contains Kernel Density Estimation (KDE) surfaces and associated shapefiles 
for cetacean distribution analysis in the Eastern Scotian Shelf study area.

## Directory Structure

```
collab_CE/
├── shapefiles/
│   └── study_area/           # Study area boundary shapefile
├── kde_surfaces/
│   ├── species/
│   │   ├── rasters/          # Individual species KDE rasters (GeoTIFF)
│   │   └── shapes/           # Individual species 0.90 contour shapefiles
│   ├── guilds/
│   │   ├── rasters/          # Guild-level KDE rasters (GeoTIFF)
│   │   └── shapes/           # Guild-level 0.90 contour shapefiles
│   └── PAM/
│       ├── overall/
│       │   ├── rasters/      # Overall PAM KDE rasters (all months)
│       │   └── shapes/       # Overall PAM 0.90 contour shapefiles
│       └── monthly/
│           ├── rasters/      # Monthly PAM KDE rasters
│           └── shapes/       # Monthly PAM 0.90 contour shapefiles
├── KDE_maps/                 # Quality check visualizations
└── kde_index.csv             # Master metadata file
```

## File Naming Convention

### Rasters (GeoTIFF)
- Species sightings: `{SPECIES_NAME}.tif`
- Guilds: `{GUILD_NAME}.tif`
- PAM overall: `{SPECIES_CODE}.tif`
- PAM monthly: `{SPECIES_CODE}_month{MM}.tif`

### Shapefiles (0.90 quantile contours)
- Corresponding shapefile: `{RASTER_NAME}_contour90.shp`
- Located in parallel `/shapes/` directories

## KDE Parameters

### Sightings Data
- **Bandwidth**: 10,000 m (10 km)
- **Kernel**: Gaussian
- **Resolution**: ~500 m
- **Weighting**: Proximity-weighted for individual species (isolated observations downweighted)
- **Normalization**: Normalized to [0, 1] by species maximum
- **CRS**: EPSG:32620 (WGS84 UTM Zone 20N)

### PAM Data
- **Bandwidth**: 10,000 m (10 km)
- **Kernel**: Gaussian
- **Resolution**: ~500 m
- **Weighting**: Detection rate (proportion of days with detections)
- **Normalization**: Normalized to [0, 1] by species maximum
- **CRS**: EPSG:32620 (WGS84 UTM Zone 20N)

## Species Codes (PAM)

- **Bb** - Sei Whale (*Balaenoptera borealis*)
- **Bm** - Blue Whale (*Balaenoptera musculus*)
- **Bp** - Fin Whale (*Balaenoptera physalus*)
- **Mn** - Humpback Whale (*Megaptera novaeangliae*)
- **Ba** - Minke Whale (*Balaenoptera acutorostrata*)
- **Eg** - North Atlantic Right Whale (*Eubalaena glacialis*)

## Metadata Index

The file `kde_index.csv` contains complete metadata for all KDE outputs including:
- Dataset type (sightings_individual, sightings_guild, PAM_monthly, PAM_overall)
- Species/group name
- Month (for monthly analyses)
- Number of input points/stations
- Raster filename
- Contour shapefile filename
- CRS
- Spatial resolution
- Bandwidth
- Weighting method

## Quality Check Maps

The `KDE_maps/` directory contains PNG visualizations of KDE surfaces with:
- Sqrt-transformed density values for improved visibility
- Coastline overlay
- Offshore Wind Area (OWA) boundaries
- PAM station locations (for PAM analyses)
- Study area boundary

## Data Sources

- **Sightings**: Combined deduplicated dataset from WSDB and aerial surveys (2015-present)
- **PAM**: Autonomous acoustic recorder deployments (baleen/ beaked whale presence)

## Contact

For questions about these outputs, contact:
Laura Joan Feyrer (ljfeyrer@dal.ca)

## Date Generated

2026-01-15

