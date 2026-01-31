# Individual Species KDEs (Sightings Data)

This directory contains KDE surfaces for individual cetacean species with ≥25 sighting records.

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Processing Details

- **Data**: Combined sightings from WSDB and aerial surveys (2015-present)
- **Weighting**: Proximity-weighted (isolated observations receive lower weight)
- **Proximity parameters**:
  - K nearest neighbors: 8
  - Reference distance: 5000 m
  - Minimum weight: 0.05
- **Bandwidth**: 10,000 m (fixed)

## Files

Each species has:
- `rasters/{SPECIES}.tif` - KDE raster surface
- `shapes/{SPECIES}_contour90.shp` - 0.90 quantile contour polygon

See ../../../kde_index.csv for complete metadata.

