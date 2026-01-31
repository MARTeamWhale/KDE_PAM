# Guild-Level KDEs (Sightings Data)

This directory contains KDE surfaces for cetacean guilds (functional groups).

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Guilds

- **Odontocetes**: All toothed whales
- **Mysticetes**: All baleen whales

## Processing Details

- **Data**: Combined sightings from WSDB and aerial surveys (2015-present)
- **Weighting**: None (guilds have naturally broad distributions)
- **Bandwidth**: 10,000 m (fixed)

## Files

Each guild has:
- `rasters/{GUILD}.tif` - KDE raster surface
- `shapes/{GUILD}_contour90.shp` - 0.90 quantile contour polygon

See ../../../kde_index.csv for complete metadata.

