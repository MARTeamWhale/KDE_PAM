# PAM Monthly KDEs

This directory contains monthly KDE surfaces for passive acoustic monitoring (PAM) detections.

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Processing Details

- **Data**: Autonomous acoustic recorder deployments
- **Weighting**: Detection rate (proportion of deployment days with species detection)
- **Bandwidth**: 10,000 m (fixed)
- **Temporal resolution**: Monthly

## File Naming

Format: `{SPECIES_CODE}_month{MM}.tif`
- MM = 01 (Jan) through 12 (Dec)

## Species Codes

See main README for species code definitions.

## Files

Each species-month combination has:
- `rasters/{SPECIES_CODE}_month{MM}.tif` - KDE raster surface
- `shapes/{SPECIES_CODE}_month{MM}_contour90.shp` - 0.90 quantile contour polygon

See ../../../../kde_index.csv for complete metadata.

