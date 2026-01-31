# PAM Overall KDEs (All Months Combined)

This directory contains KDE surfaces for passive acoustic monitoring (PAM) detections 
aggregated across all months.

## Directory Structure

- **rasters/** - GeoTIFF raster files (.tif)
- **shapes/** - 0.90 quantile contour shapefiles (.shp and associated files)

## Processing Details

- **Data**: Autonomous acoustic recorder deployments
- **Weighting**: Detection rate (proportion of deployment days with species detection)
- **Bandwidth**: 10,000 m (fixed)
- **Temporal coverage**: All months combined

## Species Codes

See main README for species code definitions.

## Files

Each species has:
- `rasters/{SPECIES_CODE}.tif` - KDE raster surface
- `shapes/{SPECIES_CODE}_contour90.shp` - 0.90 quantile contour polygon

See ../../../../kde_index.csv for complete metadata.

