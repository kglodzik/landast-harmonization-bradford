# Landsat Harmonization for Bradford Forest

This repository contains R script and reflectance data used in Deljouei et al. (in revision) to generate correction equations for harmonizing NIR and Red surface reflectance values across Landsat 5, 7, and 8, using Landsat 7 as a bridge. Although Landsat Collection 2 has improved cross-sensor calibration built in, we found implausible increases in Leaf Area Index (LAI, which is derived from NIR and Red) in our study area when transitioning from Landsat 5 to 8. The goal of this script was to enable consistent time series analysis of LAI in Bradford Forest, Florida, across the Landsat 5-8 period of record.

The 'BandValues' datasets in sample_data are NIR and Red band values at randomly generated points, plus tables of date matches between missions; these were the input datasets specific to Bradford Forest for Deljouei et al. (in revision). They are provided here as example inputs to allow users to easily test the script. Users should substitute their own data for other study areas.

## Contents

### /R_Script/
- Landsat_Harmonization_Analysis.R: R script to clean, analyze, and harmonize reflectance values from time-matched Landsat images.
### /sample_data/
- L5toL7MatchesTable_Bradford.csv: List of cloud-free image acquisition date pairs between Landsat 5 and 7 (identified via GEE script) used for Bradford Forest analysis.
- L7toL8MatchesTable_Bradford.csv: List of cloud-free image acquisition date pairs between Landsat 7 and 8 (identified via GEE script) used for Bradford Forest analysis.
- *BandValues_Bradford.csv: Extracted Red and NIR reflectance values at randomly located sample points in Bradford Forest.

## Google Earth Engine (GEE) Script

The imagery used in this analysis was first gathered and filtered using a GEE script. This script:
- Loads cloud-free Landsat 5, 7, and 8 surface reflectance imagery
- Matches acquisition dates within ~8 days
- Exports Red and NIR bands for further analysis

GEE script:  
👉 https://code.earthengine.google.com/ccea0823362ae6d5b4dc6e8b59d60317

## Reproducing the Analysis

1. Run the GEE script linked above to export Red and NIR bands for each Landsat sensor.
2. Generate random sample points within your study area (we used 350-meter spacing to reduce autocorrelation).
3. Extract Red and NIR values for those points using ArcGIS or similar GIS software.
4. Use the R script in this repository to clean, visualize, and harmonize the reflectance values.
