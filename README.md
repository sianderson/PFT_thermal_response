# PFT_thermal_response
Marine phytoplankton functional types and their responses to thermal change

These scripts are provided in the interests of open science. If you have questions or find errors, please let us know.

Contact:<br/>
Stephanie I. Anderson<br/>
Graduate School of Oceanography<br/>
University of Rhode Island<br/>
sianderson@uri.edu<br/>


## Directory structure
- [data](data/): has the raw data
- [output](output/): files output by scripts 
- [scripts](scripts/): run using data and output, produce temp, output, figures, and tables
- [figures](figures/): plotted figures
- [SST_data](SST_data/): SST rasters generated for the analyzed time periods


## Analysis  workflow breif overview
- [SST_multimodel_ensemble.R](scripts/SST_multimodel_ensemble.R) to extract SST data for specific time periods and prepare it for use in future analyses
- [Area_on_Earth.R](scripts/Area_on_Earth.R) to generate a raster of grid areas for future calculations
- [Thermal_reaction_norms.R](scripts/Thermal_reaction_norms.R) to fit exponential curves to reaction norms
    - [Curve_comparison.R](scripts/urve_comparison.R) to compare curves
- [Thermal_capacities.R](scripts/Thermal_capacities.R) to compute TSM, WT, and DGE
- [Growth_change.R](scripts/Growth_change.R) to compute proportional growth change for each PFT expected by 2100


## Script Overview

1. [Area_on_Earth.R](scripts/Area_on_Earth.R)
    1. Contains:
        1. functions to calculate lat/lon grid sizes
        2. Generates a raster of grid areas which can be used in later math
    2. Requires
        1. lat/lon inputs

2. [Curve_comparison.R](scripts/urve_comparison.R)
    1. Contains:
        1. Statistical Analyses comparing thermal reaction norms
        2. Extended Figure 2
    2. Requires
        1. Output from [Thermal_reaction_norms.R](scripts/Thermal_reaction_norms.R) ([Isolate_growth_bounds.csv](output/Isolate_growth_bounds.csv))

3. [custom_theme.R](scripts/custom_theme.R)
    1. plot characteristics for use with ggplot

4. [Growth_change.R](scripts/Growth_change.R)
    1. Contains:
        1. Predicted Growth Change by 2100
        2. Figure 5
        3. Calculations for cyanobacteria range expansion
    2. Requires
        1. SST temperature from [SST_multimodel_ensemble.R](scripts/SST_multimodel_ensemble.R)

5. [Isolation_map.R](scripts/Isolation_map.R)
    1. Map of Strain isolation locations (Extended Figure 1)

6. [nbcurve.R](scripts/nbcurve.R)
    1. thermal reaction norm function from Thomas et al 2012

7. [SST_multimodel_ensemble.R](scripts/SST_multimodel_ensemble.R)
    1. Contains:
        1. selects and prepares SST multi-model ensemble data for specific time periods
        2. computes SST change by latitude (change_by_lat.csv)
    2. Requires
        1. SST temperature data downloaded from http://climexp.knmi.nl/ 

8. [Thermal_capacities.R](scripts/Thermal_capacities.R)
    1. Contains:
        1. Modeled temperatures for each strain based on isolation locations
        2. Calculation of thermal capacities
        3. Figure 3
        4. Extended Figure 6
        5. outputs [Thermal_capacity_by_group.csv](output/Thermal_capacity_by_group.csv)

9. [Thermal_reaction_norms.R](scripts/Thermal_reaction_norms.R)
    1. Contains:
        1. Phytoplankton functional type thermal reaction norms (Figure 1)
        2. Exponential curve fits
        3. Q10 calculations (Table 1)
        4. Extended Figures 3 & 5
        5. outputs [Isolate_growth_bounds.csv](output/Isolate_growth_bounds.csv)
