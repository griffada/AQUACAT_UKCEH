#README for deconstructed texmex for large data.
## Written by Adam Griffin 2021-01-04.
## Developed for AquaCAT project in 2020.
# This code is not for widespread distribution.

This readme outlines the process for using the script 109_HeffTawn_Modelling.R, including the required data and the expected outputs.
The script uses the basis of the texmex package as on github to simulate new extreme spatially consistent events across a region given a set of observed events with magnitude and rarity.

Packages required:
texmex
extRemes
reshape2
tidyverse
# ncdf4      called by setup_script_00.R
# raster     called by setup_script_00.R
# fields     called by setup_script_00.R
# rgdal      called by setup_script_00.R

INPUTS:
thresh0 - cellwise thresholds of flow for a given region.
eventFlow - a measurement of flow for each cell for each event (one event per row, one location per column.
eventDpe - a measurement of daily exceedence probability for each cell for each event. (This is because the present work is based on mean daily flow.)
eventApe - a measurement of annual exceedence probability for each cell for each event. 
mqu, dqu: marginal and dependency quantiles: fit GPD and conditional exceedence above these quantiles. See mex documentation for more.
nSample - number of desired events.
mult - tuning parameter. If it falls over due to not enough memory, make it smaller. If it can't generate enough events, make it bigger.


Method:

setup_script_00.R looks for specific folders, and generates strings based on command line arguments for RCM, time period and region. This can be skipped, but filenames throughout may need to be modified, both input and output.

The script should run as it is. However, the functions in 07c_texmex_slimline.R rely on certain objects being in the global environment:
    mexTransform_slim needs DATA and MODELS,
    mexMonteCarlo_slim needs MODELS, COEFFS, Z, DATA, TRANSFORMED
There are also three more crucial functions:
    migpd_slim
    predict.mex_slim
    mexDependence_slim
which are defined when 07c_texmex_slimline.R is called. They all correspond to an equivalent mex function (minus the _slim suffix).

OUTPUTS:
The script saves MODELS, TRANSFORMED, DEPENDENCE, COEFFS, Z as .RDS files if the process fully completes.
It also saves MCS, a table of new events as for eventFlow, as a .csv file.

