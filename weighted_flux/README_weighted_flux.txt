# README_weighted_flux
This directory contains functions for downloading data and calculating time-integrated sea-air gas fluxes (positive flux is from the sea to the air), including weighting based on historical sea ice and wind speed data, and accounting for the mixed layer depth.

The demonstration functions are listed below:

DEMO_DOWNLOAD.M demonstrates how to download the historical data.

DEMO_WEIGHTED_GAS_FLUXES.M demonstrates how to calculate the sea-air fluxes of N2O and CH4 and compare the instantaneous flux with a 30-day weighting and a 60-day weighting using a small example dataset (available at DATA/DEMO_GAS_DATA.CSV).

The download functions require installation of GNU Wget (available at https://www.gnu.org/software/wget/ and https://eternallybored.org/misc/wget/). The executable file wget.exe should be moved to the C:/Windows/System32 directory and Matlab should be restarted before running DEMO_DOWNLOAD.M.