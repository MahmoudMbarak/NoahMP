# NoahMP Related Code

This repository contains scripts and utilities for working with NoahMP land surface model.

July 2025 Update:
NLDAS recently changed their data format from GRIB to NetCDF. The new extraction scripts  handle this transition while keeping the NoahMP workflow unchanged.

## Scripts

apcp.py - Processes precipitation data with proper accumulation metadata

init.py - Extracts initial condition variables and converts them to GRIB format

extract_nldas.perl - Extracts: Downward longwave radiation flux (W/m²), Downward shortwave radiation flux (W/m²), Surface pressure (Pa) ,Air temperature at 2m (K), Specific humidity at 2m (kg/kg), U-component wind at 10m (m/s), V-component wind at 10m (m/s)

## Usage

### Dependencies
You will need CDO module, wgrib, and python/perl for these scripts.

### Download NLDAS Data from NASA Goddard Earth Sciences Data: https://disc.gsfc.nasa.gov/datasets?keywords=NLDAS

**For Forcing Data** (apcp.py and extract_nldas.perl)
- Dataset: NLDAS_FORA0125_H_2.0
- Navigate to the NLDAS_FORA0125_H_2.0 dataset page
- Download the NetCDF files (format: NLDAS_FORA0125_H.AYYYYMMDD.HHMM.020.nc)

**For Initial Conditions** (init.py)
- Dataset: NLDAS_NOAH0125_H_2.0
- Navigate to the NLDAS_NOAH0125_H_2.0 dataset page
- Download the NetCDF file for your initialization date (format: NLDAS_NOAH0125_H.AYYYYMMDD.HHMM.020.nc)

### 1. Pre-processing NLDAS Forcing Data

#### 1.1 Extract meteorological variables from NLDAS forcing files

**extract_nldas.perl**: Extracts individual meteorological variables from NLDAS NetCDF files

In extract_nldas.perl, users need to modify these key settings:

The processed files will appear in your specified results directories, ready for use with Noah-MP create_forcing scripts.

$data_dir = "/path/to/your/NLDAS_FORA_data/";     # Input directory

$results_dir = "/path/to/output/forcing/";        # Output directory

$wgrib = "/path/to/wgrib";                        # Path to wgrib executable

Set the time period (Julian days):

$day_start = 1;     # January 1st

$day_end = 365;     # December 31st (use 366 for leap years)

Run the script: 
perl extract_nldas.perl

#### 1.2 Extract precipitation data with proper accumulation metadata

**apcp.py**: Processes precipitation data from NLDAS forcing files

In apcp.py, modify these settings:

data_dir = "/path/to/your/NLDAS_FORA_data/"      # Input directory

results_dir = "/path/to/output/forcing/"         # Output directory

wgrib = "/path/to/wgrib"                         # Path to wgrib executable

day_start = 1

day_end = 365

Run the script:
python apcp.py

### 2. Extract Initial Condition Data

**init.py**: Extract initial-only fields (SWE, CanopInt, AvgSurfT, SoilM_0-200cm, SoilT_0-200cm) from NLDAS Noah NetCDF files

Change the initialization date and paths in init.py:

data_dir = "/path/to/your/NLDAS_NOAH_data/"      # Input directory

results_dir = "/path/to/output/init/"            # Output directory

wgrib = "/path/to/wgrib"                         # Path to wgrib executable

date = "20170101"   # Set initialization date (YYYYMMDD)

hh = "01"           # Set initialization hour

Run the script:
python init.py

