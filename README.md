# ET Synchrony
This repository contains the full workflow for the analysis of ET synchrony project.
A test dataset of daily scale AmeriFlux data is provided in 00_Data

## 01_Data_processing
### 01.1 QC and pre-processing
`01.1_AMF_pre-processing.R` conducts QC for raw AmeriFlux dataset, aggregate half-hourly data to hourly, output variables include:
- **SWC**: Soil Water Content (Unit: %)
- **TA**: Air temperature (Unit: degree C)
- **VPD**: Vapor pressure deficit (Unit: hPa)
- **P_F**: Precipitation (Unit: mm)
- **WS**: Wind speed (Unit: ms-1)
- **USTAR**: Friction velocity (Unit: ms-1)
- **NETRAD**: Net radiation (Unit: Wm-2)
- **TS**: Soil temperature (Unit: degree C)
- **LE_F**: Latent heat flux (Unit: Wm-2)
### 01.2 Soil variables
`01.2_AMF_secondary_variables.R` matches soil properties with AmeriFlux sites
### 01.3 Phenology calculation
`01.4_LAI_processing.R` calculates phenology based on LAI
### 01.4 Match explanatory data
`01.7_Explanatory_variables_extraction.R` extracts explanatory variables from the remote sensing datasets
### 01.5 Summary of variables
`01.8_Variable_summary.R` calculates statistics for Shannon entropy across environmental groups

## 02_TE_implementation
### TE core codes
`TE_core_codes.R` contains the core functions for TE implementation, the functions are adapted from Edom Moges
### TE calculation across sites
`02_TE_implementation_all_sites.R` calculates TE across AmeriFlux sites

## 03_PN_construction
### Synchrony metrics extraction
`03.1_Synchrony_metrics_extraction.R` extracts synchrony metrics
### Synchrony metrics analysis
`03.2_Synchrony_metrics_exploration_V4.R` analyzes synchrony metrics
### Process network construction
`03.6_PN_plotting.R` constructs the process network
