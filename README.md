# Kachemak Bay Seagrass Monitoring

###### README LAST UPDATED: 2026/11/03

------------------------------------------------------------------------

## Overview

------------------------------------------------------------------------

This is a repository under active development to monitor the seasonal density, cover, and extent of seagrass in Kachemak Bay, Alaska.

### Objectives

- Create a year-round dataset of seagrass density, cover, and extent with a particular focus on the north side of Kachemak Bay.

- Pair presence/absence data with validation methods for further development of a Submerged Aquatic Vegetation Habitat Suitability Model for Kachemak Bay.

### Guiding Questions

- What is the seasonality of seagrass in Kachemak Bay?

- Are seagrass meadows changing in size or density through time?

### Current Status

This project is in active development and will change as additional data and analyses are included. The primary short-term goals for the project include:

- Obtain ongoing seagrass density and cover data from Mud Bay on the north side of Kachemak Bay.

- Survey seagrass extent across the Kachemak Bay intertidal using field sampling and drones.

- Create markdown summary document with major methods, inputs, and outputs.

------------------------------------------------------------------------

## Methods

------------------------------------------------------------------------

### Location

Kachemak Bay is located within lower Cook Inlet in southcentral Alaska, USA. It is a high-latitude estuary comprised of multiple sub-bays and fjords bisected by the Homer Spit that roughly divides it into an 'outer' bay and 'inner' bay. The south side of the bay is largely rocky boulder and bedrock, with some smaller sediment substrate, and abuts the Kenai Mountain Range that provide glacial outflow in the spring and summer. The north side of the bay is largely mud, sand, and cobble, and is fed by several river systems flowing down from higher elevation peatlands. It is a macro-tidal system (\>8 m total tidal range) and tends to be well-mixed with some intermittent stratification dependent upon glacial outflow and season.

### Seagrass Communities

Seagrass (*Zostera marina*) can be found throughout the bay in intertidal and shallow subtidal fine sediments both in the inner bay, and south-shore outer bay, with the largest contiguous meadow located directly east of the Homer Spit in Mud Bay, extending eastward toward Fox River Flats. Seagrass throughout the bay is a mixture of contiguous and patchy distributions, and are known to support many marine and terrestrial species including herring and migratory shorebirds. Seagrass in Kachemak Bay has been survey intermittently since at least 2005, however, these efforts have been spatially and temporally disparate.

### Transect Surveys

To capture the extent and density of seagrass in Kachemak Bay, a combination of qualitative and quantitative surveys were conducted, with a focus on the the Mud Bay seagrass meadow. The intent was to create a method for quickly identifying the extent, relative density, and estimated absolute density of seagrass with a combination of repeated transects, quadrat densities, and stratified or haphazard percent cover estimations.

In brief, transects ranging from approximately 100 m to 1500 m were run in Mud Bay and the adjacent Coal Bay at spring low tide series from May to August 2025, and in Mud Bay and Seldovia Slough from February to October 2026. Quadrats of 0.5 m^2^ for seagrass shoot density were taken at regular intervals along each transect and GPS locations for each quadrat recorded. Five haphazardly selected shoots were selected from in and around each quadrat to measure longest blade length and width. Each quadrat was also assigned a semi-quantitative cover value in the context of the surrounding seagrass cover within a \~3 m radius that included none (0%), trace (1%), sparse (2-10%), thin (11-25%), moderate (26-75%), and thick (76-100%). Regularly spaced cover values were also assigned along transects between quadrats in \~3 m radius coverage areas using the same cover values. Nominal GPS locations for each cover observation were derived as evenly-spaced points between each known quadrat GPS location. Additional quadrat densities and cover values were taken haphazardly off transect lines to expand coverage areas. When observations were taken off a transect, GPS locations were recorded for each observation.

------------------------------------------------------------------------

## Data Dictionary

------------------------------------------------------------------------

Standardized variable specifications, units of measurement, data validation rules, and ecological classification scales for intertidal eelgrass monitoring.

## 1. Study & Dataset Metadata

| Metadata Field | Value |
|:-----------------------------------|:-----------------------------------|
| **Study Context** | Intertidal eelgrass (*Zostera marina*) monitoring & shoot morphometrics |
| **Geographic Scope** | Kachemak Bay, Alaska (Lat: 59.43°N–59.66°N, Lon: -151.71°W to -151.42°W) |
| **Temporal Span** | 2025-05-29 to 2026-10-16 |
| **Total Records** | 1,612 rows |
| **Unique Transects** | 36 transects |
| **Target Species** | *Zostera marina* |
| **Observers / Custodians** | Ross Whippo (NOAA-NCCOS-KBL) |
| **Coordinate System** | WGS 84 (Decimal Degrees) |
| **Primary Sampling Unit** | Transect / Quadrat / Subquadrat |

------------------------------------------------------------------------

## 2. Variable Specifications

| Variable Name | Data Type | Unit / Format | Description | Allowable Values / Scale | Null / Missing Code | Example |
|:----------|:----------|:----------|:----------|:----------|:----------|:----------|
| `date` | Date | `m/d/YYYY` | Survey sampling date | Valid calendar dates (2025–2026) | None (Mandatory) | `5/29/2025` |
| `transect` | Text (Categorical) | Alphanumeric Code | Transect line identifier across sampling sites | `MB01`–`MB32`, `MBX, SLX, SSX, SAX` (off-transect observations end with -X) | None (Mandatory) | `MB01` |
| `quadrat` | Integer / Text | Station Index | Quadrat sampling station number along transect line | `1` – `15`, or `NA` for ad-hoc point observations | `NA` (non-transect ad-hoc points) | `1` |
| `timeBegin` | Time / Integer | `HHMM` (Military) | Field survey start time at sampling station | `0000` – `2400` (Local Alaska Time AKST/AKDT) | `NA` | `0930` |
| `timeEnd` | Time / Integer | `HHMM` (Military) | Field survey end time at sampling station | `0000` – `2400` (Local Alaska Time AKST/AKDT) | `NA` | `1030` |
| `lat` | Decimal Number | Decimal Degrees (°N) | Latitude coordinate of quadrat station (WGS84 datum) | `59.430000` to `59.660000` | `NA` | `59.634715` |
| `lon` | Decimal Number | Decimal Degrees (°W) | Longitude coordinate of quadrat station (WGS84 datum) | `-151.710000` to `-151.420000` | `NA` | `-151.491173` |
| `shoots` | Integer | Count (\# shoots) | Total vegetative eelgrass shoot count in quadrat | Integers $\ge 0$, or `NA` | `NA` (uncounted / trace) | `63` |
| `shootID` | Integer / Text | Sequential Index | Identifier for individual measured shoot within quadrat | `1`, `2`, `3`, `4`, `5`, or `NA` / `nd` | `NA` (no shoot), `nd` (no data) | `1` |
| `shootLength` | Decimal Number | Centimeters (cm) | Maximum sheath-to-tip leaf length of measured shoot | Real numbers $> 0$ | `NA` (no shoot), `nd` (no data) | `36` |
| `shootWidth` | Decimal Number | Centimeters (cm) | Leaf blade width of measured seagrass shoot | Real numbers $> 0$ (typically `0.1` – `0.5` cm) | `NA` (not measured), `nd` (no data) | `0.3` |
| `densityID` | Integer / Text | Partition Index | Subquadrat partition index for density estimation | `1` to `17`, or `NA` | `NA` (not partitioned) | `1` |
| `densityValue` | Text (Categorical) | Qualitative Rank | Qualitative eelgrass density ranking in subquadrat | `none`, `trace`, `sparse`, `thin`, `moderate`, `thick` | `NA` (unrecorded) | `thick` |
| `substrate` | Text (Categorical) | Qualitative Classification | Dominant benthic substrate sediment composition | `mud`, `sand-mud`, `gravel`, `gravel-sand`, `cobble`, `cobble-mud`, `cobble-sand`, `water-sand-mud` | `NA` (unrecorded) | `mud` |
| `collector` | Text | Observer Names | Field scientific personnel conducting observation | Unrestricted descriptive text | `NA` | `Whippo` |
| `notes` | Text | Free Text Remarks | Field observations,reproductive status, environmental notes, tidal limits | Unrestricted descriptive text | Empty cell (No remarks) | `reproductive eelgrass Q9` |

------------------------------------------------------------------------

## 3. Ecological Classification Scales

### Density Scale (`densityValue`)

Estimated over \~3 m radius from position of observer. Extrapolated densities are derived from quadrat shoot counts of 0.5 m^2^ within the area described.

| Value | Ordinal Rank | Meadow Cover Definition | Field Density Extrapolated from Associated Quadrat Shoot Counts |
|:-----------------|:-----------------:|:-----------------|:-----------------|
| **none** | 0 | Absence of eelgrass; bare substrate or macroalgae (*Ulva*) only | 0 shoots recorded, n = 59 |
| **trace** | 1 | Very rare, isolated single shoots or blade fragments (≤ 1% cover) | mean = 249, sd = 545, n = 40 |
| **sparse** | 2 | Scattered, low-density shoot clusters (\~2% – 10% cover) | mean = 758, sd = 1045, n = 24 |
| **thin** | 3 | Light canopy coverage, open meadow patches (\~11% – 25% cover) | mean = 3950, sd = 2114, n = 21 |
| **moderate** | 4 | Well-established eelgrass meadow, medium density (\~26% – 75% cover) | mean = 5756, sd = 3011, n = 32 |
| **thick** | 5 | Dense, continuous eelgrass carpet (\> 75% cover) | mean = 11574, sd = 6898, n = 38 |

------------------------------------------------------------------------

## 4. QA/QC Codes & Survey Conventions

| Code / Syntax | Applicable Variables | Interpretation & Protocol |
|:-----------------------|:-----------------------|:-----------------------|
| `NA` | All variables | **Not Applicable / Not Sampled**: Station or partition where this specific metric was not evaluated (e.g., ad-hoc observation points without quadrat stations). |
| `nd` | `shootID`, `shootLength`, `shootWidth` | **No Data / Not Determined**: Shoots were present within the quadrat, but individual blade dimensions were not measured. |
| `0` | `shoots` | **True Biological Zero**: Quadrat was thoroughly inspected and confirmed to contain 0 eelgrass shoots. |
| `none` | `densityValue` | Confirmed absence of eelgrass in subquadrat partition (substrate bare or *Ulva* macroalgae only). |
| `MBX` , `SLX` , `SBX` , `SAX` | `transect` | Exploratory reconnaissance waypoint sampled outside the established transect grid. |
| `1`–`17` | `densityID` | Partition grid cell index within the sampling frame used for spatial density mapping. |
| `shootLength` | Morphometrics | Sheath meristem to leaf apex blade length measured in centimeters. |
| `shootWidth` | Morphometrics | Mid-blade leaf width measured in centimeters. |
| `lat` / `lon` | Spatial Coordinates | Field handheld GPS coordinates referenced to datum WGS 84. |

\`\`\`
