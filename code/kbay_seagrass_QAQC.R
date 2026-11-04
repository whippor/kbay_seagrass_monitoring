#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             
# KBay Seagrass QAQC                                                   
# Script created 2026-10-22                                                  
# Data source: NOAA-NCCOS-KBL                                                 
# R code prepared by Ross Whippo                                              
# Last updated 2026-11-03                                                   
#                                                                             
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:
# Script for the QAQC of seagrass density, cover, and extent data
# collected from Kachemak Bay, AK beginning in 2025.

# Required Files (check that script is loading latest version):
# data/seagrass_data.csv

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                         
#                                                                              
# LOAD PACKAGES    
# GPS INTERPOLATION FUNCTION
# READ IN AND PREPARE DATA                                                    
# PROCESS DATA                                                              
#                                                                              
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES ####                                                              
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(purrr)
library(viridis)
library(terra)
library(sf)
library(basemaps)
library(tidyterra)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GPS INTERPOLATION FUNCTION ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


interpolate_quadrats_with_measured <- function(df_group) {
  group_date <- df_group$date[1]
  group_transect <- df_group$transect[1]
  
  # Step 1: Quadrat coordinates
  quadrats <- df_group |>
    distinct(quadrat, lat, lon) |>
    arrange(quadrat) |>
    mutate(q_order = row_number())
  
  # Step 2: Separate measured and interpolated densities
  measured_df <- df_group |>
    filter(densityID == 1) |>
    distinct(date, transect, quadrat, lat, lon, densityID, densityValue)
  
  interpolated_df <- df_group |>
    filter(densityID != 1) |>
    select(quadrat, densityID, densityValue) |>
    arrange(quadrat, densityID) |>
    group_by(quadrat) |>
    mutate(density_index = row_number()) |>
    ungroup()
  
  # Step 3: Interpolated values between quadrats
  interpolated <- map2_dfr(
    1:(nrow(quadrats) - 1),
    2:nrow(quadrats),
    ~{
      q1 <- quadrats[.x, ]
      q2 <- quadrats[.y, ]
      q_start <- q1$quadrat
      q_end <- q2$quadrat
      
      d <- interpolated_df |> filter(quadrat == q_start)
      if (nrow(d) == 0) return(NULL)
      
      tibble(
        date = group_date,
        transect = group_transect,
        quadrat_start = q_start,
        quadrat_end = q_end,
        densityID = d$densityID,
        densityValue = d$densityValue,
        interp_index = d$density_index,
        lat = seq(q1$lat, q2$lat, length.out = nrow(d) + 2)[-c(1, nrow(d) + 2)],
        lon = seq(q1$lon, q2$lon, length.out = nrow(d) + 2)[-c(1, nrow(d) + 2)],
        source = "interpolated"
      )
    }
  )
  
  # Step 4: Add measured points back (always densityID == 1)
  measured <- measured_df |>
    mutate(
      quadrat_start = quadrat,
      quadrat_end = NA,
      interp_index = 1,
      source = "measured"
    ) |>
    select(date, transect, quadrat_start, quadrat_end, densityID, densityValue,
           interp_index, lat, lon, source)
  
  # Step 5: Combine both
  bind_rows(measured, interpolated)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# read in seagrass data table
seagrass_raw <- read_csv("data/seagrass_data.csv", 
                         col_types = cols(date = col_date(format = "%m/%d/%Y"), 
                                          timeBegin = col_time(format = "%H%M"), 
                                          timeEnd = col_time(format = "%H%M")))

# Main pipeline to apply across all date + transect combinations
interpolated_all <- seagrass_raw |>
  group_by(date, transect) |>
  group_split() |>
  map_dfr(interpolate_quadrats_with_measured) 

interp_vect <- vect(interpolated_all, geom=c("lon","lat"), crs = "epsg:4326")
interp_vect <- project(interp_vect, "epsg:3338")
interp_vect$densityValue <- factor(interp_vect$densityValue, levels = c("thick",
                                                                        "moderate",
                                                                        "thin",
                                                                        "sparse",
                                                                        "trace",
                                                                        "none"))

crdref <- "+proj=longlat +datum=WGS84"
interpcrs <- crs(interp_vect)

# make density values as df
interp_latlon <- project(interp_vect, crdref)
interp_df <- data.frame(crds(interp_latlon))
interp_df$type <- "seagrass"
interp_df$qual_dens <- interp_vect$densityValue
interp_df <- data.frame(interp_df)
interp_df <- interp_df |>
  mutate(lat = y, lon = x) |>
  select(-x, -y)

write_csv(interp_df, "data/seagrass_density.csv")

#### END OF SCRIPT ####
