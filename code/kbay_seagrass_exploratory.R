#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             ##
# KBay Seagrass Exploratory                                                   ##
# Script created 2025-06-04                                                   ##
# Data source: NOAA-NCCOS-KBL                                                 ##
# R code prepared by Ross Whippo                                              ##
# Last updated 2025-07-15                                                     ##
#                                                                             ##
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:
# Script for the exploration of seagrass density, cover, and extent data
# collected from Kachemak Bay, AK beginning in 2025.

# Required Files (check that script is loading latest version):
# data/seagrass_data.csv

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABLE OF CONTENTS                                                         ####
#                                                                              +
# LOAD PACKAGES                                                                +
# READ IN AND PREPARE DATA                                                     +
# MANIPULATE DATA                                                              +
#                                                                              +
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LOAD PACKAGES                                                             ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

library(tidyverse)
library(purrr)
library(viridis)
library(terra)
library(sf)
library(basemaps)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# GPS INTERPOLATION FUNCTION                                                ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


interpolate_quadrats_with_measured <- function(df_group) {
  group_date <- df_group$date[1]
  group_transect <- df_group$transect[1]
  
  # Step 1: Quadrat coordinates
  quadrats <- df_group %>%
    distinct(quadrat, lat, lon) %>%
    arrange(quadrat) %>%
    mutate(q_order = row_number())
  
  # Step 2: Separate measured and interpolated densities
  measured_df <- df_group %>%
    filter(densityID == 1) %>%
    distinct(date, transect, quadrat, lat, lon, densityID, densityValue)
  
  interpolated_df <- df_group %>%
    filter(densityID != 1) %>%
    select(quadrat, densityID, densityValue) %>%
    arrange(quadrat, densityID) %>%
    group_by(quadrat) %>%
    mutate(density_index = row_number()) %>%
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
      
      d <- interpolated_df %>% filter(quadrat == q_start)
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
  measured <- measured_df %>%
    mutate(
      quadrat_start = quadrat,
      quadrat_end = NA,
      interp_index = 1,
      source = "measured"
    ) %>%
    select(date, transect, quadrat_start, quadrat_end, densityID, densityValue,
           interp_index, lat, lon, source)
  
  # Step 5: Combine both
  bind_rows(measured, interpolated)
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# READ IN AND PREPARE DATA                                                  ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# read in seagrass data table
seagrass_raw <- read_csv("data/seagrass_data.csv", 
                         col_types = cols(date = col_date(format = "%m/%d/%Y"), 
                                          timeBegin = col_time(format = "%H%M"), 
                                          timeEnd = col_time(format = "%H%M")))

# turn points into spatial object
seagrass_shape <- vect(seagrass_raw)
# visual check 
plot(seagrass_shape)


# Main pipeline to apply across all date + transect combinations
interpolated_all <- seagrass_raw %>%
  group_by(date, transect) %>%
  group_split() %>%
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

# export density values as df
interp_latlon <- project(interp_vect, crdref)
interp_df <- data.frame(crds(interp_latlon))
interp_df$type <- "seagrass"
interp_df$qual_dens <- interp_vect$densityValue
interp_df <- data.frame(interp_df)
interp_df <- interp_df |>
  mutate(lat = y, lon = x) |>
  select(-x, -y)

write_csv(interp_df, "data/seagrass_density.csv")

# import west mud bay meadow polygon
meadow_lines <- vect("data/Seagrasspathpolygon.kml")
meadow_poly <- as.polygons(meadow_lines)
meadow_poly <- project(meadow_poly, interpcrs)
expanse(meadow_poly)

plot(interp_dense)
plot(meadow_poly, add = TRUE)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# VISUALIZE PLOTS                                                           ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# plot, no baselayer
plot(interp_vect, "densityValue", col=viridis(6, option = "mako"), sort = c("thick",
                                                                            "moderate",
                                                                            "thin",
                                                                            "sparse",
                                                                            "trace",
                                                                            "none"),
     plg=list(x="topright", title="Density", bty = "o"), main="", axes=FALSE)

# plot with leaflet map
plet(interp_vect, "densityValue",col=viridis(6, option = "mako"))

interp_df |>
  ggplot() +
  geom_point(aes(x = lon, y = lat, color = qual_dens))

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# INTERPOLATE DENSITY                                                       ####
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# read in bathymetry layer
bathymetry <- terra::rast("../kbay_SAV-HSI_model/data/b_intermediate_data/bathymetry/bathymetry.grd")

# read in ROI
ROI <- vect("../kbay_SAV-HSI_model/data/b_intermediate_data/roi")

# read in mud bay
mudbay <- vect("data/homerSpitEast.kml")
mudbay <- project(mudbay, "epsg:3338")


# clip bathymetry layer to ROI and mud bay
bath_mask <- mask(bathymetry, ROI)
bath_mask <- mask(bath_mask, mudbay)

# correct depth for NAVD88/MLLW offset
# Alaska Tidal Datum Portal. 
# Alaska Tidal Datum Reference Table. 
# Retrieved from 
# https://dggs.alaska.gov/hazards/coastal/ak-tidal-datum-portal.html  
bath_corr <- bath_mask + 1.573
bath_clamp <- clamp(bath_corr, lower = -6, upper = 3, values = FALSE)
bath_clamp <- interpIDW(bath_clamp,
                        ROI,
                        field = "KBL-bathymetry_GWA-area_50m_EPSG3338",
                        radium = 50)
bath_clamp <- crop(bath_clamp, mudbay)
plet(bath_clamp)

# assign values to qualitative densities
densityValue <- interp_vect$densityValue
characterDens <- as.character(densityValue)
densDF <- data.frame(characterDens)
numericDens <- densDF %>%
  mutate(densityNum = case_when(characterDens == "none" ~ 0,
                                characterDens == "trace" ~ 1,
                                characterDens == "sparse" ~ 2,
                                characterDens == "thin" ~ 3,
                                characterDens == "moderate" ~ 4,
                                characterDens == "thick" ~ 5)) %>%
  select(densityNum)
interp_vect$densityNum <- numericDens$densityNum

# interpolated depth for missing bathymetry
plot(bath_clamp)


# interpolated density
interp_dense <- interpIDW(bath_clamp, 
                          interp_vect, 
                          field = "densityNum", 
                          radius = 500,
                          smooth = 100)
interp_mask <- mask(interp_dense, bath_clamp)


plet(interp_mask)
plet(interp_vect)

plot(bath_clamp)
plot(interp_vect, add = TRUE)

plot(test, col = viridis(12, option = "magma"))
plot(interp_mask, add = TRUE)


plot(interp_vect, "densityValue", col=viridis(6, option = "viridis", begin = 1, end = 0), sort = c("thick",
                                                                            "moderate",
                                                                            "thin",
                                                                            "sparse",
                                                                            "trace",
                                                                            "none"),
     plg=list(x="topright", title="Density", bty = "o"), main="", axes=FALSE,
     add = TRUE)


adjusted <- bathymetry + 1.546
test <- crop(adjusted, mudbay)
plot(test)

test <- crop(bathymetry, mudbay)
test <- clamp(test, lower = -6, upper = 3, values = TRUE)
plot(bath_clamp) 
plot(bath_clamp, add = TRUE)
plot(interp_vect, add = TRUE)
############### SUBSECTION HERE
####
#<<<<<<<<<<<<<<<<<<<<<<<<<<END OF SCRIPT>>>>>>>>>>>>>>>>>>>>>>>>#

# SCRATCH PAD ####


# rasterize plots
plot_rast <- rasterize(interp_vect, bath_corr, field = "densityValue", fun = "mean")


# writeVector(interp_vect, "seagrass_points.kmz", filetype = "KML")


test <- vect("data/mudBay.kml")

r <- rast(xmin = xmin(interp_vect),
          xmax = xmax(interp_vect),
          ymin = ymin(interp_vect),
          ymax = ymax(interp_vect),
          ncols = 200,
          nrows = 200)
test <- rasterize(interp_vect, r, fields = "densityValue")
plot(test)














polyseagrass <- vect("C:/Users/Ross.Whippo/Desktop/Seagrasspathpolygon.kml")
expanse(polyseagrass)
mudbay <- project(mudbay, "epsg:3338")



polyseagrass




# TEST FOR INTERPOLATION

ordinal_map <- c("none" = 1, 
                 "trace" = 2, 
                 "sparse" = 3,
                 "thin" = 4,
                 "moderate" = 5,
                 "thick" = 6)

interp_vect$ordinal <- ordinal_map[interp_vect$denseValue]
crs(interp_vect) <- "EPSG:4326"


r_template <- rast(ext(ROI))

r_points <- rasterize(interp_vect, r_template, field = "ordinal", fun = "first")

dist_ras <- distance(interp_vect)

inv_dist <- 1 / (dist_ras + 1e-6)  # add small number to avoid divide-by-zero

interp <- interpolate(r_template, interp_vect, "ordinal", method = "idw", idp = 2)

