#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                                                             ##
# KBay Seagrass Exploratory                                                   ##
# Script created 2025-06-04                                                   ##
# Data source: NOAA-NCCOS-KBL                                                 ##
# R code prepared by Ross Whippo                                              ##
# Last updated 2025-12-03                                                     ##
#                                                                             ##
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# SUMMARY:
# Script for the exploration of seagrass density, cover, and extent data
# collected from Kachemak Bay, AK beginning in 2025.

# Required Files (check that script is loading latest version):
# data/seagrass_data.csv
# data/MarineGEO_3M_Seagrass_Blades.csv
# data/MarineGEO_3M_Seagrass_Cover.csv
# data/MarineGEO_3M_Seagrass_Density.csv
# data/MarineGEO_3M_Seagrass_Shoots.csv

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
library(tidyterra)

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
     plg=list(x="topright", title="Cover", bty = "o"), main="", axes=FALSE)

# plot with leaflet map
new_interp <- interp_vect |>
  rename(QualitativeCover = densityValue)
  
plet(new_interp, "QualitativeCover",col=viridis(6, option = "mako",
                                             begin = .6, end = 1)) 

interp_df |>
  ggplot() +
  geom_point(aes(x = lon, y = lat, color = qual_dens)) +
  scale_color_brewer(direction = -1, palette = 2) 

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
# Seldovia NAVD88 0 m is 1.573 m above MLLW
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


adjusted <- bathymetry - 1.546
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

# test figs for SBB AMSS Poster


## x = depth, y = density/cover 

###combine datasets so depth and density are in the same df, reproject bath to 
### get lat lon
bath_dd <- project(bathymetry + 1.573, "EPSG:4326")

### Extract raster values at data frame locations
### The extract function will find the raster cell corresponding to each point
### and return its value.
extracted_values <- extract(bath_dd, interpolated_all[, c("lon", "lat")])
seagrass_depths <- cbind(interpolated_all, depth = extracted_values[,2])

seagrass_depths |>
  mutate(densityNum = case_when(densityValue == "none" ~ 0,
                                densityValue == "trace" ~ 1,
                                densityValue == "sparse" ~ 10,
                                densityValue == "thin" ~ 25,
                                densityValue == "moderate" ~ 75,
                                densityValue == "thick" ~ 100)) |> 
ggplot() +
  geom_jitter(aes(x = depth, y = densityNum, color = month(date, label = TRUE)), #shape = source),
              width = 0, height = 5, size = 4, alpha = 0.75) +
  theme_bw() +
  scale_color_viridis_d() +
  labs(x = "Depth (m)", 
       y = "Maximum Cover (%)",
       color = "Month") +
       #shape = "GPS location") +
  geom_smooth(aes(x = depth, y = densityNum),
              method = "lm", se = FALSE) +
  scale_y_continuous(limits = c(-5, 105))

# x = density, y = reproductive
# Add the reproductive binary column
seagrass_depths |>
  mutate(Reproductive = as.numeric(str_detect(notes, "reproductive"))) |>
  mutate(densityNum = case_when(densityValue == "none" ~ 0,
                                densityValue == "trace" ~ 1,
                                densityValue == "sparse" ~ 10,
                                densityValue == "thin" ~ 25,
                                densityValue == "moderate" ~ 75,
                                densityValue == "thick" ~ 100)) |>
  left_join(seagrass_raw |>
              select(notes))
  #filter(densityNum > 0) |>
  ggplot() +
  geom_jitter(aes(x = densityNum, y = Reproductive, color = month(date, label = TRUE)), 
              height = 0.05, width = 4, alpha = 0.75, size = 4) +
  stat_smooth(aes(x = densityNum, y = Reproductive),
              method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  labs(x = "Maximum Cover (%)", 
       y = "Reproductive Plants Detected",
       color = "Month",
       shape = "GPS location") +
  theme_bw()


# INCLUDE IN POSTER  
  # reproductives by cover
seagrass_raw |>
  filter(densityID == 1) |>
  mutate(Reproductive = as.numeric(str_detect(notes, "reproductive"))) |>
  mutate(densityNum = case_when(densityValue == "none" ~ 0,
                                densityValue == "trace" ~ 1,
                                densityValue == "sparse" ~ 10,
                                densityValue == "thin" ~ 25,
                                densityValue == "moderate" ~ 75,
                                densityValue == "thick" ~ 100)) |>
  ggplot() +
  geom_jitter(aes(x = densityNum, y = Reproductive, color = month(date, label = TRUE)), 
              height = 0.05, width = 4, alpha = 0.75, size = 4) +
  stat_smooth(aes(x = densityNum, y = Reproductive),
              method = "glm", method.args = list(family = "binomial"), se = FALSE) +
  labs(x = "Maximum Cover (%)", 
       y = "Reproductive Plants Detected",
       color = "Month",
       shape = "GPS location") +
  theme_bw()
  
 
# INCLUDE IN POSTER
# leaflet of all surveys
plet(new_interp, "QualitativeCover",col=viridis(6, option = "magma",
                                                begin = .6, end = 1)) 

# dryweights
weights <- read_csv("data/MarineGEO_3M_Seagrass_Shoots.csv")
calc_weights <- weights |>
  mutate(epibionts = epibiont_drymass - epibiont_empty) |>
  mutate(blades = blade_drymass - blade_empty) |>
  select(epibionts, blades) |>
  pivot_longer(epibionts:blades, names_to = "category", values_to = "mass (g)") |>
  mutate(`mass (g)` = case_when(`mass (g)` < 0 ~ 0,
                                .default = `mass (g)`))

# INCLUDE IN POSTER
# dryweights
ggplot(calc_weights, aes(x = category, y = `mass (g)`, fill = category)) +
  labs(y = "dry mass (g)") +
  geom_boxplot() +
  theme_bw() +
  scale_fill_viridis(discrete = TRUE,
                     begin = 0.9, 
                     end = 0.4)

# 


# Check depth accuracy of mud bay
# 2. Define the classification breaks and new values
# The matrix should have 3 columns: from value, to value, new value
# This example creates three bins: 0-30, 30-70, 70-100
reclass_matrix <- matrix(
  c(-4, -3, 1,
    -3, -2, 2,
    -2, -1, 3,
    -1, 0, 4,
    0, 1, 5,
    1, 2, 6),
  ncol = 3,
  byrow = TRUE
)

# 3. Apply the classification
r_binned <- classify(bath_clamp, reclass_matrix, include.lowest=TRUE)
# include.lowest=TRUE ensures the lowest value (0) is included in the first bin.

# 4. View results (optional)
plot(r_binned)
freq(r_binned) # Check the frequency of each bin

# PROBLEM : the bathymetry product seems to have a break at 0 m. This may
# be a result of the way the actual bathymetry product is stitched together.
# The subtidal portions may not have been corrected for the diff between 
# NAVD88 and MLLW in Homer/Seldovia. 








# Create example data
# Create example data
df <- data.frame(
  x = 1:20,
  y = 2 * (1:20)^2 - 5 * (1:20) + 10 + rnorm(20, sd = 20)
)

# Create a scatter plot with a parabolic smooth line
ggplot(df, aes(x = x, y = y)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = TRUE)




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


writeVector(interp_vect, "C:/Users/Ross.Whippo/Desktop/seagrass.kml", filetype = "KML")



























































# Vectorized version
mandelbrot_vectorized <- function(xmin=-2, xmax=2, nx=500,
                                  ymin=-1.5, ymax=1.5, ny=500,
                                  n=100, showplot=TRUE,
                                  cols=colorRampPalette(c("blue","yellow","red","black"))(11)) 
{
  
  # variables
  x <- seq(xmin, xmax, length.out=nx)
  y <- seq(ymin, ymax, length.out=ny)
  c <- outer(x,y*1i,FUN="+")
  z <- matrix(0.0, nrow=length(x), ncol=length(y))
  k <- matrix(0.0, nrow=length(x), ncol=length(y))
  
  for (rep in 1:n) { 
    index <- which(Mod(z) < 2)
    z[index] <- z[index]^2 + c[index]
    k[index] <- k[index] + 1
  }
  
  if (showplot==TRUE) { image(x,y,k,col=cols, xlab="Re(c)", ylab="Im(c)")}
  
  return(k)
  
}
mandelbrot_vectorized(n = 2000)
