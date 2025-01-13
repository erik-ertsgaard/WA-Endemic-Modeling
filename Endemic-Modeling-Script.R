##### WASHINGTON ENDEMIC MODELING PROJECT #####
#### Gjording et al.

# 1.0 LOAD ITEMS -----------------------------------------------------------

# 1.1 Load Packages ----

install.packages("biomod2")
install.packages("tidyverse")
install.packages("rinat")
#install.packages("rgbif")
#install.packages("rvest")
#install ClimateNAr from register.climatena.ca
install.packages("terra")
install.packages("sp")

library(biomod2)
library(tidyverse)
library(rinat)
#library(rgbif)
#library(rvest)
library(ClimateNAr)
library(terra)
library(sp)

# 1.2 Load Functions ----

# 1.3 Load Data ----

# Create Shapefiles and SpatVectors for Study Extent

WenatcheeExtentCoords <- data.frame(x = c(-120.670, -120.670, -121.025, -121.025),
                                    y = c(47.565, 47.340, 47.565, 47.340))

WenatcheeExtentPolygon <- Polygon(WenatcheeExtentCoords) %>% 
  list() %>% 
  Polygons(ID = "Wenatchee Study Extent") %>% 
  list() %>% 
  SpatialPolygons(proj4string = CRS("WGS84"))

directory <- "C:/..."

shapefile(x = WenatcheeExtentPolygon, file = paste0(directory, "/WenatcheeExtent.shp"))

SpatVectorWenatchee <- vect(WenatcheeExtentPolygon)



RainierExtentCoords <- data.frame(x = c(-121.955, -121.955, -121.540, -121.540),
                                    y = c(46.990, 46.735, 46.990, 46.735))

RainierExtentPolygon <- Polygon(RainierExtentCoords) %>% 
  list() %>% 
  Polygons(ID = "Rainier Study Extent") %>% 
  list() %>% 
  SpatialPolygons(proj4string = CRS("WGS84"))

directory <- "C:/..."

shapefile(x = RainierExtentPolygon, file = paste0(directory, "/RainierExtent.shp"))

SpatVectorRainier <- vect(RainierExtentPolygon)



# Digital Elevation Models

WenatcheeW <- rast("C:/Users/nicgj/Downloads/WenatcheeDEM W GeoTiff USGS One Third Arc Second n48w122 20230307.tif")
WenatcheeE <- rast("C:/Users/nicgj/Downloads/WenatcheeDEM E GeoTiff USGS One Third Arc Second n48w121 20240617.tif")
Rainier <- rast("C:/Users/nicgj/Downloads/RainierDEM GeoTiff USGS One Third Arc Second n47w122 20220919.tif")

Wenatchee <- mosaic(WenatcheeW, WenatcheeE) %>% 
  crop(SpatVectorWenatchee) %>% 
  project("+proj=longlat +datum=WGS84")

Rainier <- crop(Rainier, SpatVectorRainier) %>% 
  project("+proj=longlat +datum=WGS84")

directory <- "C:/Users/nicgj/OneDrive/Documents/50 Peaks Plus/Wenatchee Climate Modeling/2025Work/"

writeRaster(Wenatchee, filename = paste0(directory, "WenatcheeDEM.tif"))

writeRaster(Rainier, filename = paste0(directory, "RainierDEM.tif"))


# Climate Data

climateNAr(inputFile = paste0(directory, "WenatcheeDEM.tif"),
           periodList = "Normal_1961_1990.nrm",
           varList = "YSM",
           outDir = directory)

# Response Data - Occurrences

inat.all <- rbind(get_inat_obs(taxon_name = "Pedicularis rainierensis"),
                  get_inat_obs(taxon_name = "Castilleja cryptantha"),
                  get_inat_obs(taxon_name = "Tauschia stricklandii"),
                  get_inat_obs(taxon_name = "Chaenactis thompsonii"),
                  get_inat_obs(taxon_name = "Oreocarya thompsonii"),
                  get_inat_obs(taxon_name = "Lomatium cuspidatum")) %>%
  filter(quality_grade == "research") %>%
  filter(positional_accuracy < 30) %>%
  filter(!is.na(positional_accuracy)) %>%
  filter(coordinates_obscured == "false")

pera.inat <- get_inat_obs(taxon_name = "Pedicularis rainierensis",
                          quality = "research",
                          geo = TRUE) %>%
  filter(positional_accuracy < 30) %>%
  filter(coordinates_obscured == "false")

cacr.inat <- get_inat_obs(taxon_name = "Castilleja cryptantha",
                          quality = "research",
                          geo = TRUE) %>%
  filter(positional_accuracy < 30) %>%
  filter(coordinates_obscured == "false")

tast.inat <- get_inat_obs(taxon_name = "Tauschia stricklandii",
                          quality = "research",
                          geo = TRUE) %>%
  filter(positional_accuracy < 30) %>%
  filter(coordinates_obscured == "false")

cath.inat <- get_inat_obs(taxon_name = "Chaenactis thompsonii",
                          quality = "research",
                          geo = TRUE) %>%
  filter(positional_accuracy < 30) %>%
  filter(coordinates_obscured == "false")

orth.inat <- get_inat_obs(taxon_name = "Oreocarya thompsonii",
                          quality = "research",
                          geo = TRUE) %>%
  filter(positional_accuracy < 30) %>%
  filter(coordinates_obscured == "false")

locu.inat <- get_inat_obs(taxon_name = "Lomatium cuspidatum",
                          quality = "research",
                          geo = TRUE) %>%
  filter(positional_accuracy < 30) %>%
  filter(coordinates_obscured == "false")

# 2.0 Data Adjustments -----------------------------------------------------

# 2.1 Filtering Occurrence Data ----

library(terra)
# 2.2 Adjusting Predictor Variables ----

# 2.3 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
