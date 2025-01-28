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
install.packages("raster")
install.packages("dismo")
install.packages("stats")

library(biomod2)
library(tidyverse)
library(rinat)
#library(rgbif)
#library(rvest)
library(ClimateNAr)
library(terra)
library(sp)
library(raster)
library(dismo)
library(stats)

# 1.2 Load Functions ----

# 1.3 Load Data ----

## Set working directory

directory <- "C:/.../"

## Create Shapefiles and SpatVectors for Study Extent

WenatcheeExtentCoords <- data.frame(x = c(-120.670, -120.670, -121.025, -121.025),
                                    y = c(47.565, 47.340, 47.340, 47.565))

WenatcheeExtentPolygon <- Polygon(WenatcheeExtentCoords) %>% 
  list() %>% 
  Polygons(ID = "Wenatchee Study Extent") %>% 
  list() %>% 
  SpatialPolygons(proj4string = CRS("WGS84"))

shapefile(x = WenatcheeExtentPolygon, file = paste0(directory, "WenatcheeExtent.shp"))

SpatVectorWenatchee <- vect(WenatcheeExtentPolygon)


RainierExtentCoords <- data.frame(x = c(-121.955, -121.955, -121.540, -121.540),
                                    y = c(46.990, 46.735, 46.735, 46.990))

RainierExtentPolygon <- Polygon(RainierExtentCoords) %>% 
  list() %>% 
  Polygons(ID = "Rainier Study Extent") %>% 
  list() %>% 
  SpatialPolygons(proj4string = CRS("WGS84"))

shapefile(x = RainierExtentPolygon, file = paste0(directory, "RainierExtent.shp"))

SpatVectorRainier <- vect(RainierExtentPolygon)

## Digital Elevation Models

WenatcheeW <- rast(paste0(directory, "WenatcheeDEM W GeoTiff USGS One Third Arc Second n48w122 20230307.tif"))
WenatcheeE <- rast(paste0(directory, "WenatcheeDEM E GeoTiff USGS One Third Arc Second n48w121 20240617.tif"))
Rainier <- rast(paste0(directory, "RainierDEM GeoTiff USGS One Third Arc Second n47w122 20220919.tif"))

WenatcheeDEM <- mosaic(WenatcheeW, WenatcheeE) %>% 
  crop(SpatVectorWenatchee) %>% 
  project("+proj=longlat +datum=WGS84")

RainierDEM <- crop(Rainier, SpatVectorRainier) %>% 
  project("+proj=longlat +datum=WGS84")

writeRaster(WenatcheeDEM, filename = paste0(directory, "WenatcheeDEM.tif"))

writeRaster(RainierDEM, filename = paste0(directory, "RainierDEM.tif"))

## Climate Data

###Function for making 15 tiles in each study area. Smaller input area needed for climateNAr
ExtentTiles <- function(a,b,c) {
  grid <- rast(a, nrows = 3, ncols = 5)
  makeTiles(b, grid, filename = paste0(c, "DEMTile.tif"))
}

ExtentTiles(SpatVectorRainier, RainierDEM, "Rainier")
ExtentTiles(SpatVectorWenatchee, WenatcheeDEM, "Wenatchee")

###NOTE:Long run time (1hr+) Downscales climate data to ~10m resolution. Outputs folders containing .tifs of PPT01-12, Tmin01-12, Tmax01-12.
ClimateNAVars <- function(d,e) {
  for (i in 1:15) {
    climateNAr(inputFile = paste0(directory, d, "DEMTile", as.character(i), ".tif"),
               periodList = e,
               varList = c(paste0(rep("PPT0", 9), seq(1, 9)), "PPT10", "PPT11", "PPT12",
                           paste0(rep("Tmin0", 9), seq(1, 9)), "Tmin10", "Tmin11", "Tmin12",
                           paste0(rep("Tmax0", 9), seq(1, 9)), "Tmax10", "Tmax11", "Tmax12"),
               
               outDir = directory)
  }
}

ClimateNAVars("Rainier", "Normal_1961_1990.nrm")
ClimateNAVars("Wenatchee", "Normal_1961_1990.nrm")


###Reads .tif files from ClimateNAVars(). Creates lists of RasterStack objects for PPT, Tmin, Tmax.
MakeRasterStacks <- function(L,f,g) {
  
  PrecStacks <- list()
  TminStacks <- list()
  TmaxStacks <- list()
  
  for (i in 1:15) {
    PrecStacks[[length(PrecStacks) + 1]] <- c(paste0(directory, L, "Tile", as.character(i), "/", g, "/", c(paste0(rep("PPT0", 9), seq(1, 9)), "PPT10", "PPT11", "PPT12"), ".tif"))
    TminStacks[[length(TminStacks) + 1]] <- c(paste0(directory, L, "Tile", as.character(i), "/", g, "/", c(paste0(rep("Tmin0", 9), seq(1, 9)), "Tmin10", "Tmin11", "Tmin12"), ".tif"))
    TmaxStacks[[length(TmaxStacks) + 1]] <- c(paste0(directory, L, "Tile", as.character(i), "/", g, "/", c(paste0(rep("Tmax0", 9), seq(1, 9)), "Tmax10", "Tmax11", "Tmax12"), ".tif"))
  }
  
  assign(paste0(f, "_PrecStacks"), lapply(PrecStacks, raster::stack), envir = .GlobalEnv)
  assign(paste0(f, "_TminStacks"), lapply(TminStacks, raster::stack), envir = .GlobalEnv)
  assign(paste0(f, "_TmaxStacks"), lapply(TmaxStacks, raster::stack), envir = .GlobalEnv)
}

MakeRasterStacks("Rainier", "R_1961_90", "Normal_1961_1990")
MakeRasterStacks("Wenatchee", "W_1961_90", "Normal_1961_1990")

###NOTE:Long run time (30min+). Passes RasterStacks through biovars() to generate 19 Bioclimatic variables. Outputs list of Bioclimatic data tiles. 
GenerateBioclimaticVariables <- function(h,j,k,m) {
  
  BiovarsList <- list()
  
  for (i in 1:15) {
    
    BiovarsList[[length(BiovarsList) + 1]] <- biovars(j[[i]],k[[i]],m[[i]])
  }
  
  assign(paste0(h, "_BiovarsList"), lapply(BiovarsList, terra::rast), envir = .GlobalEnv)
}

GenerateBioclimaticVariables("Rainier_1961_90", R_1961_90_PrecStacks, R_1961_90_TminStacks, R_1961_90_TmaxStacks)
GenerateBioclimaticVariables("Wenatchee_1961_90", W_1961_90_PrecStacks, W_1961_90_TminStacks, W_1961_90_TmaxStacks)

###Creates mosaic of full study area after resampling to standardize resolution of bioclimatic variable tiles.
###Output will be SpatRaster object with the name of the input, minus "List".
BiovarsMosaic <- function(p) {
  reference <- list()
  
  for (i in 1:15) {
    reference[[length(reference) + 1]] <- rast(paste0(directory, sub("_.*", "", as.character(substitute(p))), "DEMTile", as.character(i), ".tif"))
  }
  
  Biovars <- list()
  for (i in 1:15) {
    Biovars[[length(Biovars) + 1]] <- terra::resample(p[[i]], reference[[i]], method = "bilinear")
  }
  
  name <- as.character(substitute(p))
  BiovarsSprc <- sprc(Biovars)
  assign(substr(name, 1, nchar(name) - 4), terra::mosaic(BiovarsSprc), envir = .GlobalEnv)
  
}

BiovarsMosaic(Rainier_1961_90_BiovarsList)
BiovarsMosaic(Wenatchee_1961_90_BiovarsList)

###Write to a .tif file-- Make sure to name this how you want
terra::writeRaster(Rainier_1961_90_Biovars, "Rainier_1961_90_Biovars.tif")

## Response Data (Species Occurrences)

# pulling iNaturalist observations for each species into one table using the 'rinat' package
inat.all <- rbind(get_inat_obs(taxon_name = "Pedicularis rainierensis"),
                  get_inat_obs(taxon_name = "Castilleja cryptantha"),
                  get_inat_obs(taxon_name = "Tauschia stricklandii"),
                  get_inat_obs(taxon_name = "Chaenactis thompsonii"),
                  get_inat_obs(taxon_name = "Oreocarya thompsonii"),
                  get_inat_obs(taxon_name = "Lomatium cuspidatum")) %>%
  filter(quality_grade == "research") %>%
  filter(positional_accuracy < 30) %>%
  filter(!is.na(positional_accuracy)) %>%
  filter(coordinates_obscured == "false") #applying filters to ensure identification and spatial accuracy

# creating a file for reference in the repository, may be replaced when more response data is appended
write_csv(inat.all, file = "Data/inat-response-data.csv")


# 2.0 Data Adjustments -----------------------------------------------------

# 2.1 Filtering Occurrence Data ----

library(terra)
# 2.2 Adjusting Predictor Variables ----

# 2.3 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
