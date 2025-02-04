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

## Set working directory (any desired working directory for writing/reading files)

directory <- "C:/.../" #eg. "C:/Users/username/WA-Endemic-Modeling/"
setwd(directory)

## Create SpatVectors and Shapefiles (optional) for Study Extent

WenatcheeExtentPolygon <- data.frame(x = c(-120.600, -120.600, -121.135, -121.135),
                                     y = c(47.635, 47.340, 47.340, 47.635)) %>% 
  Polygon() %>% 
  list() %>% 
  Polygons(ID = "Wenatchee Study Extent") %>% 
  list() %>% 
  SpatialPolygons(proj4string = CRS("WGS84"))

SpatVectorWenatchee <- vect(WenatcheeExtentPolygon) #Convert polygon to SpatVector class for terra package

shapefile(x = WenatcheeExtentPolygon, file = paste0(directory, "WenatcheeExtent.shp")) #Optional-Write as shapefile


RainierExtentPolygon <- data.frame(x = c(-121.955, -121.955, -121.540, -121.540),
                                   y = c(46.990, 46.735, 46.735, 46.990)) %>% 
  Polygon() %>% 
  list() %>% 
  Polygons(ID = "Rainier Study Extent") %>% 
  list() %>% 
  SpatialPolygons(proj4string = CRS("WGS84"))

SpatVectorRainier <- vect(RainierExtentPolygon) #Convert polygon to SpatVector class for terra package

shapefile(x = RainierExtentPolygon, file = paste0(directory, "RainierExtent.shp")) #Optional-Write as shapefile

## Digital Elevation Models

WenatcheeW <- rast(paste0(directory, "WenatcheeDEM_W_GeoTiff_USGS_One_Third_Arc_Second_n48w122_20230307.tif"))
WenatcheeE <- rast(paste0(directory, "WenatcheeDEM_E_GeoTiff_USGS_One_Third_Arc_Second_n48w121_20240617.tif"))
Rainier <- rast(paste0(directory, "RainierDEM_GeoTiff_USGS_One_Third_Arc_Second_n47w122_20220919.tif"))

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
ClimateNAVars("Rainier", "8GCMs_ensemble_ssp245_2071-2100.gcm")
ClimateNAVars("Wenatchee", "Normal_1961_1990.nrm")
ClimateNAVars("Wenatchee", "8GCMs_ensemble_ssp245_2071-2100.gcm")


###Reads .tif files from ClimateNAVars(). Creates lists of RasterStack objects for PPT, Tmin, Tmax.
MakeRasterStacks <- function(L,f,g) {
  
  PrecStacks <- list()
  TminStacks <- list()
  TmaxStacks <- list()
  
  for (i in 1:15) {
    PrecStacks[[length(PrecStacks) + 1]] <- c(paste0(directory, L, "DEMTile", as.character(i), "/", g, "/", c(paste0(rep("PPT0", 9), seq(1, 9)), "PPT10", "PPT11", "PPT12"), ".tif"))
    TminStacks[[length(TminStacks) + 1]] <- c(paste0(directory, L, "DEMTile", as.character(i), "/", g, "/", c(paste0(rep("Tmin0", 9), seq(1, 9)), "Tmin10", "Tmin11", "Tmin12"), ".tif"))
    TmaxStacks[[length(TmaxStacks) + 1]] <- c(paste0(directory, L, "DEMTile", as.character(i), "/", g, "/", c(paste0(rep("Tmax0", 9), seq(1, 9)), "Tmax10", "Tmax11", "Tmax12"), ".tif"))
  }
  
  assign(paste0(f, "_PrecStacks"), lapply(PrecStacks, raster::stack), envir = .GlobalEnv)
  assign(paste0(f, "_TminStacks"), lapply(TminStacks, raster::stack), envir = .GlobalEnv)
  assign(paste0(f, "_TmaxStacks"), lapply(TmaxStacks, raster::stack), envir = .GlobalEnv)
}

MakeRasterStacks("Rainier", "R_1961_1990", "Normal_1961_1990")

#Check folder name in directory/RainierTile1/ for third argument 
MakeRasterStacks("Rainier", "R_2071_2100", "8GCMs_ensemble_ssp245_2071-2100")

MakeRasterStacks("Wenatchee", "W_1961_1990", "Normal_1961_1990")

#Check folder name in directory/WenatcheeTile1/ for third argument 
MakeRasterStacks("Wenatchee", "W_2071_2100", "8GCMs_ensemble_ssp245_2071-2100")

###NOTE:Long run time (30min+). Passes RasterStacks through biovars() to generate 19 Bioclimatic variables. Outputs list of Bioclimatic data tiles. 
GenerateBioclimaticVariables <- function(h,j,k,m) {
  
  BiovarsList <- list()
  
  for (i in 1:15) {
    
    BiovarsList[[length(BiovarsList) + 1]] <- biovars(j[[i]],k[[i]],m[[i]])
  }
  
  assign(paste0(h, "_BiovarsList"), lapply(BiovarsList, terra::rast), envir = .GlobalEnv)
}

GenerateBioclimaticVariables("Rainier_1961_1990", R_1961_1990_PrecStacks, R_1961_1990_TminStacks, R_1961_1990_TmaxStacks)
GenerateBioclimaticVariables("Rainier_2071_2100", R_2071_2100_PrecStacks, R_2071_2100_TminStacks, R_2071_2100_TmaxStacks)
GenerateBioclimaticVariables("Wenatchee_1961_1990", W_1961_1990_PrecStacks, W_1961_1990_TminStacks, W_1961_1990_TmaxStacks)
GenerateBioclimaticVariables("Wenatchee_2071_2100", W_2071_2100_PrecStacks, W_2071_2100_TminStacks, W_2071_2100_TmaxStacks)

###Creates mosaic of full study area after resampling to standardize resolution of bioclimatic variable tiles.
###Output will be SpatRaster object with the name of the input, minus "List".
BiovarsMosaic <- function(p,q) {
  Biovars <- list()
  for (i in 1:15) {
    Biovars[[length(Biovars) + 1]] <- terra::resample(p[[i]], q, method = "bilinear")
  }
  
  name <- as.character(substitute(p))
  BiovarsSprc <- sprc(Biovars)
  assign(substr(name, 1, nchar(name) - 4), terra::mosaic(BiovarsSprc), envir = .GlobalEnv)
  
}

BiovarsMosaic(Rainier_1961_1990_BiovarsList, RainierDEM)
BiovarsMosaic(Rainier_2071_2100_BiovarsList, RainierDEM)
BiovarsMosaic(Wenatchee_1961_1990_BiovarsList, WenatcheeDEM)
BiovarsMosaic(Wenatchee_2071_2100_BiovarsList, WenatcheeDEM)

###Write to a .tif file-- Make sure to name this how you want
terra::writeRaster(Rainier_1961_1990_Biovars, "Rainier_1961_1990_Biovars.tif")
terra::writeRaster(Rainier_2071_2100_Biovars, "Rainier_2071_2100_Biovars.tif")
terra::writeRaster(Wenatchee_1961_1990_Biovars, "Wenatchee_1961_1990_Biovars.tif")
terra::writeRaster(Wenatchee_2071_2100_Biovars, "Wenatchee_2071_2100_Biovars.tif")

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

# 2.2 Adjusting Predictor Variables ----

# 2.3 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
