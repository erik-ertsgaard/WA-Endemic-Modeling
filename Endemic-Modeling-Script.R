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
install.packages("sf")
install.packages("rnaturalearth")
install.packages("ggspatial")
install.packages("prettymapr")
install.packages("vegan")
install.packages("hexbin")
install.packages("ggnewscale")
install.packages("irlba")
install.packages("ggrepel")

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
library(sf)
library(rnaturalearth)
library(ggspatial)
library(prettymapr)
library(vegan)
library(hexbin)
library(ggnewscale)
library(irlba)
library(ggrepel)

# 1.2 Load Functions ----

reproject_to_wgs84 <- function(df) {
  
  df$Geodetic.Datum <- gsub(" ", "", df$Geodetic.Datum)
  
  df <- df %>%
    mutate(is_wgs84 = ifelse(`Geodetic.Datum` == "WGS84", TRUE, FALSE))
  
  data_wgs84 <- df %>% filter(is_wgs84 == TRUE)
  data_not_wgs84 <- df %>% filter(is_wgs84 == FALSE)
  
  # Handle reprojection based on the original CRS for NAD83 and NAD27
  if (nrow(data_not_wgs84) > 0) {
    # Check the existing CRS and reproject accordingly
    data_not_wgs84 <- st_transform(data_not_wgs84, crs = st_crs(4326))  # Default to WGS 84 (EPSG:4326)
    
    # Handle specific cases for NAD83 and NAD27
    # NAD83 is commonly represented as EPSG:4269
    # NAD27 is commonly represented as EPSG:4267
    
    # Detect if the CRS is NAD83 (EPSG:4269) or NAD27 (EPSG:4267) and reproject accordingly
    data_not_wgs84 <- data_not_wgs84 %>%
      mutate(
        is_nad83 = st_crs(.) == st_crs(4269),
        is_nad27 = st_crs(.) == st_crs(4267)
      ) %>%
      # Reproject to WGS 84 from NAD83 or NAD27
      mutate(
        geometry = case_when(
          is_nad83 ~ st_transform(geometry, st_crs(4326)),  # Reproject NAD83 to WGS 84
          is_nad27 ~ st_transform(geometry, st_crs(4326)),  # Reproject NAD27 to WGS 84
          TRUE ~ geometry  # Keep unchanged if it's already in WGS84
        )
      ) %>%
      # Remove the helper columns
      select(-is_nad83, -is_nad27)
    
    df <- bind_rows(data_wgs84, data_not_wgs84)
  }
  
  return(df)
}

remove_nearby_points <- function(df, dist_threshold = 10) {
  
  # Create an empty list to store the indices of points to keep
  keep_indices <- vector("list", length = nrow(df))
  
  # Loop over all points and compare distances
  for (i in 1:nrow(df)) {
    # Calculate distance between the i-th point and all other points
    distances <- st_distance(df[i, ], df)
    
    # Convert the distances to numeric values (in meters)
    distances_numeric <- as.numeric(distances)
    
    # Check if any point is within the distance threshold (excluding the point itself)
    close_points <- which(distances_numeric < dist_threshold & distances_numeric > 0)
    
    # If no close points, keep this point (it will be added later)
    if (length(close_points) == 0) {
      keep_indices[[i]] <- i
    }
  }
  
  # Return the subset of points to keep
  return(df[unlist(keep_indices), ])
}

# 1.3 Load Data ----

# 1.3a Load Predictor Data ----

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
MakeRasterStacks("Rainier", "R_2071_2100", "8GCMs_ensemble_ssp245_2071-2100")
MakeRasterStacks("Wenatchee", "W_1961_1990", "Normal_1961_1990")
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

## Lithology Data
Lith <- st_read(paste0(directory, "lithologyclip.shp")) %>% 
  vect() %>% 
  project("+proj=longlat +datum=WGS84")

categories <- unique(Lith$Map_Unit)

Lith$Map_Unit_ID <- as.numeric(factor(Lith$Map_Unit, levels = categories))

R.lith <- rasterize(Lith, RainierDEM, field = "Map_Unit_ID")
names(R.lith) <- "Lithology"
W.lith <- rasterize(Lith, WenatcheeDEM, field = "Map_Unit_ID")
names(W.lith) <- "Lithology"

levels(R.lith) <- data.frame(ID = seq_along(categories), Category = categories)
levels(W.lith) <- data.frame(ID = seq_along(categories), Category = categories)

## Soils Data
R.soils <- st_read(paste0(directory, "Rainier_Soil_Types.shp")) %>% 
  mutate(comp_tx_86 = ifelse(is.na(comp_tx_86), "other", comp_tx_86)) %>%
  vect() %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  rasterize(RainierDEM, field = "comp_tx_86")
names(R.soils) <- "Soil Type"

R.categorical <- c(R.soils, R.lith)

W.soils <- st_read(paste0(directory, "Wenatchees_Soil_Types.shp")) %>% 
  vect() %>% 
  project("+proj=longlat +datum=WGS84")
W.soils$comp_tx_86[2] <- "other"
W.soils$comp_tx_86[11] <- "other"
W.soils$comp_tx_86[12] <- "Andic Haplocryods"
W.soils$comp_tx_86[13] <- "other"
W.soils$comp_tx_86[14] <- "Andic Haplocryods"
W.soils$comp_tx_86[19] <- "Andic Dystrocryepts"
W.soils$comp_tx_86[22] <- "other"
W.soils <- rasterize(W.soils, WenatcheeDEM, field = "comp_tx_86")
names(W.soils) <- "Soil Type"

W.categorical <- c(W.soils, W.lith)
##

## Tree Canopy Cover Data (Current Data)
R.Canopy <- rast(paste0(directory, "USANLCDTreeCanopyCover_Ra.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=5) %>% 
  resample(RainierDEM, method = "bilinear")
names(R.Canopy) <- "Tree Canopy Cover"

W.Canopy <- rast(paste0(directory, "USANLCDTreeCanopyCover_We.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=5) %>% 
  resample(WenatcheeDEM, method="bilinear")
names(W.Canopy) <- "Tree Canopy Cover"

## Tree Canopy Cover Data (GLM Model to predict Future Canopy Cover)


## Topological Data
R.TPI.20 <- rast(paste0(directory, "RainierDEM10_TPI_20m.tif"))
R.TPI.30 <- rast(paste0(directory, "RainierDEM10_TPI_30m.tif"))
R.TPI.50 <- rast(paste0(directory, "RainierDEM10_TPI_50m.tif"))
R.TPI.100 <- rast(paste0(directory, "RainierDEM10_TPI_100m.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
names(R.TPI.100) <- "TPI"
R.TPI.300 <- rast(paste0(directory, "RainierDEM10_TPI_300m.tif"))
R.TPI.500 <- rast(paste0(directory, "RainierDEM10_TPI_500m.tif"))


W.TPI.100 <- rast(paste0(directory, "WenatcheesDEM10_TPI100m.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")
names(W.TPI.100) <- "TPI"


R.Slope <- rast(paste0(directory, "Slope_Rainier.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
W.Slope <- rast(paste0(directory, "Slope_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")

R.Northness <- rast(paste0(directory, "Northness_Rainier.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
names(R.Northness) <- "Northness"

R.Eastness <- rast(paste0(directory, "Eastness_Rainier.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
names(R.Eastness) <- "Eastness"

R.NorthEastness <- rast(paste0(directory, "NorthEastness_Rainier.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
names(R.NorthEastness) <- "NorthEastness"

W.Northness <- rast(paste0(directory, "Northness_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")
names(W.Northness) <- "Northness"

W.Eastness <- rast(paste0(directory, "Eastness_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")
names(W.Eastness) <- "Eastness"

W.NorthEastness <- rast(paste0(directory, "NorthEastness_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")
names(W.NorthEastness) <- "NorthEastness"

names(RainierDEM) <- "Elevation"
names(WenatcheeDEM) <- "Elevation"

##Generate Principal Components for Model Input

R_1961_1990_PCAraster <- c(Rainier_1961_1990_Biovars, RainierDEM, R.Canopy, R.TPI.100, R.Slope, R.Northness, R.Eastness, R.NorthEastness)
R_1961_1990_PCAmatrix <- as.matrix(values(R_1961_1990_PCAraster))

R_2071_2100_PCAraster <- c(Rainier_2071_2100_Biovars, RainierDEM, R.Canopy, R.TPI.100, R.Slope, R.Northness, R.Eastness, R.NorthEastness)
R_2071_2100_PCAmatrix <- as.matrix(values(R_2071_2100_PCAraster))

W_1961_1990_PCAraster <- c(Wenatchee_1961_1990_Biovars, WenatcheeDEM, W.Canopy, W.TPI.100, W.Slope, W.Northness, W.Eastness, W.NorthEastness)
W_1961_1990_PCAmatrix <- as.matrix(values(W_1961_1990_PCAraster))

W_2071_2100_PCAraster <- c(Wenatchee_2071_2100_Biovars, WenatcheeDEM, W.Canopy, W.TPI.100, W.Slope, W.Northness, W.Eastness, W.NorthEastness)
W_2071_2100_PCAmatrix <- as.matrix(values(W_2071_2100_PCAraster))

###PCA
Wenatchee_1961_1990_PCA <- prcomp(na.omit(W_1961_1990_PCAmatrix), scale. = TRUE)
Wenatchee_2071_2100_PCA <- prcomp(na.omit(W_2071_2100_PCAmatrix), scale. = TRUE)
Rainier_1961_1990_PCA <- prcomp(na.omit(R_1961_1990_PCAmatrix), scale. = TRUE)
Rainier_2071_2100_PCA <- prcomp(na.omit(R_2071_2100_PCAmatrix), scale. = TRUE)

saveRDS(Wenatchee_1961_1990_PCA, "Wenatchee_1961_1990_PCA.rds")
saveRDS(Wenatchee_2071_2100_PCA, "Wenatchee_2071_2100_PCA.rds")
saveRDS(Rainier_1961_1990_PCA, "Rainier_1961_1990_PCA.rds")
saveRDS(Rainier_2071_2100_PCA, "Rainier_2071_2100_PCA.rds")

# Load later with:
WenatcheeRDS <- readRDS("Wenatchee_1961_1990_PCA.rds")

WenPCAscores <- scores(Wenatchee_1961_1990_PCA)

WenEnvfit <- envfit(Wenatchee_1961_1990_PCA, Wenatchee_1961_1990_PCAMatrix, WenPCAscores, permutations = 0)

summary(Wenatchee_1961_1990_PCA)

# Project PCA result back onto a raster with a layer for each PC
Wenatchee_1961_1990_PC1.6_scores <- predict(W_1961_1990_PCAraster, Wenatchee_1961_1990_PCA, index = 1:6)

Rainier_1961_1990_PC1.6_scores <- predict(R_1961_1990_PCAraster, Rainier_1961_1990_PCA, index = 1:6)


# Combine First 6 Principal Components, Lithology, and Soils layers into single raster

R.expl.var <- c(Rainier_1961_1990_PC1.6_scores, R.lith, R.soils)
W.expl.var <- c(Wenatchee_1961_1990PC1.6_scores, W.lith, W.soils)

# If exporting files, it is best to write the PCs as a .tif and the categorical variables as a .rds
writeRaster(Rainier_1961_1990_PC1.6_scores, paste0(directory, "Rainier-PCs-1-6.tif"))
writeRaster(Wenatchee_1961_1990PC1.6_scores, paste0(directory, "Wenatchee-PCs-1-6.tif"))


writeRaster(R.expl.var, paste0(directory, "Rainier-expl-var-raster.tif"), datatype = c("FLT4S", "INT1U"), overwrite = TRUE)

###
# 1.3b Load Response Data ----

## Response Data (Species Occurrences)

# pulling iNaturalist observations for each species into one table using the 'rinat' package
inat.all <- rbind(get_inat_obs(taxon_name = "Pedicularis rainierensis"),
                  get_inat_obs(taxon_name = "Castilleja cryptantha"),
                  get_inat_obs(taxon_name = "Tauschia stricklandii"),
                  get_inat_obs(taxon_name = "Chaenactis thompsonii"),
                  get_inat_obs(taxon_name = "Oreocarya thompsonii"),
                  get_inat_obs(taxon_name = "Lomatium cuspidatum"),
                  get_inat_obs(taxon_name = "Androsace nivalis"),
                  get_inat_obs(taxon_name = "Claytonia megarhiza nivalis")) %>% #synonym for Wenatchee Mts disjunct
  filter(quality_grade == "research") %>%
  filter(positional_accuracy < 30) %>%
  filter(!is.na(positional_accuracy)) %>%
  filter(coordinates_obscured == "false") #applying filters to ensure identification and spatial accuracy

# creating a file for reference in the repository, may be replaced when more response data is appended
write_csv(inat.all, file = "Data/inat-response-data.csv")

# converting inat.all to spatial data frame
inat.all <- st_as_sf(inat.all,
                     coords = c("longitude", "latitude"), 
                     crs = 4326) 

# reading in CPNWH occurrence data downloaded (and cleaned) from <https://www.pnwherbaria.org/data/search.php>
cpnwh.all <- read.csv("Data/cpnwh_response-data-cleaned.csv") 

#filtering for records that meet occurrence standards and creating a spatial data frame
cpnwh.filtered <- filter(cpnwh.all, Valid.LatLng == "Y") %>%
  filter(Coordinate.Uncertainty.in.Meters <= 1000) %>%
  filter(Decimal.Latitude > 45) %>%
  st_as_sf(coords = c("Decimal.Longitude", "Decimal.Latitude"), 
           crs = 4326) 

# reprojecting points not already in WGS84 to the project's CRS
cpnwh.wgs84 <- reproject_to_wgs84(cpnwh.filtered) %>%
  rename(was_wgs84 = is_wgs84,
         Original.Geodetic.Datum = Geodetic.Datum)
  
st_crs(cpnwh.wgs84) #WGS84 as expected

# selecting relevant columns and renaming to match iNat
cpnwh.wgs84 <- rename(cpnwh.wgs84,
                      scientific_name = Accepted.Name,
                      uncertainty = Coordinate.Uncertainty.in.Meters) 

cpnwh.wgs84$uncertainty <- as.integer(cpnwh.wgs84$uncertainty)

inat.all <- rename(inat.all, uncertainty = positional_accuracy)

# combining all records into one large dataframe
occurrence_data_raw <- bind_rows(inat.all, cpnwh.wgs84) %>%
  st_as_sf(coords = geometry, 
           crs = 4326)

# clipping points' extent to study area
StudyAreaExtentsMerged <- read_sf("Data/StudyAreaExtentsMerged.shp")

occurrence_data_clipped <- st_intersection(occurrence_data_raw, StudyAreaExtentsMerged)

# correcting names to be consistent with current accepted taxonomy 
unique(occurrence_data_clipped$scientific_name)

occurrence_data_clipped$scientific_name[occurrence_data_clipped$scientific_name == "Douglasia nivalis"] <- "Androsace nivalis"
occurrence_data_clipped$scientific_name[occurrence_data_clipped$scientific_name == "Claytonia megarhiza nivalis"] <- "Claytonia megarhiza"

# removing duplicative occurrences within 1 meter of each other
occurrence_data_cleaned <- occurrence_data_clipped %>%
  group_by(scientific_name) %>%
  group_modify(~ remove_nearby_points(.x, dist_threshold = 5)) %>%
  ungroup() %>%
  st_as_sf()

## preliminary summary and plots of response data
summary_response <- occurrence_data_cleaned %>%
  group_by(scientific_name) %>%
  summarise(
    num_occurrences = n())

ggplot() +
  annotation_map_tile(type = "osm") +
  geom_sf(data = occurrence_data_cleaned, aes(color = scientific_name), size = 1, alpha = 0.7) +
  theme_minimal() +
  theme(
    legend.position = "right",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()) +
  coord_sf(datum = NA) +  
  labs(
    title = "Endemic Species Distributions",
    color = "Species"
  ) +
  facet_wrap(~ scientific_name, ncol = 2)  # Facet by species

# 1.3c Load Background Data ----
## Generating background data from iNaturalist to match bias of response data

# vascular plants in Mount Rainier study extent
rainier.background.data <- get_inat_obs(taxon_name = "Phylum Tracheophyta", #vascular plants
                                        geo = TRUE, #must be georeferenced
                                        quality = "research",
                                        bounds = RainierExtentPolygon,
                                        maxresults = 10000) %>%
  filter(positional_accuracy < 30) %>%
  filter(!is.na(positional_accuracy)) %>%
  filter(coordinates_obscured == "false") %>%
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326) 

# vascular plants in Wenatchee Mts study extent
wenatchees.background.data <- get_inat_obs(taxon_name = "Phylum Tracheophyta", #vascular plants
                                           geo = TRUE, #must be georeferenced
                                           quality = "research",
                                           bounds = WenatcheeExtentPolygon,
                                           maxresults = 10000) %>%
  filter(positional_accuracy < 30) %>%
  filter(!is.na(positional_accuracy)) %>%
  filter(coordinates_obscured == "false") %>%
  st_as_sf(coords = c("longitude", "latitude"), 
           crs = 4326) 

# 2.0 Data Adjustments -----------------------------------------------------

# 2.1 Adjusting Predictor Variables ----

# 2.2 Adjusting Response Variables ----
## creating a data frame of occurrence and background data compatible with 'biomod2'

# converting data to a SpatVector for better compatibility with biomod
occurrence_data_cleaned <- vect(occurrence_data_cleaned)

# adding columns to represent 10 repetitions of background data, with 500 randomly-selected background points each
for (i in 1:10) {
  pseudo_absences <- rep(FALSE, nrow(rainier.background.data))  # starting with all rows as false
  
  pseudo_absences[sample(1:nrow(rainier.background.data), 500)] <- TRUE  # Set 500 random points to TRUE
  
  rainier.background.data[[paste0("PA", i)]] <- pseudo_absences  # Add as new column (e.g., PA_1, PA_2, ...)
}

# repeating for wenatchees data
for (i in 1:10) {
  pseudo_absences <- rep(FALSE, nrow(wenatchees.background.data))  # starting with all rows as false
  
  pseudo_absences[sample(1:nrow(wenatchees.background.data), 500)] <- TRUE  # Set 500 random points to TRUE
  
  wenatchees.background.data[[paste0("PA", i)]] <- pseudo_absences  # Add as new column (e.g., PA_1, PA_2, ...)
}


# inputting data into the biomod package's formatting function
BIOMOD_FormatingData(resp.name = "Pedicularis rainierensis",
                     resp.var = occurrence_data_cleaned[occurrence_data_cleaned$scientific_name == "Pedicularis rainierensis", "geometry"],
                     expl.var = ,
                     PA.strategy = "user.defined",
                     PA.nb.rep = 10,
                     PA.user.table = rainier.background.data[, c("geometry", "PA1", "PA2", "PA3", "PA4", "PA5", "PA6", "PA7", "PA8", "PA9", "PA10")])

?bm_CrossValidation()

# 2.3 Principle Coordinate Analysis (PCA) ----

# Load Occurrence Data
Occurrence <- read_csv(paste0(directory, "all-occurence-data_cleaned_latlon.csv")) %>% 
  mutate(Latitude=as.numeric(Latitude)) %>% 
  mutate(Longitude=as.numeric(Longitude)) %>% 
  dplyr::select(scientific_name, Longitude, Latitude) %>% 
  na.omit()

# Rainier Climate PCA

R_CurrentClimate_matrix <- as.matrix(Rainier_1961_1990_Biovars)
R_FutureClimate_matrix <- as.matrix(Rainier_2071_2100_Biovars)

R_ClimateMatrix <- rbind(R_CurrentClimate_matrix, R_FutureClimate_matrix)

# Create a grouping factor: first raster labeled as "Current", second as "Future (2071-2100)"
R_Time_period_labels <- factor(c(rep("Current", nrow(R_CurrentClimate_matrix)), rep("Future (2071-2100)", nrow(R_FutureClimate_matrix))))

#Run PCA (This function only calculates the first 2 PCs)
R_Climate_PCA <- prcomp_irlba(R_ClimateMatrix, n = 2, center = TRUE, scale. = TRUE)

# Extract first two principal components as a dataframe
R_Climate_PCA_df <- data.frame(
  PC1 = R_Climate_PCA$x[,1],
  PC2 = R_Climate_PCA$x[,2],
  Source = R_Time_period_labels  # Label data source
)

# Occurrence Data
R_Occurrence <- filter(Occurrence, scientific_name == "Pedicularis rainierensis" | 
                                   scientific_name == "Castilleja cryptantha" | 
                                   scientific_name == "Tauschia stricklandii")

R_cell_indices <- as.matrix(R_Occurrence[, c("Longitude", "Latitude")]) %>% # Convert to matrix for cell lookup
  cellFromXY(Rainier_1961_1990_Biovars, .)

R_Occurrence_pcs <- as.data.frame(R_Climate_PCA$x[R_cell_indices, 1:2]) # Extract PCA values (PC1 and PC2) for these cells
R_Occurrence_pcs$scientific_name <- R_Occurrence$scientific_name  # Keep species info


# Plot
RainierClimatePCAplot <- ggplot() +
  # Hexbin for Current Climate (red gradient)
  geom_hex(data = subset(R_Climate_PCA_df, Source == "Current"),  
           aes(x = PC1, y = PC2, fill = after_stat(count)),  
           bins = 50, alpha = 0.5) +
  scale_fill_gradient(name = "Climate Period", low = "lightgray", high = "black", guide = "none") +
  
  # Reset the fill scale before adding Future Climate
  new_scale_fill() +
  
  # Hexbin for Future Climate (gray gradient)
  geom_hex(data = subset(R_Climate_PCA_df, Source == "Future (2071-2100)"),  
           aes(x = PC1, y = PC2, fill = after_stat(count)),  
           bins = 50, alpha = 0.5) +
  scale_fill_gradient(name = "Climate Period", low = "lightpink", high = "darkred", guide = "none") +
  
  # Manual Colors for Taxa
  scale_color_manual(name = "Taxon",
                     values = c("blue", "darkgreen", "orange")) +
  
  # Points for Rainier Occurrence (highlighted taxa)
  geom_point(data = R_Occurrence_pcs,  
             aes(x = PC1, y = PC2, color = scientific_name, shape = scientific_name), size = 3, stroke = 1) +
  scale_shape_manual(values = c(1, 0, 2)) + #add 3 and 8 for plus sign and star
  # #Ellipses for Highlighted Taxa
  stat_ellipse(data = R_Occurrence_pcs,
               aes(x = PC1, y = PC2, color = scientific_name),
               level = 0.95, size = 1.2) +
  # stat_chull(data = R_Occurrence_pcs,
  #            geom = "polygon",
  #            aes(x = PC1, y = PC2, color = scientific_name),
  #            alpha = 0.001,
  #            show.legend = FALSE) +
  # 
  new_scale_color() +
  
  #Ellipses for Current and Future Climate Envelopes
  # stat_ellipse(data = R_Climate_PCA_df,
  #              aes(x = PC1, y = PC2, color = Source),
  #              level = 0.95, size = 1.2) +
  # scale_color_manual(name = "Climate Period",
  #                    values = c("Current" = "black", "Future (2071-2100)" = "darkred")) +
  # stat_chull(data = R_Climate_PCA_df,
  #            geom = "polygon",
  #            aes(x = PC1, y = PC2, color = Source),
  #            alpha = 0.2, 
  #            show.legend = FALSE) + 
  # scale_color_manual(name = "Climate Period",
  #                    values = c("Current" = "black", "Future (2071-2100)" = "darkred")) +
  geom_segment(aes(x = -13, xend = 10, y = -6.5, yend = -6.5), 
               arrow = arrow(type = "closed", ends = "both", length = unit(0.3, "cm")),
               linewidth = 1) +
  geom_segment(aes(x = -15, xend = -15, y = -5, yend = 5), 
               arrow = arrow(type = "closed", ends = "both", length = unit(0.3, "cm")),
               linewidth = 1) +
  annotate("text", x = -12, y = 5.5, label = "More Variable Temp./Precip.", size = 3) +
  annotate("text", x = -12, y = -5.5, label = "More Stable Temp./Precip.", size = 3) +
  annotate("text", x = 9, y = -7, label = "Warmer", size = 3) +
  annotate("text", x = -12, y = -7, label = "Wetter", size = 3)+
  theme_minimal() +
  
  # new_scale_fill() +
  # scale_fill_manual(name = "Climate Period",  
  #                   values = c("Current" = "gray", "Future" = "darkred"),
  #                   labels = c("Current", "Future")) + 
  
  theme(legend.position = "top",
        legend.box = "vertical",  # Stack different legend groups vertically
        legend.box.just = "left",
        panel.grid = element_blank()) +  # Removes all grid lines
  guides(
    color = guide_legend(order = 1),   # Ensures "Taxa" appears first
    shape = guide_none()
    # fill = guide_legend(order = 2)
    # linetype = guide_legend(order = 2) # Ensures "Climate Period" appears below
  ) +
  labs(title = "Rainier Target Taxa Plotted Within Current and Projected Climate Envelopes",
       x = "PC1", y = "PC2")

# Wenatchee Climate PCA

W_CurrentClimate_matrix <- as.matrix(Wenatchee_1961_1990_Biovars)
W_FutureClimate_matrix <- as.matrix(Wenatchee_2071_2100_Biovars)

W_ClimateMatrix <- rbind(W_CurrentClimate_matrix, W_FutureClimate_matrix)

# Create a grouping factor: first raster labeled as "Current", second as "Future (2071-2100)"
W_Time_period_labels <- factor(c(rep("Current", nrow(W_CurrentClimate_matrix)), rep("Future (2071-2100)", nrow(W_FutureClimate_matrix))))

#Run PCA (This function only calculates the first 2 PCs)
W_Climate_PCA <- prcomp_irlba(W_ClimateMatrix, n = 2, center = TRUE, scale. = TRUE)

# Extract first two principal components as a dataframe
W_Climate_PCA_df <- data.frame(
  PC1 = W_Climate_PCA$x[,1],
  PC2 = W_Climate_PCA$x[,2],
  Source = W_Time_period_labels  # Label data source
)

# Occurrence Data
W_Occurrence <- filter(Occurrence, scientific_name == "Androsace nivalis" | 
                                   scientific_name == "Chaenactis thompsonii" | 
                                   scientific_name == "Claytonia megarhiza" | 
                                   scientific_name == "Lomatium cuspidatum" | 
                                   scientific_name == "Oreocarya thompsonii")

W_cell_indices <- as.matrix(W_Occurrence[, c("Longitude", "Latitude")]) %>% # Convert to matrix for cell lookup
  cellFromXY(Wenatchee_1961_1990_Biovars, .)

W_Occurrence_pcs <- as.data.frame(W_Climate_PCA$x[W_cell_indices, 1:2]) %>% # Extract PCA values (PC1 and PC2) for these cells and bring in taxon names
  mutate(scientific_name = W_Occurrence$scientific_name)

# Plot
WenatcheeClimatePCAplot <- ggplot() +
  # Hexbin for Current Climate (red gradient)
  geom_hex(data = subset(W_Climate_PCA_df, Source == "Current"),  
           aes(x = PC1, y = PC2, fill = after_stat(count)),  
           bins = 50, alpha = 0.5) +
  scale_fill_gradient(name = "Climate Period", low = "lightgray", high = "black", guide = "none") +
  
  # Reset the fill scale before adding Future Climate
  new_scale_fill() +
  
  # Hexbin for Future Climate (gray gradient)
  geom_hex(data = subset(W_Climate_PCA_df, Source == "Future (2071-2100)"),  
           aes(x = PC1, y = PC2, fill = after_stat(count)),  
           bins = 50, alpha = 0.5) +
  scale_fill_gradient(name = "Climate Period", low = "lightpink", high = "darkred", guide = "none") +
  
  # Manual Colors for Taxa
  scale_color_manual(name = "Taxon",
                     values = c("purple", "cyan", "blue","orange","darkgreen")) +
  
  # Points for Rainier Occurrence (highlighted taxa)
  geom_point(data = W_Occurrence_pcs,  
             aes(x = PC1, y = PC2, color = scientific_name, shape = scientific_name), size = 3, stroke = 1) +
  scale_shape_manual(values = c(8, 1, 0, 2, 3)) + #add 3 and 8 for plus sign and star
  # #Ellipses for Highlighted Taxa
  stat_ellipse(data = W_Occurrence_pcs,
               aes(x = PC1, y = PC2, color = scientific_name),
               level = 0.95, size = 1.2) +
  # stat_chull(data = R_Occurrence_pcs,
  #            geom = "polygon",
  #            aes(x = PC1, y = PC2, color = scientific_name),
  #            alpha = 0.001,
  #            show.legend = FALSE) +
  # 
  new_scale_color() +
  
  #Ellipses for Current and Future Climate Envelopes
  # stat_ellipse(data = R_Climate_PCA_df,
  #              aes(x = PC1, y = PC2, color = Source),
  #              level = 0.95, size = 1.2) +
  # scale_color_manual(name = "Climate Period",
  #                    values = c("Current" = "black", "Future (2071-2100)" = "darkred")) +
  # stat_chull(data = R_Climate_PCA_df,
  #            geom = "polygon",
  #            aes(x = PC1, y = PC2, color = Source),
  #            alpha = 0.2, 
  #            show.legend = FALSE) + 
  # scale_color_manual(name = "Climate Period",
  #                    values = c("Current" = "black", "Future (2071-2100)" = "darkred")) +
  geom_segment(aes(x = -10, xend = 10, y = -5.5, yend = -5.5), 
               arrow = arrow(type = "closed", ends = "both", length = unit(0.3, "cm")),
               linewidth = 1) +
  geom_segment(aes(x = -11, xend = -11, y = -4, yend = 5), 
               arrow = arrow(type = "closed", ends = "both", length = unit(0.3, "cm")),
               linewidth = 1) +
  annotate("text", x = -9, y = 5.5, label = "More Stable Temp./Precip.", size = 3) +
  annotate("text", x = -9, y = -4.5, label = "More Variable Temp./Precip.", size = 3) +
  annotate("text", x = 9, y = -6, label = "Wetter", size = 3) +
  annotate("text", x = -10.5, y = -6, label = "Warmer", size = 3)+
  theme_minimal() +
  
  # new_scale_fill() +
  # scale_fill_manual(name = "Climate Period",  
  #                   values = c("Current" = "gray", "Future" = "darkred"),
  #                   labels = c("Current", "Future")) + 
  
  theme(legend.position = "top",
        legend.box = "horizontal",  # Stack different legend groups vertically
        legend.box.just = "left",
        legend.text = element_text(size = 8),  # Adjust text size to fit everything
        legend.title = element_text(size = 9),  # Slightly larger title for readability
        legend.justification = c(0, 0.5),
        legend.margin = margin(0, -50, 0, 0),
        panel.grid = element_blank()) +  # Removes all grid lines
  guides(
    color = guide_legend(order = 1, nrow = 2, byrow = TRUE),   # Ensures "Taxa" appears first
    shape = guide_none()
    # fill = guide_legend(order = 2)
    # linetype = guide_legend(order = 2) # Ensures "Climate Period" appears below
  ) +
  labs(title = "Wenatchee Target Taxa Plotted Within Current and Projected Climate Envelopes",
       x = "PC1", y = "PC2")

# Biplot for interpretation

W_loadings <- W_Climate_PCA$rotation

# Extract the first two components of the loadings for plotting
W_loadings_df <- as.data.frame(W_loadings[, 1:2])
W_loadings_df$Variable <- rownames(W_loadings_df)

# Scale loadings for better visualization (optional)
W_loadings_df$PC1 <- W_loadings_df$PC1 * 5
W_loadings_df$PC2 <- W_loadings_df$PC2 * 5

# Plot vectors with better label placement
W_Climate_Biplot <- ggplot(W_loadings_df, aes(x = 0, y = 0, xend = PC1, yend = PC2, label = Variable)) +
  geom_segment(arrow = arrow(length = unit(0.2, "cm")), color = "blue") +
  geom_text_repel(aes(x = PC1, y = PC2), size = 4, max.overlaps = 30) +
  labs(title = "Wenatchee Climate Biplot",
       x = "PC1", y = "PC2") +
  theme_minimal()



# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
