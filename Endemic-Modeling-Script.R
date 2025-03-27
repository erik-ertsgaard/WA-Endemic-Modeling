##### WASHINGTON ENDEMIC MODELING PROJECT #####
#### Gjording et al.

# 1.0 LOAD ITEMS -----------------------------------------------------------

# 1.1 Load Packages ----

install.packages("biomod2")
install.packages("tidyverse")
install.packages("rinat")
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
install.packages("xgboost")
install.packages("randomForest")
install.packages("maxnet")
install.packages("mda")
install.packages("earth")
install.packages("devtools")
install.packages("vegan")
install.packages("tidyterra")
install.packages("ggtext")

library(biomod2)
library(tidyverse)
library(rinat)
#library(ClimateNAr)
library(terra)
library(sp)
library(raster)
library(dismo)
library(stats)
library(sf)
library(rnaturalearth)
library(ggspatial)
library(prettymapr)
library(devtools)
library(vegan)
library(tidyterra)
library(ggtext)


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
      dplyr::select(-is_nad83, -is_nad27)
    
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

prepare_biomod_data <- function(presence_data, species_name, pseudo_absence_data, env.var) {
  
  # Convert the SpatVector to an sf object
  occurrence_sf <- st_as_sf(presence_data)
  
  # Filter the data for the species of interest
  occurrence_data <- occurrence_sf %>%
    filter(scientific_name == species_name) %>%
    dplyr::select(geometry)  # Keep only the geometry column
  
  # Extract coordinates (longitude and latitude) from the geometry column
  coords <- st_coordinates(occurrence_data)
  coords <- as.data.frame(coords)  # Convert to data.frame
  colnames(coords) <- c("Longitude", "Latitude")
  
  # Add a "Presence" column with 1 values for all occurrence records
  coords$Presence <- 1
  
  # Extract coordinates from the PA dataset
  coords.pa <- st_coordinates(pseudo_absence_data$geometry)
  pseudo_absence_data$Longitude <- coords.pa[, 1]  
  pseudo_absence_data$Latitude <- coords.pa[, 2]
  
  pseudo_absence_data$Presence <- NA  # Add a Presence column to pseudo-absence data with NA values
  
  # Combine the occurrence data with pseudo-absence data
  combined_data <- bind_rows(coords, pseudo_absence_data)
  combined_data <- as.data.frame(combined_data)
  
  # Format for BIOMOD_FormattingData
  resp.var <- combined_data$Presence  # Response variable (1 for presence, NA for pseudo-absence)
  resp.xy <- combined_data[, c("Longitude", "Latitude")]  # Coordinates for response variable
  PA.user.table <- combined_data[, (colnames(combined_data) %in% c("PA1", "PA2", "PA3", "PA4", "PA5", "PA6", "PA7", "PA8", "PA9", "PA10"))]  
  PA.user.table <- PA.user.table %>%
    mutate_all(~ ifelse(is.na(.), TRUE, .))
  PA.user.table <- as.data.frame(PA.user.table)
  
  # Input into biomod package's formatting function
  formatted_data <- BIOMOD_FormatingData(
    resp.var = resp.var,             # Species presence/absence
    expl.var = env.var,              # Environmental data 
    PA.nb.rep = 10,                   # 10 repitions in pseudo-absence table
    PA.strategy = "user.defined",    # Use your own pseudo-absence table
    PA.user.table = PA.user.table,   # User-defined pseudo-absence table
    resp.xy = resp.xy,              # Coordinates of species presence points
    resp.name = species_name,      # Species name
    filter.raster = TRUE,           # removing occurrences in the same cell
    na.rm = TRUE,
    dir.name = "Modeling/")
  
  # Return the formatted data for further use
  return(formatted_data)
}

biomod_single_models <- function(bm.data) { 

single.models <-  BIOMOD_Modeling(bm.format = bm.data,      # data from biomod_modeling 
                                  models = c("GLM", "RF", "ANN", "MARS"), # using the six most ubiquitous models in SDM
                                  CV.strategy = "kfold",  # using k-fold cross-validation instead of random calibration splits
                                  CV.k = 3,               # with 3 folds
                                  CV.nb.rep = 2,          # and 2 repititions
                                  CV.do.full.models = TRUE, # computing full models
                                  metric.eval = c("TSS", "ROC"),  # evaluating by standard model-performance metrics, TSS and ROC
                                  OPT.strategy = "bigboss",   # tuning parameters set by biomod2 authors, optimized for each model
                                  OPT.data.type = "binary",   # default data type (i.e., no abundance)
                                  var.import = 3)             # 3 permutations to test variable importance
  
return(single.models)

}

biomod_ensemble_models <- function(bm.modeling.output) {
  
ensemble.models <- BIOMOD_EnsembleModeling(bm.mod = bm.modeling.output,    # data from single.modes step above
                                           models.chosen = "all",     # keeping all models
                                           em.by = "all",             # combining all model runs
                                           em.algo = 'EMwmean',       # combining by mean values weighted by model performance
                                           metric.select = "TSS",     # excluding single models on the basis of TSS & ROC
                                           metric.select.thresh = 0.6,  # which both must be above 0.8
                                           metric.eval = "TSS",       # weights for WM combination correspond to TSS values
                                           var.import = 3,            # 3 permutations to test variable importance
                                           EMci.alpha = 0.05,         # alpha value of 0.05 for confidence interval
                                           EMwmean.decay = 'proportional')      # "weights are assigned to each model proportionally to their evaluation scores"

return(ensemble.models)

}

biomod_ensemble_forecast <- function(bm.ensemble.output, proj.name, new.env)  {
  
ensemble.forecast <- BIOMOD_EnsembleForecasting(bm.em = bm.ensemble.output,
                                                proj.name = proj.name,
                                                new.env = new.env,
                                                models.chosen = 'all',
                                                metric.binary = 'all',
                                                metric.filter = 'all',
                                                na.rm = FALSE,
                                                nb.cpu = 1)

return(ensemble.forecast)
  
}

# 1.3 Load Data ----

# 1.3a Load Predictor Data ----

## Set working directory (any desired working directory for writing/reading files)

directory <- "C:/.../" #eg. "C:/Users/username/WA-Endemic-Modeling/"
setwd(directory)

## Create Shapefiles and SpatVectors for Study Extent

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
a
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
W.lith <- rasterize(Lith, WenatcheeDEM, field = "Map_Unit_ID")

levels(R.lith) <- data.frame(ID = seq_along(categories), Category = categories)
levels(W.lith) <- data.frame(ID = seq_along(categories), Category = categories)

## Tree Canopy Cover Data
R.Canopy <- rast(paste0(directory, "USANLCDTreeCanopyCover_Ra.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=5) %>% 
  resample(RainierDEM, method = "bilinear")

W.Canopy <- rast(paste0(directory, "USANLCDTreeCanopyCover_We.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=5) %>% 
  resample(WenatcheeDEM, method="bilinear")

## Topological Data
R.TPI.20 <- rast(paste0(directory, "RainierDEM10_TPI_20m.tif"))
R.TPI.30 <- rast(paste0(directory, "RainierDEM10_TPI_30m.tif"))
R.TPI.50 <- rast(paste0(directory, "RainierDEM10_TPI_50m.tif"))
R.TPI.100 <- rast(paste0(directory, "RainierDEM10_TPI_100m.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
R.TPI.300 <- rast(paste0(directory, "RainierDEM10_TPI_300m.tif"))
R.TPI.500 <- rast(paste0(directory, "RainierDEM10_TPI_500m.tif"))


W.TPI.100 <- rast(paste0(directory, "WenatcheesDEM10_TPI100m.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")


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
R.Eastness <- rast(paste0(directory, "Eastness_Rainier.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
R.NorthEastness <- rast(paste0(directory, "NorthEastness_Rainier.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(RainierDEM, method="bilinear")
W.Northness <- rast(paste0(directory, "Northness_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")
W.Eastness <- rast(paste0(directory, "Eastness_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")
W.NorthEastness <- rast(paste0(directory, "NorthEastness_Wenatchees.tif")) %>% 
  project("+proj=longlat +datum=WGS84") %>% 
  disagg(fact=3) %>% 
  resample(WenatcheeDEM, method="bilinear")

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
                  get_inat_obs(taxon_name = "Claytonia megarhiza nivalis")) %>% 
  filter(quality_grade == "research") %>%
  filter(positional_accuracy < 30) %>%
  filter(!is.na(positional_accuracy)) %>%
  filter(coordinates_obscured == "false") #applying filters to ensure identification and spatial accuracy

# creating a file for reference in the repository, may be replaced when more response data is appended
write_csv(inat.all, file = "Data/inat-response-data.csv")

inat.all <- read_csv("Data/inat-response-data.csv")

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
occurrence_data_clipped$scientific_name[occurrence_data_clipped$scientific_name == "Claytonia megarhiza"] <- "Claytonia megarhiza nivalis"

# removing duplicative occurrences within 1 meter of each other
occurrence_data_cleaned <- occurrence_data_clipped %>%
  group_by(scientific_name) %>%
  group_modify(~ remove_nearby_points(.x, dist_threshold = 5)) %>%
  ungroup() %>% 
  filter(uncertainty <= 1000) %>%
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
###PCA
RainierClimateMatrix <- as.matrix(values(Rainier_1961_1990_Biovars))

RainierPCA <- prcomp(na.omit(RainierClimateMatrix), scale. = TRUE)

summary(RainierPCA)

Rainierpcs_1_3_scores <- predict(RainierBiovarsAll, RainierPCA, index = 1:3)
###


# 2.1 Adjusting Predictor Variables ----

# loading Rainier variables from prior work, combining continuous and discrete raster layers
Rainier_PCs_current <- rast("Data/Large-Files/Rainier_3s_expl_PCs1-6.tif")

Rainier_geology <- read_rds("Data/Large-Files/Rainier-categorical-vars.rds") %>%
  unwrap()

# removing spaces in column names
names(Rainier_geology) <- gsub(" ", "", names(Rainier_geology))

# Removing factor levels in categorical layers that are rare on the landscape(in our response data) and would impact modeling efficacy
rainier.background.data1 <- vect(rainier.background.data) # converting background data to a SpatVector

levelsdf.r <- extract(Rainier_geology, rainier.background.data1) # extracting variable identities for every point

mask_ice <- Rainier_geology$Lithology == "ice" # creating this mask for later

# Finding which levels occur for less than 5% of the background data
soiltype_counts.r <- levelsdf.r %>%
  group_by(SoilType) %>%
  summarize(count = n(), proportion = n()/nrow(rainier.background.data1))

# Conerting the name of all rare levels to "other
rare_soils.r <- soiltype_counts.r$SoilType[soiltype_counts.r$proportion < 0.05]

soil_levels.r <- levels(Rainier_geology$SoilType)[[1]]

soil_levels.r$comp_tx_86[soil_levels.r$comp_tx_86 %in% rare_soils.r] <- "other"

levels(Rainier_geology$SoilType) <- soil_levels.r  

Rainier_geology$SoilType <- ifel(
  Rainier_geology$SoilType %in% c(4, 6),  # These are  duplicate "other" levels
  2,                                           # Reassign them all to ID 2 ("other")
  Rainier_geology$SoilType                      # Keep everything else as is
) %>%
  droplevels()


# Repeating for Lithology layer
lithology_counts <- levelsdf.r %>%
  group_by(Lithology) %>%
  summarize(count = n(), proportion = n()/nrow(rainier.background.data1))

rare_lithologies <- lithology_counts$Lithology[lithology_counts$proportion < 0.05]

lith_levels <- levels(Rainier_geology$Lithology)[[1]]

lith_levels$Category[lith_levels$Category %in% rare_lithologies] <- "other"

levels(Rainier_geology$Lithology) <- lith_levels  

Rainier_geology$Lithology <- ifel(
  Rainier_geology$Lithology %in% c(4, 5, 7, 8, 10, 11, 12, 13),  # These are  duplicate "other" levels
  1,                                           # Reassign them all to ID 1 ("other")
  Rainier_geology$Lithology                      # Keep everything else as is
) %>%
  droplevels()

# combining explanatory data into one raster stack
Rainier_expl_current <- c(Rainier_geology, Rainier_PCs_current)

# repeating steps for ssp245 (moderate emissions scenario)
Rainier_PCs_ssp245 <- rast("Data/Large-Files/Rainier_3s_ssp245_expl_PCs1-6.tif") 

Rainier_expl_ssp245 <- c(Rainier_geology, Rainier_PCs_ssp245)

# repeating steps for ssp585 (business as usual emissions scenario)
Rainier_PCs_ssp585 <- rast("Data/Large-Files/Rainier_3s_ssp585_expl_PCs1-6.tif") 

Rainier_expl_ssp585 <- c(Rainier_geology, Rainier_PCs_ssp585)

# Removing data from all layers of areas covered by ice, with mask layer from earlier
Rainier_expl_current <- mask(x = Rainier_expl_current,
                             mask = mask_ice,
                             maskvalues = TRUE)

Rainier_expl_ssp245 <- mask(x = Rainier_expl_ssp245,
                             mask = mask_ice,
                             maskvalues = TRUE)

Rainier_expl_ssp585 <- mask(x = Rainier_expl_ssp585,
                            mask = mask_ice,
                            maskvalues = TRUE)

# loading Wenatchee variables from prior work, combining continuous and discrete raster layers, and removing spaces in column names
Wenatchee_PCs <- rast("Data/Large-Files/Wenatchee-3s-current-expl-PCs1-6.tif")

Wenatchee_soil <- read_rds("Data/Large-Files/Wenatchee-categorical-vars.rds") %>%
  unwrap()

names(Wenatchee_soil) <- gsub(" ", "", names(Wenatchee_soil))

# Removing factor levels in categorical layers that are rare on the landscape(in our response data) and would impact modeling efficacy
wenatchees.background.data1 <- vect(wenatchees.background.data) # converting background data to a SpatVector

levelsdf <- extract(Wenatchee_soil, wenatchees.background.data1) # extracting variable identities for every point

# Finding which levels occur for less than 5% of the background data
soiltype_counts <- levelsdf %>%
  group_by(SoilType) %>%
  summarize(count = n(), proportion = n()/nrow(wenatchees.background.data1))

# Conerting the name of all rare levels to "other
rare_soils <- soiltype_counts$SoilType[soiltype_counts$proportion < 0.05]

soil_levels <- levels(Wenatchee_soil$SoilType)[[1]]

soil_levels$comp_tx_86[soil_levels$comp_tx_86 %in% rare_soils] <- "other"

levels(Wenatchee_soil$SoilType) <- soil_levels  

Wenatchee_soil$SoilType <- ifel(
  Wenatchee_soil$SoilType %in% c(5, 6, 9, 10),  # These are  duplicate "other" levels
  0,                                           # Reassign them all to ID 0 ("other")
  Wenatchee_soil$SoilType                      # Keep everything else as is
) %>%
  droplevels()



# Repeating for Category (Lithology) layer
category_counts <- levelsdf %>%
  group_by(Category) %>%
  summarize(count = n(), proportion = n()/nrow(wenatchees.background.data1))

rare_categories <- category_counts$Category[category_counts$proportion < 0.05]

cat_levels <- levels(Wenatchee_soil$Category)[[1]]

cat_levels$Category[cat_levels$Category %in% rare_categories] <- "other"

levels(Wenatchee_soil$Category) <- cat_levels  


Wenatchee_soil$Category <- ifel(
  Wenatchee_soil$Category %in% c(4, 10, 12, 14, 15, 18, 19),  
  3,                                           
  Wenatchee_soil$Category                      
) %>%
  droplevels()

# Combining PCs and categorical variables into one raster stack
Wenatchee_expl_current <- c(Wenatchee_soil, Wenatchee_PCs)


# Repeating steps for ssp245
Wenatchee_PCs_ssp245 <- rast("Data/Large-Files/Wenatchee-3s-ssp245-expl-PCs1-6.tif")

Wenatchee_expl_ssp245 <- c(Wenatchee_soil, Wenatchee_PCs_ssp245)


# Repeating steps for ssp585
Wenatchee_PCs_ssp585 <- rast("Data/Large-Files/Wenatchee-3s-ssp585-expl-PCs1-6.tif")

Wenatchee_expl_ssp585 <- c(Wenatchee_soil, Wenatchee_PCs_ssp585)

# 2.2 Adjusting Response Variables ----
## creating a data frame of occurrence and background data compatible with 'biomod2'

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



# 2.3 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# Preparing data for biomod2 formatting requirements
bm.pera <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Pedicularis rainierensis",
                               pseudo_absence_data = rainier.background.data,
                               env.var = Rainier_expl_current)

bm.tast <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Tauschia stricklandii",
                               pseudo_absence_data = rainier.background.data,
                               env.var = Rainier_expl_current)

bm.cacr <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Castilleja cryptantha",
                               pseudo_absence_data = rainier.background.data,
                               env.var = Rainier_expl_current)

bm.clme <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Claytonia megarhiza nivalis",
                               pseudo_absence_data = wenatchees.background.data,
                               env.var = Wenatchee_expl_current)

bm.chth <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Chaenactis thompsonii",
                               pseudo_absence_data = wenatchees.background.data,
                               env.var = Wenatchee_expl_current)

bm.anni <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Androsace nivalis",
                               pseudo_absence_data = wenatchees.background.data,
                               env.var = Wenatchee_expl_current)

bm.locu <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Lomatium cuspidatum",
                               pseudo_absence_data = wenatchees.background.data,
                               env.var = Wenatchee_expl_current)

bm.orth <- prepare_biomod_data(presence_data = occurrence_data_cleaned,
                               species_name = "Oreocarya thompsonii",
                               pseudo_absence_data = wenatchees.background.data,
                               env.var = Wenatchee_expl_current)

# Running biomod2 single models 
pera.sms <- biomod_single_models(bm.data = bm.pera)

#tast.sms <- biomod_single_models(bm.data = bm.tast)

cacr.sms <- biomod_single_models(bm.data = bm.cacr)

clme.sms <- biomod_single_models(bm.data = bm.clme)

anni.sms <- biomod_single_models(bm.data = bm.anni)

chth.sms <- biomod_single_models(bm.data = bm.chth)

locu.sms <- biomod_single_models(bm.data = bm.locu)

orth.sms <- biomod_single_models(bm.data = bm.orth)


# Running biomod2 ensemble models 
pera.em <- biomod_ensemble_models(bm.modeling.output = pera.sms)

#tast.em <- biomod_ensemble_models(bm.modeling.output = tast.sms)

cacr.em <- biomod_ensemble_models(bm.modeling.output = cacr.sms)

clme.em <- biomod_ensemble_models(bm.modeling.output = clme.sms)

anni.em <- biomod_ensemble_models(bm.modeling.output = anni.sms)

chth.em <- biomod_ensemble_models(bm.modeling.output = chth.sms)

locu.em <- biomod_ensemble_models(bm.modeling.output = locu.sms)

orth.em <- biomod_ensemble_models(bm.modeling.output = orth.sms)


# Projecting current suitable habitat
pera.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = pera.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Rainier_expl_current)

tast.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = tast.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Rainier_expl_current)

cacr.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = cacr.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Rainier_expl_current)

clme.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = clme.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_current)

anni.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = anni.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_current)

chth.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = chth.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_current)

locu.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = locu.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_current)

orth.emproj.current <- biomod_ensemble_forecast(bm.ensemble.output = orth.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_current)


# Projecting future suitable habitat in different climate change scenarios
pera.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = pera.em,
                                               proj.name = "ssp245EM", 
                                               new.env = Rainier_expl_ssp245)

tast.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = tast.em,
                                               proj.name = "ssp245EM", 
                                               new.env = Rainier_expl_ssp245)

cacr.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = cacr.em,
                                               proj.name = "ssp245EM", 
                                               new.env = Rainier_expl_ssp245)

clme.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = clme.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp245)

anni.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = anni.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp245)

chth.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = chth.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp245)

locu.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = locu.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp245)

orth.emproj.ssp245 <- biomod_ensemble_forecast(bm.ensemble.output = orth.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp245)


# Repeating for ssp585
pera.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = pera.em,
                                               proj.name = "ssp585EM", 
                                               new.env = Rainier_expl_ssp585)

tast.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = tast.em,
                                               proj.name = "ssp585EM", 
                                               new.env = Rainier_expl_ssp585)

cacr.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = cacr.em,
                                               proj.name = "ssp585EM", 
                                               new.env = Rainier_expl_ssp585)

clme.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = clme.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp585)

anni.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = anni.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp585)

chth.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = chth.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp585)

locu.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = locu.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp585)

orth.emproj.ssp585 <- biomod_ensemble_forecast(bm.ensemble.output = orth.em,
                                                proj.name = "CurrentEM", 
                                                new.env = Wenatchee_expl_ssp585)


######## Move to figures script ##########
# 4.0 Figures --------------------------------------------------------------

# Variable Importance boxplots (Figure S2)
p.pe <- bm_PlotVarImpBoxplot(bm.out = pera.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Pedicularis rainierensis")
p.ca <- bm_PlotVarImpBoxplot(bm.out = cacr.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Castilleja cryptantha")
p.ta <- bm_PlotVarImpBoxplot(bm.out = tast.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Tauschia stricklandii")
p.cl <- bm_PlotVarImpBoxplot(bm.out = clme.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Claytonia megarhiza var. nivalis")
p.an <- bm_PlotVarImpBoxplot(bm.out = anni.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Androsace nivalis")
p.ch <- bm_PlotVarImpBoxplot(bm.out = chth.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Chaenactis thompsonii")
p.lo <- bm_PlotVarImpBoxplot(bm.out = locu.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Lomatium cuspidatum")
p.or <- bm_PlotVarImpBoxplot(bm.out = orth.em, group.by = c('expl.var', 'algo', 'merged.by.run'), main = "Oreocarya thompsonii")

install.packages("cowplot")
library(cowplot)
plot_grid(plotlist = c(p.ca, p.ta, p.cl, p.an, p.ch, p.lo, p.or, p.ca), 
          ncol = 2, nrow = 4)



install.packages("basemaps")
library(basemaps)
get_maptypes()
basemap(ext = RainierExtentPolygon,
        map_service = "esri",
        map_type = "world_hillshade_dark")