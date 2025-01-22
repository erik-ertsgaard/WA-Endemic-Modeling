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
install.packages("sf")

library(biomod2)
library(tidyverse)
library(rinat)
#library(rgbif)
#library(rvest)
library(ClimateNAr)
library(terra)
library(sp)
library(raster)
library(sf)


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


# 1.3 Load Data ----

# 1.3a Load Predictor Data ----

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

Wenatchee <- mosaic(WenatcheeW, WenatcheeE) %>% 
  crop(SpatVectorWenatchee) %>% 
  project("+proj=longlat +datum=WGS84")

Rainier <- crop(Rainier, SpatVectorRainier) %>% 
  project("+proj=longlat +datum=WGS84")

writeRaster(Wenatchee, filename = paste0(directory, "WenatcheeDEM.tif"))

writeRaster(Rainier, filename = paste0(directory, "RainierDEM.tif"))

## Climate Data

climateNAr(inputFile = paste0(directory, "WenatcheeDEM.tif"),
           periodList = "Normal_1961_1990.nrm",
           varList = "YSM",
           outDir = directory)

# 1.3b Loading Response Data ----

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

inat.all <- rename(inat.all, uncertainty = positional_accuracy)

# combining all records into one large dataframe
occurrence_data_raw <- bind_rows(inat.all, cpnwh.wgs84) 

# 2.0 Data Adjustments -----------------------------------------------------

# 2.1 Filtering Occurrence Data ----

library(terra)
# 2.2 Adjusting Predictor Variables ----

# 2.3 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
