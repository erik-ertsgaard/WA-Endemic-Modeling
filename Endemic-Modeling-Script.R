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
install.packages("rnaturalearth")
install.packages("ggspatial")
install.packages("prettymapr")


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
library(rnaturalearth)
library(ggspatial)
library(prettymapr)

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

bm_PseudoAbsences()

# 2.0 Data Adjustments -----------------------------------------------------

# 2.1 Adjusting Predictor Variables ----

# 2.2 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
