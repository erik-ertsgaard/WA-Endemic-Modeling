install.packages("tidyterra")
install.packages("basemaps")
install.packages("tidyverse")
install.packages("biomod2")
install.packages("terra")
install.packages("cowplot")

library(tidyterra)
library(basemaps)
library(tidyverse)
library(biomod2)
library(terra)
library(cowplot)

wen.bm <- basemap_terra(ext = WenatcheeExtentPolygon,
                        map_service = "esri",
                        map_type = "world_hillshade_dark",
                        map_res = 10) %>%
  project("epsg:4326")
  
Wen.ref.points <- data.frame(
  name = c("Leavenworth", "Mt. Stuart"),
  lon = c(-120.661, -120.902372),
  lat = c(47.596, 47.475139))

Wen.ref.sf <- st_as_sf(Wen.ref.points, coords = c("lon", "lat"), crs = 4326)

# Oreocarya thompsonii -----

load("Modeling/Oreocarya.thompsonii/Oreocarya.thompsonii.CurrentEM.ensemble.projection.out")

orth.emproj.current <- Oreocarya.thompsonii.CurrentEM.ensemble.projection.out
plot(orth.emproj.current)

orth.binproj <- get_predictions(orth.emproj.current, metric.binary = "TSS")
plot(orth.binproj)

orth.binproj.f <- terra::ifel(orth.binproj == 0, NA, 1) 

orth.binproj.f <- as.data.frame(orth.binproj.f, xy = TRUE) 

names(orth.binproj.f)[3] <- "Current-Suitable-Habitat"

orth.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM_Oreocarya.thompsonii_ensemble_TSSbin.tif")

orth.binproj.ssp245.f <- terra::ifel(orth.binproj.ssp245 == 0, NA, 1) 

orth.binproj.ssp245.f <- as.data.frame(orth.binproj.ssp245.f, xy = TRUE) 

names(orth.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"



orth.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM_Oreocarya.thompsonii_ensemble_TSSbin.tif")

orth.binproj.ssp585.f <- terra::ifel(orth.binproj.ssp585 == 0, NA, 1) 

orth.binproj.ssp585.f <- as.data.frame(orth.binproj.ssp585.f, xy = TRUE) 

names(orth.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

orth.points <- filter(occurrence_data_cleaned, scientific_name == "Oreocarya thompsonii")

orth.binproj.f$category <- "Current"
orth.binproj.ssp245.f$category <- "SSP245"
orth.binproj.ssp585.f$category <- "SSP585"
orth.points$category <- "Recorded Species Locations"
Wen.ref.sf$category <- "Reference Locations"


orth.plot <- ggplot() + 
  geom_spatraster_rgb(data = wen.bm) +
  geom_tile(data = orth.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = orth.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.5) +
  geom_tile(data = orth.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.5) +
  geom_sf(data = orth.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Wen.ref.sf, aes(color = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink", 
                               "SSP585" = "magenta")) +
  scale_color_manual(values = c("Recorded Species Locations" = "forestgreen", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Wen.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Oreocarya thompsonii",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.25, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

orth.rs245 <- BIOMOD_RangeSize(proj.current = orth.binproj,
                               proj.future = orth.binproj.ssp245)

orth.rs585 <- BIOMOD_RangeSize(proj.current = orth.binproj,
                               proj.future = orth.binproj.ssp585)



# Androsace nivalis ----  
  
anni.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Androsace.nivalis_ensemble_TSSbin.tif")

anni.binproj.f <- terra::ifel(anni.binproj == 0, NA, 1) 

anni.binproj.f <- as.data.frame(anni.binproj.f, xy = TRUE) 

names(anni.binproj.f)[3] <- "Current-Suitable-Habitat"


anni.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Androsace.nivalis_ensemble_TSSbin.tif")

anni.binproj.ssp245.f <- terra::ifel(anni.binproj.ssp245 == 0, NA, 1) 

anni.binproj.ssp245.f <- as.data.frame(anni.binproj.ssp245.f, xy = TRUE) 

names(anni.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


anni.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Androsace.nivalis_ensemble_TSSbin.tif")

anni.binproj.ssp585.f <- terra::ifel(anni.binproj.ssp585 == 0, NA, 1) 

anni.binproj.ssp585.f <- as.data.frame(anni.binproj.ssp585.f, xy = TRUE) 

names(anni.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

anni.points <- filter(occurrence_data_cleaned, scientific_name == "Androsace nivalis")

anni.binproj.f$category <- "Current"
anni.binproj.ssp245.f$category <- "SSP245"
anni.binproj.ssp585.f$category <- "SSP585"
anni.points$category <- "Recorded Species Locations"

anni.plot <- ggplot() + 
  geom_spatraster_rgb(data = wen.bm) +
  geom_tile(data = anni.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = anni.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.5) +
  geom_tile(data = anni.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.15) +
  geom_sf(data = anni.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Wen.ref.sf, aes(color = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink", 
                               "SSP585" = "magenta")) +
  scale_color_manual(values = c("Recorded Species Locations" = "purple", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Wen.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Androsace nivalis",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.25, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

anni.rs245 <- BIOMOD_RangeSize(proj.current = anni.binproj,
                               proj.future = anni.binproj.ssp245)

anni.rs585 <- BIOMOD_RangeSize(proj.current = anni.binproj,
                               proj.future = anni.binproj.ssp585)


# Chaenactis thompsonii ----  

chth.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Chaenactis.thompsonii_ensemble_TSSbin.tif")

chth.binproj.f <- terra::ifel(chth.binproj == 0, NA, 1) 

chth.binproj.f <- as.data.frame(chth.binproj.f, xy = TRUE) 

names(chth.binproj.f)[3] <- "Current-Suitable-Habitat"


chth.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Chaenactis.thompsonii_ensemble_TSSbin.tif")

chth.binproj.ssp245.f <- terra::ifel(chth.binproj.ssp245 == 0, NA, 1) 

chth.binproj.ssp245.f <- as.data.frame(chth.binproj.ssp245.f, xy = TRUE) 

names(chth.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


chth.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Chaenactis.thompsonii_ensemble_TSSbin.tif")

chth.binproj.ssp585.f <- terra::ifel(chth.binproj.ssp585 == 0, NA, 1) 

chth.binproj.ssp585.f <- as.data.frame(chth.binproj.ssp585.f, xy = TRUE) 

names(chth.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

chth.points <- filter(occurrence_data_cleaned, scientific_name == "Chaenactis thompsonii")

chth.binproj.f$category <- "Current"
chth.binproj.ssp245.f$category <- "SSP245"
chth.binproj.ssp585.f$category <- "SSP585"
chth.points$category <- "Recorded Species Locations"

chth.plot <- ggplot() + 
  geom_spatraster_rgb(data = wen.bm) +
  geom_tile(data = chth.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = chth.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.5) +
  geom_tile(data = chth.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.3) +
  geom_sf(data = chth.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Wen.ref.sf, aes(color = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink", 
                               "SSP585" = "magenta")) +
  scale_color_manual(values = c("Recorded Species Locations" = "steelblue1", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Wen.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Chaenactis thompsonii",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.25, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

chth.rs245 <- BIOMOD_RangeSize(proj.current = chth.binproj,
                               proj.future = chth.binproj.ssp245)

chth.rs585 <- BIOMOD_RangeSize(proj.current = chth.binproj,
                               proj.future = chth.binproj.ssp585)

# Lomatium cuspidatum ----  

locu.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Lomatium.cuspidatum_ensemble_TSSbin.tif")

locu.binproj.f <- terra::ifel(locu.binproj == 0, NA, 1) 

locu.binproj.f <- as.data.frame(locu.binproj.f, xy = TRUE) 

names(locu.binproj.f)[3] <- "Current-Suitable-Habitat"


locu.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Lomatium.cuspidatum_ensemble_TSSbin.tif")

locu.binproj.ssp245.f <- terra::ifel(locu.binproj.ssp245 == 0, NA, 1) 

locu.binproj.ssp245.f <- as.data.frame(locu.binproj.ssp245.f, xy = TRUE) 

names(locu.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


locu.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Lomatium.cuspidatum_ensemble_TSSbin.tif")

locu.binproj.ssp585.f <- terra::ifel(locu.binproj.ssp585 == 0, NA, 1) 

locu.binproj.ssp585.f <- as.data.frame(locu.binproj.ssp585.f, xy = TRUE) 

names(locu.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

locu.points <- filter(occurrence_data_cleaned, scientific_name == "Lomatium cuspidatum")

locu.binproj.f$category <- "Current"
locu.binproj.ssp245.f$category <- "SSP245"
locu.binproj.ssp585.f$category <- "SSP585"
locu.points$category <- "Recorded Species Locations"

locu.plot <- ggplot() + 
  geom_spatraster_rgb(data = wen.bm) +
  geom_tile(data = locu.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = locu.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.5) +
  geom_tile(data = locu.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.3) +
  geom_sf(data = locu.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Wen.ref.sf, aes(color = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink", 
                               "SSP585" = "magenta")) +
  scale_color_manual(values = c("Recorded Species Locations" = "orange", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Wen.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Lomatium cuspidatum",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.25, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

locu.rs245 <- BIOMOD_RangeSize(proj.current = locu.binproj,
                               proj.future = locu.binproj.ssp245)
locu.rs245

locu.rs585 <- BIOMOD_RangeSize(proj.current = locu.binproj,
                               proj.future = locu.binproj.ssp585)
locu.rs585




# Claytonia megarhiza var. nivalis ----  

clme.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Claytonia.megarhiza.nivalis_ensemble_TSSbin.tif")

clme.binproj.f <- terra::ifel(clme.binproj == 0, NA, 1) 

clme.binproj.f <- as.data.frame(clme.binproj.f, xy = TRUE) 

names(clme.binproj.f)[3] <- "Current-Suitable-Habitat"


clme.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Claytonia.megarhiza.nivalis_ensemble_TSSbin.tif")

clme.binproj.ssp245.f <- terra::ifel(clme.binproj.ssp245 == 0, NA, 1) 

clme.binproj.ssp245.f <- as.data.frame(clme.binproj.ssp245.f, xy = TRUE) 

names(clme.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


clme.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Claytonia.megarhiza.nivalis_ensemble_TSSbin.tif")

clme.binproj.ssp585.f <- terra::ifel(clme.binproj.ssp585 == 0, NA, 1) 

clme.binproj.ssp585.f <- as.data.frame(clme.binproj.ssp585.f, xy = TRUE) 

names(clme.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

clme.points <- filter(occurrence_data_cleaned, scientific_name == "Claytonia megarhiza nivalis")

clme.binproj.f$category <- "Current"
clme.binproj.ssp245.f$category <- "SSP245"
clme.binproj.ssp585.f$category <- "SSP585"
clme.points$category <- "Recorded Species Locations"

clme.plot <- ggplot() + 
  geom_spatraster_rgb(data = wen.bm) +
  geom_tile(data = clme.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = clme.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.6) +
  geom_tile(data = clme.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.2) +
  geom_sf(data = clme.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Wen.ref.sf, aes(color = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink", 
                               "SSP585" = "magenta")) +
  scale_color_manual(values = c("Recorded Species Locations" = "darkblue", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Wen.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Claytonia megarhiza var. nivalis",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.25, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

clme.rs245 <- BIOMOD_RangeSize(proj.current = clme.binproj,
                               proj.future = clme.binproj.ssp245)
clme.rs245

clme.rs585 <- BIOMOD_RangeSize(proj.current = clme.binproj,
                               proj.future = clme.binproj.ssp585)
clme.rs585

### Rainier species -----

rain.bm <- basemap_terra(ext = RainierExtentPolygon,
                        map_service = "esri",
                        map_type = "world_hillshade_dark",
                        map_res = 10) %>%
  project("epsg:4326")

Rain.ref.points <- data.frame(
  name = c("Paradise", "Sunrise"),
  lon = c(-121.7424528, -121.643463),
  lat = c(46.78617778, 46.914454))

Rain.ref.sf <- st_as_sf(Rain.ref.points, coords = c("lon", "lat"), crs = 4326)


mask.ice.df <- as.data.frame(mask_ice, xy = TRUE)


# Pedicularis rainierensis ----  

pera.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Pedicularis.rainierensis_ensemble_TSSbin.tif")

pera.binproj.f <- terra::ifel(pera.binproj == 0, NA, 1) 

pera.binproj.f <- as.data.frame(pera.binproj.f, xy = TRUE) 

names(pera.binproj.f)[3] <- "Current-Suitable-Habitat"


pera.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Pedicularis.rainierensis_ensemble_TSSbin.tif")

pera.binproj.ssp245.f <- terra::ifel(pera.binproj.ssp245 == 0, NA, 1) 

pera.binproj.ssp245.f <- as.data.frame(pera.binproj.ssp245.f, xy = TRUE) 

names(pera.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


pera.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Pedicularis.rainierensis_ensemble_TSSbin.tif")

pera.binproj.ssp585.f <- terra::ifel(pera.binproj.ssp585 == 0, NA, 1) 

pera.binproj.ssp585.f <- as.data.frame(pera.binproj.ssp585.f, xy = TRUE) 

names(pera.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

pera.points <- filter(occurrence_data_cleaned, scientific_name == "Pedicularis rainierensis")

pera.binproj.f$category <- "Current"
pera.binproj.ssp245.f$category <- "SSP245"
pera.binproj.ssp585.f$category <- "SSP585"
pera.points$category <- "Recorded Species Locations"
Rain.ref.sf$category <- "Reference Locations"


pera.plot <- ggplot() + 
  geom_spatraster_rgb(data = rain.bm) +
  geom_tile(data = pera.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = pera.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.5) +
  geom_tile(data = pera.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.3) +
  #geom_tile(data = mask_ice, aes(x = x, y = y), fill = "white", alpha = 0.9) +
  geom_sf(data = pera.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Rain.ref.sf, aes(colour = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink", 
                               "SSP585" = "magenta")) +
  scale_color_manual(values = c("Recorded Species Locations" = "forestgreen", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Rain.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Pedicularis rainierensis",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.3, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

pera.rs245 <- BIOMOD_RangeSize(proj.current = pera.binproj,
                               proj.future = pera.binproj.ssp245)
pera.rs245

pera.rs585 <- BIOMOD_RangeSize(proj.current = pera.binproj,
                               proj.future = pera.binproj.ssp585)
pera.rs585


# Tauschia stricklandii ----  

tast.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Tauschia.stricklandii_ensemble_TSSbin.tif")

tast.binproj.f <- terra::ifel(tast.binproj == 0, NA, 1) 

tast.binproj.f <- as.data.frame(tast.binproj.f, xy = TRUE) 

names(tast.binproj.f)[3] <- "Current-Suitable-Habitat"


tast.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Tauschia.stricklandii_ensemble_TSSbin.tif")

tast.binproj.ssp245.f <- terra::ifel(tast.binproj.ssp245 == 0, NA, 1) 

tast.binproj.ssp245.f <- as.data.frame(tast.binproj.ssp245.f, xy = TRUE) 

names(tast.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


tast.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Tauschia.stricklandii_ensemble_TSSbin.tif")

tast.binproj.ssp585.f <- terra::ifel(tast.binproj.ssp585 == 0, NA, 1) 

tast.binproj.ssp585.f <- as.data.frame(tast.binproj.ssp585.f, xy = TRUE) 

names(tast.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

tast.points <- filter(occurrence_data_cleaned, scientific_name == "Tauschia stricklandii")

tast.binproj.f$category <- "Current"
tast.binproj.ssp245.f$category <- "SSP245"
#tast.binproj.ssp585.f$category <- "SSP585"
tast.points$category <- "Recorded Species Locations"


tast.plot <- ggplot() + 
  geom_spatraster_rgb(data = rain.bm) +
  geom_tile(data = tast.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  geom_tile(data = tast.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  #geom_tile(data = tast.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.3) +
  #geom_tile(data = mask_ice, aes(x = x, y = y), fill = "white", alpha = 0.9) +
  geom_sf(data = tast.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Rain.ref.sf, aes(colour = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey", 
                               "SSP245" = "pink")) +
  scale_color_manual(values = c("Recorded Species Locations" = "orange", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Rain.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Tauschia stricklandii",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.3, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

tast.rs245 <- BIOMOD_RangeSize(proj.current = tast.binproj,
                               proj.future = tast.binproj.ssp245)
tast.rs245

tast.rs585 <- BIOMOD_RangeSize(proj.current = tast.binproj,
                               proj.future = tast.binproj.ssp585)
tast.rs585

# Castilleja cryptantha ----  

cacr.binproj <- rast("Modeling/Final-Projections/proj_currentEM.final_Castilleja.cryptantha_ensemble_TSSbin.tif")

cacr.binproj.f <- terra::ifel(cacr.binproj == 0, NA, 1) 

cacr.binproj.f <- as.data.frame(cacr.binproj.f, xy = TRUE) 

names(cacr.binproj.f)[3] <- "Current-Suitable-Habitat"


cacr.binproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Castilleja.cryptantha_ensemble_TSSbin.tif")

cacr.binproj.ssp245.f <- terra::ifel(cacr.binproj.ssp245 == 0, NA, 1) 

cacr.binproj.ssp245.f <- as.data.frame(cacr.binproj.ssp245.f, xy = TRUE) 

names(cacr.binproj.ssp245.f)[3] <- "ssp245-Suitable-Habitat"


cacr.binproj.ssp585 <- rast("Modeling/Final-Projections/proj_ssp585EM.final_Castilleja.cryptantha_ensemble_TSSbin.tif")

cacr.binproj.ssp585.f <- terra::ifel(cacr.binproj.ssp585 == 0, NA, 1) 

cacr.binproj.ssp585.f <- as.data.frame(cacr.binproj.ssp585.f, xy = TRUE) 

names(cacr.binproj.ssp585.f)[3] <- "ssp585-Suitable-Habitat"

cacr.points <- filter(occurrence_data_cleaned, scientific_name == "Castilleja cryptantha")

cacr.binproj.f$category <- "Current"
#cacr.binproj.ssp245.f$category <- "SSP245"
#cacr.binproj.ssp585.f$category <- "SSP585"
cacr.points$category <- "Recorded Species Locations"


cacr.plot <- ggplot() + 
  geom_spatraster_rgb(data = rain.bm) +
  geom_tile(data = cacr.binproj.f, aes(x = x, y = y, fill = category), alpha = 0.8) +
  #geom_tile(data = cacr.binproj.ssp245.f, aes(x = x, y = y, fill = category), alpha = 0.65) +
  #geom_tile(data = cacr.binproj.ssp585.f, aes(x = x, y = y, fill = category), alpha = 0.3) +
  #geom_tile(data = mask_ice, aes(x = x, y = y), fill = "white", alpha = 0.9) +
  geom_sf(data = cacr.points, aes(colour = category), size = 3, shape = 10) +
  geom_sf(data = Rain.ref.sf, aes(colour = category), size = 3, shape = 20) + 
  scale_fill_manual(values = c("Current" = "lightgrey")) +
  scale_color_manual(values = c("Recorded Species Locations" = "darkblue", 
                                "Reference Locations" = "black")) +
  geom_sf_label(data = Rain.ref.sf, aes(geometry = geometry, label = name), nudge_y = 0.007, nudge_x = 0.005) +
  labs(title = "Castilleja cryptantha",
       x = "Longitude", y = "Latitude", fill = "Modeled Suitable Habitat", color = "Points") + 
  annotation_scale(location = "bl", width_hint = 0.3, 
                   text_col = "black", line_col = "black", 
                   height = unit(0.25, "cm"), pad_x = unit(0.6, "cm"), 
                   pad_y = unit(0.1, "cm"), bar_cols = c("white", "gray")) +  
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering, 
                         height = unit(0.7, "cm"), width = unit(0.4, "cm")) +  
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold.italic"))

cacr.rs245 <- BIOMOD_RangeSize(proj.current = cacr.binproj,
                               proj.future = cacr.binproj.ssp245)
cacr.rs245

cacr.rs585 <- BIOMOD_RangeSize(proj.current = cacr.binproj,
                               proj.future = cacr.binproj.ssp585)
cacr.rs585


load("Modeling/Final-Projections/Castilleja.cryptantha.ssp245EM.final.ensemble.projection.out")

cacr.emproj.ssp245 <- rast("Modeling/Final-Projections/proj_ssp245EM.final_Castilleja.cryptantha_ensemble.tif")
plot(cacr.emproj.ssp245)



### All Plots Together ---------

plot_grid(pera.plot + theme(legend.position = "none"),
          tast.plot + theme(legend.position = "none"),
          cacr.plot + theme(legend.position = "none"),
          clme.plot + theme(legend.position = "none"),
          orth.plot + theme(legend.position = "none"),
          anni.plot + theme(legend.position = "none"),
          chth.plot + theme(legend.position = "none"),
          locu.plot + theme(legend.position = "none"),
          get_legend(clme.plot),
          nrow = 3, byrow = TRUE)







