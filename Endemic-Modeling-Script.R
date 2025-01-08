##### WASHINGTON ENDEMIC MODELING PROJECT #####
#### Gjording et al.

# 1.0 LOAD ITEMS -----------------------------------------------------------

# 1.1 Load Packages ----

install.packages("biomod2")
install.packages("tidyverse")
install.packages("rinat")
#install.packages("rgbif")
#install.packages("rvest")

library(biomod2)
library(tidyverse)
library(rinat)
#library(rgbif)
#library(rvest)

# 1.2 Load Functions ----

# 1.3 Load Data ----

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

# 2.2 Adjusting Predictor Variables ----

# 2.3 Principle Coordinate Analysis (PCA) ----

# 3.0 Species Distribution Model (SDM) -------------------------------------

# 4.0 Figures --------------------------------------------------------------
