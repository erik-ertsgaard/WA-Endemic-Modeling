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
