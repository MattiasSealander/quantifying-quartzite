library(tidyverse)
library(prospectr)

#Read in raw NIR measurements with meatadata
NIR <- as.data.frame(read.csv2("analysis/data/raw_data//NIR_raw_data.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Prepare filtering
artefact <- c("Point", "Point fragment", "Preform")
material <- c("Quartzite", "Breccie quartz", "Breccie quartzite")

#Filter data to focus on quartzite points, point fragments and preforms
Points <-
  NIR %>%
    filter(material %in% material,
           type %in% artefact)

#Select visible range of data
vis_nir <-
  Points %>%
    select(unique_id, X400.0:X700.0)

#mean center measurements
vis_nir_center <-
  vis_nir %>%
    mutate(across(c(2:302), scale, scale = FALSE))


