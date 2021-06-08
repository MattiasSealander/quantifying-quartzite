library(tidyverse)
library(ggsci)

#Read in metadata file for points
points <- as.data.frame(read.csv2("analysis/data/raw_data//Samples_GIS_file.csv", sep = ";", dec = ",", header = TRUE))


#Prepare filters for data
artefact <- c("Point", "Point fragment", "Preform")
material <- c("Quartzite", "Breccie quartz", "Breccie quartzite", "Felsitic porphyric quartz", "Quartz")


#Filter data to focus on quartzite points, point fragments and preforms made from quartz/quartzite materials
Points <-
  points %>%
  filter(Material %in% material,
         Type %in% artefact)

# If you need to flip the order (because you've flipped the orientation)
# call position_stack() explicitly:
ggplot(Points, aes(y = Material)) +
  geom_bar(stat = "count", fill = "black") +
  geom_text(stat='count', aes(label=..count..), vjust = -0.5) +
  theme_minimal() +
  ggtitle("Analyzed Points by Material") +
  xlab("Count") +
  theme(plot.title = element_text(size = 22,
                                  hjust = 0.5,
                                  face = "bold"),
        axis.title = element_text(size = 14,
                                  face = "bold"),
        axis.text.x = element_text(face = "bold"),
        axis.title.x = element_blank(),
        legend.title = element_text(size = 12,
                                    face = "bold"),
        legend.position = "bottom") +
  coord_flip()
