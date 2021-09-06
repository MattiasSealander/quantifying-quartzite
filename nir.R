library(tidyverse)
library(ggplot2)
library(ggbiplot)
library(prospectr)
library(patchwork)
library(factoextra)

#Import raman data
nir.csv <-
  read.csv2("C:/Users//masj0062/Documents/DOKTORANDSPROJEKT/VERSION CONTROL/MethodPaper/analysis/data/derived_data//NIR_master_averaged_20210827.csv", sep = ";",
            dec = ",", header = TRUE, na = c("","NA","NULL",NULL))

#Filter to focus on quartz/quartzite material that are points and preforms
Points.nir <-
  nir.csv %>%
  filter(material == "Quartz" | material == "Quartzite" | material == "Breccie quartz" | material == "Breccie quartzite" | material  == "Felsitic porphyric quartz",
         type == "Point" | type == "Point fragment" | type == "Preform")

#Mark the samples that lack a munsell hue as transclucent
Points.nir[8][is.na(Points[8])] <- "Translucent"

#perform PCA with SNV normalization and mean-center
nir.pca <-
  prcomp(Points.nir[,c(678:2178)], center = TRUE, scale = FALSE)

#prepare a basic score plot with fviz_pca_ind using PC1 and PC2
basic_plot1 <-
  fviz_pca_ind(nir.pca, axes = c(1,2), label="none")

#prepare labels for PCs
pc1var <- round(summary(nir.pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(nir.pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")


#prepare color/fill and symbols
colors <- c("#E31CFF", "#718EFF", "#0080FF", "#D5DF9E", "#C6DF39", "#FFFFC1", "#FFFFA1", "#FFFF00", "#FFFF38", "#FFB468", "#000000", "#202020", "#404040", "#606060", "#808080", "#9F9F9F", "#DFDFDF", "#FFFFFF")
shapes <- c(24,25,23,22,21)

#bind the basic fviz plot for PC 1 and 2, and use as basis for a more customizeable plot in ggpplot
ggplot(cbind(basic_plot1$data, Points[, c(8,10)]),
             aes(x=x, y=y, shape = material, fill = munsell_hue)) +
  geom_point(size=3) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = pc1lab,
       y = pc2lab) +
  scale_shape_manual(name = "Material",
                     values=shapes) +
  scale_fill_manual(name = "Munsell hue",
                    breaks = c("Translucent", "2.5PB", "10B", "2.5GY", "5GY", "2.5Y", "5Y", "7.5Y", "10Y", "5YR", "N1", "N2", "N3", "N4", "N5", "N6", "N8", "N9"),
                    values=colors) +
  guides(shape = guide_legend(override.aes = list(fill = "black"), ncol=2, title.position="top", title.hjust = 0.5),
         fill = guide_legend(override.aes = list(shape = 22,
                                                 fill = colors,
                                                 color = "black",
                                                 size=3),
                             title.position="top", title.hjust = 0.5))+
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black", vjust = - 10, hjust = 0.02),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")
