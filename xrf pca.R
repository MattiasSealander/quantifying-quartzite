library(tidyverse)
library(ggplot2)
library(missMDA)
library(ggbiplot)
library(factoextra)



#Import xrf data
XRF.csv <-
  read.csv2("C:/Users//masj0062/Documents/DOKTORANDSPROJEKT/VERSION CONTROL/MethodPaper/analysis/data/raw_data//XRF/XRF_quantitative.csv", sep = ";",
            dec = ".", header = TRUE)

Points <-
  XRF.csv %>%
  filter(material == "Quartz" | material == "Quartzite" | material == "Breccie quartz" | material == "Breccie quartzite" | material  == "Felsitic porphyric quartz",
         type == "Point" | type == "Point fragment" | type == "Preform")

#select elements to include in pca
xrf <-
  Points %>%
    select(reading_no, sample_id, munsell_hue, material, X12Mg, X13Al, X14Si, X15P, X16S, X19K, X20Ca, X23V, X25Mn, X26Fe, X30Zn, X38Sr, X39Y, X40Zr, X56Ba, X81Ti)

#log transform data
xrf.log <- log(xrf[,5:20]+1)

#impute missing vales with scale
xrf.imp <- imputePCA(xrf.log, ncp = 2, scale = TRUE, method = "Regularized", maxiter = 1000)
#perform pca with mean centering
xrf.pca <- prcomp(xrf.imp$completeObs, center = TRUE, scale.=FALSE)

#prepare labels for PCs
pc1var <- round(summary(xrf.pca)$importance[2,1]*100, digits=1)
pc2var <- round(summary(xrf.pca)$importance[2,2]*100, digits=1)
pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")

#prepare color and symbols
colors <- c("#E31CFF", "#0080FF", "#FFFF00", "#D5DF9E", "#718EFF", "#FFFFC1", "#C6DF39", "#FFFFC1", "#FFB468", "#FFFF38", "#000000", "#202020", "#404040", "#606060", "#808080", "#9F9F9F", "#DFDFDF", "#FFFFFF")
shapes <- c(24,25,23,22,21)

#score plot, need to specify to use ggbiplot package
ggbiplot::ggbiplot(xrf.pca,
         obs.scale = 1,
         var.scale = 1,
         group=xrf$munsell_hue,
         varname.size = 4,
         labels.size=3,
         ellipse = FALSE,
         circle = FALSE) +
  theme_minimal() +
  geom_point(aes(shape=xrf$material, fill = xrf$munsell_hue), size = 4) +
  scale_shape_manual(name = "Material", values=shapes) +
  scale_color_manual(name = "Munsell hue", values=colors) +
  scale_fill_manual(values=colors) +
  guides(fill = "none",
         shape = guide_legend(override.aes = list(fill = "black"),
                              title.position="top",
                              title.hjust = 0.5),
         color = guide_legend(override.aes = list(shape = 22,
                                                  fill = colors,
                                                  color = "black",
                                                  size=3))) +
  theme(legend.background = element_rect(linetype = "solid", color = "black"))
