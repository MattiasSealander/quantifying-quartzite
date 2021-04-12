library(tidyverse)
library(prospectr)
library(FactoMineR)
library(factoextra)

#Read in raw NIR measurements with meatadata
NIR <- as.data.frame(read.csv2("analysis/data/raw_data//NIR_raw_data.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Prepare filtering
artefact <- c("Point", "Point fragment", "Preform")
material <- c("Quartzite", "Breccie quartz", "Breccie quartzite", "Felsitic porphyric quartz", "Quartz")

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

vn <- prcomp(vis_nir[,c(2:302)], center = TRUE, scale = FALSE)

#centered  PCA
VNIR.pca <- PCA(vis_nir_center[,c(2:302)], graph = FALSE, scale.unit = FALSE)

#centered and SNV normalized  PCA
VNIR_snv.pca <- PCA(standardNormalVariate(vis_nir_center[,c(2:302)]), graph = FALSE)

#Screeplot showing eigenvalues of visNIR PCA
#dashed line corresponds to value 1
fviz_eig(VNIR.pca,
         geom = "bar",
         choice = "eigenvalue",
         ncp = 3,
         barfill = "black",
         barcolor = "black",
         title = "Eigenvalues of PCA",
         xlab = "PCs") +
  geom_hline(yintercept = 1, col="red", lty=5) +
  theme_minimal(base_size = 10)

#Screeplot of contributing variables to PC1
#dashed reference line corresponds to the expected value if the contribution was uniform.
fviz_contrib(VNIR.pca,
             choice = "var",
             fill = "black",
             color = "black",
             axes = 1,
             top = 60,
             title = "Contribution of variables to PC1") +
        theme_minimal(base_size = 20) +
        theme(axis.title.x=element_blank(), text = element_text(size=10),
                axis.text.x = element_text(angle=45, hjust=1))


#Screeplot of contributing variables to PC2
#dashed reference line corresponds to the expected value if the contribution was uniform.
fviz_contrib(VNIR.pca,
             choice = "var",
             fill = "black",
             color = "black",
             axes = 2,
             top = 60,
             title = "Contribution of variables to PC2")+
       theme_minimal(base_size = 20) +
       theme(axis.title.x=element_blank(), text = element_text(size=10),
             axis.text.x = element_text(angle=45, hjust=1))

#Screeplot of contributing variables to PC2
#dashed reference line corresponds to the expected value if the contribution was uniform.
fviz_contrib(VNIR.pca,
             choice = "var",
             fill = "black",
             color = "black",
             axes = 3,
             top = 60,
             title = "Contribution of variables to PC3")+
      theme_minimal(base_size = 20) +
      theme(axis.title.x=element_blank(), text = element_text(size=10),
            axis.text.x = element_text(angle=45, hjust=1))


#PCA of centered data
fviz_pca_ind(VNIR.pca, axes = c(1, 2),
                geom = c("point"),
                geom.ind = c("point"), # show points only (but not "text")
                label = "var",
                palette = c("#BC3C29FF", "#0072B5FF", "#EE4C97FF", "#E18727FF"),
                col.ind = Points$hue, #color by sample hue
                fill.ind = Points$hue, #fill point by sample hue
                addEllipses = FALSE,
                repel = TRUE,
                obs.scale = 1, var.scale = 1,
                legend.title = "Hue",
                title = "PCA - Quartzite points") +
      theme_minimal(base_size = 20) +
      scale_shape_manual(values=c(21,22,23,24)) +
      labs(x = "PC1 (98%)",
           y = "PC2 (1.9%)")




