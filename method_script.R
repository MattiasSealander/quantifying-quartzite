library(dplyr)
library(tidyverse)
library(prospectr)
library(FactoMineR)
library(factoextra)
library(ChemoSpec)
library(MASS)
library(reshape2)
library(cowplot)
library(ggplot2)
library(ggsci)

#Read in raw NIR measurements with meatadata
NIR <- as.data.frame(read.csv2("analysis/data/derived_data//nir_mean_20210507.csv", sep = ";", dec = ",", header = TRUE, encoding = "UTF-8"))

#Read in raw NIR measurements for 1000-2500 nm range
csv <- read.csv2("analysis/data/raw_data//NIR_only_raw_data_1000_2500.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8")

#make first column row names
csv2<- csv %>%
  remove_rownames %>%
  column_to_rownames(var="unique_id") %>%
  as.data.frame()

#SNV correction for each row(frequency)
snv <- standardNormalVariate(csv2)

#Write csv file for SNV treated data
write.csv2(Points, "a.csv", fileEncoding = "UTF-8")

#Prepare filters for data
#Didnt work as intended, missed some scrapers etc.
artefact <- c("Point", "Point fragment", "Preform")
material <- c("Quartzite", "Breccie quartz", "Breccie quartzite", "Felsitic porphyric quartz", "Quartz")

#Filter data to focus on quartzite points, point fragments and preforms made from quartz/quartzite materials
Points <-
  NIR %>%
    filter(material == "Quartz" | material == "Quartzite" | material == "Breccie quartz" | material == "Breccie quartzite" | material  == "Felsitic porphyric quartz",
           type == "Point" | type == "Point fragment" | type == "Preform")


#If average is supposed to be done on repeated measurements
#NIR_averaged <- aggregate(NIR[, 26:2176], list(NIR$Sample_ID), mean)

#Make spectra object between 400-2500 nm
spec_400 <- matrix2SpectraObject(gr.crit = c("Dark", "Light", "Translucent", "White"),
                             gr.cols = c("auto"),
                             freq.unit = "Wavelength nm",
                             int.unit = "Absorbance",
                             descrip = "NIR measurements on quartz/quartzite points between 400-2500 nm",
                             in.file = "analysis/data/raw_data//NIR_only_raw_data.csv",
                             chk = TRUE,
                             dec = ".",
                             sep = ";"
)

#Make spectra object between 1000-2500 nm
spec_1000 <- matrix2SpectraObject(gr.crit = c("Dark", "Light", "Translucent", "White"),
                              gr.cols = c("auto"),
                              freq.unit = "Wavelength nm",
                              int.unit = "Absorbance",
                              descrip = "NIR measurements on quartz/quartzite points between 100-2500 nm",
                              in.file = "analysis/data/raw_data//NIR_only_raw_data_1000_2500.csv",
                              chk = TRUE,
                              dec = ".",
                              sep = ";"
)


#Make spectra object between 1000-2500 nm which has been SNV corrected
spec_snv <- matrix2SpectraObject(gr.crit = c("Dark", "Light", "Translucent", "White"),
                             gr.cols = c("auto"),
                             freq.unit = "Wavelength nm",
                             int.unit = "Absorbance",
                             descrip = "NIR measurements on quartz/quartzite points between 100-2500 nm",
                             in.file = "analysis/data/derived_data//NIR_snv_1000_2500.csv",
                             chk = TRUE,
                             dec = ",",
                             sep = ";"
)


#inspect data
sumSpectra(spec_snv)

#normalization using TotInt (y is scaled to 1)
res1 <- normSpectra(spec_1000, method = "TotInt")

#compare normalized spectra to raw
plotSpectra(spec_1000, which = c(1:50))

plotSpectra(spec_snv, which = c(1:50))

#Calculate PCA using mean centering and autoscale
pca <-
  c_pcaSpectra(spec_1000, choice = "noscale", cent = TRUE)

#Calculate PCA on normalized data using mean centering and no scaling
pca_snv <-
  c_pcaSpectra(spec_snv, choice = "Pareto", cent = TRUE)

#Screeplot of PCs
plotScree(pca_snv, style = "alt")

#Plot the loadings of PC1 to PC3 with sample 1 as reference
plotLoadings(spec2, pca1, loads = c(1,2,3), ref = 10)

#Plot the PCA without SNV
plotScores(spec_1000,
           pca,
           pcs = c(1,2),
           ellipse = "none",
           tol = "none")

#Plot the PCA with SNV
plotScores(spec_snv,
           pca_snv,
           pcs = c(1,2),
           ellipse = "none",
           tol = "none")




#Select visible range of data
vis_nir <-
  Points %>%
    select(sample_id, X400.0:X750.0) %>%
      remove_rownames %>%
      column_to_rownames(var="sample_id") %>%
      as.data.frame()

VNIR_snv <- standardNormalVariate(vis_nir)

#mean center measurements
vis_nir_center <-
  vis_nir %>%
    mutate(across(c(1:351), scale, scale = FALSE))

#centered  PCA
VNIR.pca <- PCA(vis_nir_center, graph = FALSE, scale.unit = FALSE)

#centered and SNV normalized  PCA
VNIR_snv.pca <- PCA(VNIR_snv_center, graph = FALSE)

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
windows();fviz_pca_ind(VNIR.pca, axes = c(1, 2),
                geom = c("point"),
                geom.ind = c("point"), # show points only (but not "text")
                label = "var",
                palette = c("black", "#0072B5FF", "forestgreen", "#E18727FF"),
                col.ind = Points$hue, #color by sample hue
                fill.ind = Points$hue, #fill point by sample hue
                addEllipses = TRUE,
                repel = TRUE,
                obs.scale = 1, var.scale = 1,
                legend.title = "Hue",
                pointsize = 5,
                title = "PCA - Quartzite points",
                invisible = "quali") +
      theme_minimal(base_size = 20) +
      scale_shape_manual(values=c(3,22,24,21)) +
      scale_size_manual(values=c(5,3,3,3)) +
      labs(x = "PC1 (98.4%)",
           y = "PC2 (1.5%)")


# extract pc scores for first two component and add to dat dataframe
Points$pc1 <- VNIR.pca$ind$coord[, 1] # indexing the first column

Points$pc2 <- VNIR.pca$ind$coord[, 2]  # indexing the second column


#We also need to extract the data for the variable contributions to each of the pc axes.
pca.vars <- VNIR.pca$var$coord %>% data.frame

pca.vars$vars <- rownames(pca.vars)

pca.vars.m <- melt(pca.vars, id.vars = "vars")


#By convention, the variable contribution plot has a circle around the variables that has a radius of 1. Here’s some code to make one.
circleFun <- function(center = c(0,0),diameter = 1, npoints = 100){
  r = diameter / 2
  tt <- seq(0,2*pi,length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

circ <- circleFun(c(0,0),2,npoints = 500)


p1 = ggplot(data = Points,
                 aes(x = pc1,
                     y = pc2,
                     color = material,
                     shape = hue)) +
          geom_hline(yintercept = 0, lty = 2) +
          geom_vline(xintercept = 0, lty = 2) +
          guides(color = guide_legend(title = "Material"), shape = guide_legend(title = "Hue")) +
          geom_point(alpha = 0.8, size = 3) +
          scale_shape_manual(values=c(3,15,17,16)) +
          scale_color_nejm()

#pass two plots as objects and then arrange them in one figure
p1_nejm = p1 + scale_color_nejm()
p2_nejm = p2 + scale_fill_nejm()
grid.arrange(p1_nejm, p2_nejm, ncol = 2)








##NIR PCA OF 1000-2500 nm RANGE
#Select NIR range of data
nir_1000 <-
  Points %>%
  select(unique_id, X750.0:X2500.0) %>%
  remove_rownames %>%
  column_to_rownames(var="unique_id") %>%
  as.data.frame()

nir_snv <- standardNormalVariate(nir)

#mean center measurements
nir_snv1 <-
  nir_snv %>%
  mutate(across(c(1:1751), scale, scale = FALSE))


#Do PCA
NIR_snv.pca <- PCA(nir_snv1, graph = FALSE)


#Screeplot showing eigenvalues of visNIR PCA
#dashed line corresponds to value 1
fviz_eig(NIR_snv.pca,
         geom = "bar",
         choice = "eigenvalue",
         ncp = 5,
         barfill = "black",
         barcolor = "black",
         title = "Eigenvalues of PCA",
         xlab = "PCs") +
  geom_hline(yintercept = 1, col="red", lty=5) +
  theme_minimal(base_size = 10)

#dashed reference line corresponds to the expected value if the contribution was uniform.
fviz_contrib(NIR_snv.pca,
             choice = "var",
             fill = "black",
             color = "black",
             axes = 1,
             top = 100,
             title = "Contribution of variables to PC1") +
  theme_minimal(base_size = 20) +
  theme(axis.title.x=element_blank(), text = element_text(size=10),
        axis.text.x = element_text(angle=-45, vjust=1)) +
  coord_flip()

#PCA of centered data
fviz_pca_ind(NIR_snv.pca, axes = c(1, 2),
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
