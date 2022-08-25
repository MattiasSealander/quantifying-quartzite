suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(rgr))
suppressPackageStartupMessages(library(robCompositions))
suppressPackageStartupMessages(library(tidyverse))

#Import xrf data, set empty fields to NA
xrf.csv <-
  read.csv2("./analysis/data/raw_data/XRF/xrf_quantitative_data_20220407.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL"))

#Import descriptive metadata
metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata_20220510.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#merge XRF data with metadata
xrf.merged <- 
  as.data.frame(merge(metadata.csv, xrf.csv, by='sample_id'))

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
Points.xrf <-
  xrf.merged %>%
  filter(site_id == "Vilhelmina 1069" | site_id == "Vilhelmina 109" | site_id == "Vilhelmina 112" | site_id == "Vilhelmina 1124" | site_id == "Vilhelmina 1127" | site_id == "Vilhelmina 114" |
           site_id == "Vilhelmina 115" | site_id == "Vilhelmina 117" | site_id == "Vilhelmina 118" | site_id == "Vilhelmina 1254" | site_id == "Vilhelmina 216" | site_id == "Vilhelmina 235" |
           site_id == "Vilhelmina 240" | site_id == "Vilhelmina 245" | site_id == "Vilhelmina 252" | site_id == "Vilhelmina 263" | site_id == "Vilhelmina 335" | site_id == "Vilhelmina 356" |
           site_id == "Vilhelmina 399" | site_id == "Vilhelmina 411" | site_id == "Vilhelmina 419" | site_id == "Vilhelmina 439" | site_id == "Vilhelmina 444" | site_id == "Vilhelmina 450" |
           site_id == "Vilhelmina 458" | site_id == "Vilhelmina 539" | site_id == "Vilhelmina 542" | site_id == "Vilhelmina 611" | site_id == "Vilhelmina 619" | site_id == "Vilhelmina 636" | 
           site_id == "Vilhelmina 637" | site_id == "Vilhelmina 643" | site_id == "Vilhelmina 769" | site_id == "Vilhelmina 949" | site_id == "Vilhelmina 95" | site_id == "Åsele 101" | 
           site_id == "Åsele 107" | site_id == "Åsele 115" | site_id == "Åsele 117" | site_id == "Åsele 119" | site_id == "Åsele 129" | site_id == "Åsele 182" | site_id == "Åsele 188" |
           site_id == "Åsele 393" | site_id == "Åsele 56" | site_id == "Åsele 91" | site_id == "Åsele 92" | site_id == "Åsele 99", 
         !sample_id %in% c('183', '193', '259', '359', '400'),
         type == "Point" | type == "Point fragment" | type == "Preform", 
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") 

#fill the NA fields in the munsell hue column to mark them as colourless/translucent material
Points.xrf[8][is.na(Points.xrf[8])] <- "Colourless"

#filter XRF data on sample size (only include samples with smallest dimension >= 10mm), 
#select elements to include in PCA, as well as columns to group by
xrf <-
  Points.xrf %>%
  filter(max_length_mm >= 10 & max_width_mm >= 10) %>% 
  dplyr::select(reading_no, site_id, sample_id, hue, material, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`)

#impute missing data using least trimmed squares regression
xrf.imp <- 
  impAll(xrf[,6:12])

#transform data using centered log-ratio 
xrf.clr <- 
  clr(xrf.imp, ifclose = TRUE)

#perform pca with mean centering and unit variance scaling
xrf.pca <- 
  prcomp(xrf.clr, center = TRUE, scale.=TRUE)

#prepare labels for PCs
xrf.pc1var <- round(summary(xrf.pca)$importance[2,1]*100, digits=2)
xrf.pc2var <- round(summary(xrf.pca)$importance[2,2]*100, digits=2)
xrf.pc3var <- round(summary(xrf.pca)$importance[2,3]*100, digits=2)
xrf.pc4var <- round(summary(xrf.pca)$importance[2,4]*100, digits=2)
xrf.pc1lab <- paste0("PC1 (",as.character(xrf.pc1var),"%)")
xrf.pc2lab <- paste0("PC2 (",as.character(xrf.pc2var),"%)")
xrf.pc3lab <- paste0("PC3 (",as.character(xrf.pc3var),"%)")
xrf.pc4lab <- paste0("PC4 (",as.character(xrf.pc4var),"%)")


#prepare color/fill and symbols for score plots
pca.colors <- c("#9467BDFF", "#7F7F7FFF", "#FF7F0EFF", "#1F77B4FF")
pca.hue <- c("Colourless", "Dark", "Light", "White")

#Extract loadings from xrf.pca
loadings <- as.data.frame(xrf.pca$rotation[,1:4])
#As elements are stored as rownames, set them as first column so it can be used as a variable
loadings <- setDT(loadings, keep.rownames = "Element")[]
#Melt loadings from wide to long format
loadings.melted <- melt(loadings, id.vars="Element")
#Set names of columns
loadings <- setNames(loadings.melted, c("Element", "PC", "Value"))

#Show eigenvalue scree plot
#fviz_eig(xrf.pca, choice = "eigenvalue")

#keep PCs above eigenvalue 1
xrf.transform = as.data.frame(-xrf.pca$x[,1:4])

#Determine no. of clusters with elbow graph
#set.seed(100)
#fviz_nbclust(xrf.transform, kmeans, method = 'wss')

#Perform cluster analysis on the Principal Components
kmeans.xrf <- 
  kmeans(xrf.transform, centers = 3, nstart = 50)

#Prepare result of cluster analysis, so that they can be used in biplot
xrf$cluster <- 
  factor(kmeans.xrf$cluster)

#pivot data so boxplot can be generated for each cluster and element
xrf_long <- xrf %>% 
  dplyr::select(cluster, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`) %>% 
  pivot_longer(-cluster, names_to = "variable", values_to = "value")

#Prepare for facet labels
Elements <- c(
  'Al' = "Al",
  'Si' = "Si",
  'S' = "S",
  'K' = "K",
  'Ca' ="Ca",
  'Fe' = "Fe",
  'Zr' = "Zr",
  'Ti' = "Ti"
)  

#Determine no. PC scree plot
load.1 <- 
  fviz_eig(xrf.pca, 
           choice = "eigenvalue",
           barfill = "#999999",
           barcolor = "black",
           linecolor = "black",
           ncp = 10,
           bar_width=0.5) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")


#bar graph of the loadings (P1-P2) for the XRF data
load.2 <- 
  loadings %>%
  ggplot(aes(x=Element, y=Value)) +
  geom_bar(
    stat = "identity", position = "identity",
    color = "black", fill = "#999999",
    width = 0.5) +
  facet_wrap( ~ PC, scales = "free") +
  theme_bw() +
  theme(plot.margin = margin(10,10,10,10),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10, face="bold", colour = "black"),
        axis.text.x = element_text(size=10, face="bold", colour = "black"),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size=14, face="bold", colour="black"))


#Layout the plots in one figure
fig1 <-
  ggpubr::ggarrange(load.1,load.2, 
                    ncol = 1, 
                    nrow = 2,
                    labels = c("A", "B"))


ggsave("009-xrf-pca-load-screeNOSR.png",
       fig1,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=25, 
       height=20,
       units = "cm",
       dpi = 300)


#Biplot with ellipses showing the results of the k-means cluster analysis 
fig2a <- 
  fviz_pca_biplot(xrf.pca,
                axes = c(1,2),
                geom = "text",
                geom.var = c("point", "text"),
                label = "var",
                labelsize = 5,
                col.ind = xrf$hue,
                col.var = "black",
                fill.var = "red", 
                repel = TRUE,
                pointshape = 21, 
                pointsize = 4) +
  theme_bw() +
  labs(x = xrf.pc1lab,
       y = xrf.pc2lab) +
  geom_text(aes(label=xrf$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(aes(fill = xrf$hue, 
                 shape = factor(xrf$material)), 
             size = 3) +
  scale_shape_manual(name = "Material", 
                     values=c(24,22,21),
                     guide = guide_legend(override.aes = list(fill = "black"),
                                          title.position="top",
                                          title.hjust = 0.5,
                                          order = 1)) +
  scale_color_manual(values = pca.colors,
                    guide = "none") +
  scale_fill_manual(name = "Hue",
                    values = pca.colors,
                    guide = guide_legend(override.aes = list(shape = 21,
                                                             fill = pca.colors,
                                                             color = "black",
                                                             size=3),
                                         title.position="top",
                                         title.hjust = 0.5,
                                         order = 2)) +
  ggnewscale::new_scale_fill() +
  stat_chull(aes(fill = xrf$cluster),
             alpha = 0.3, 
             geom = "polygon") + 
  scale_fill_manual(name = "Cluster", 
                    values = c( "#0072B2", "#D55E00", "#CC79A7"),
                    guide_legend(order = 1)) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

#Biplot with ellipses showing the results of the k-means cluster analysis 
fig2b <- 
  fviz_pca_biplot(xrf.pca,
                  axes = c(3,4),
                  geom = "text",
                  geom.var = c("point", "text"),
                  label = "var",
                  labelsize = 5,
                  col.ind = xrf$hue,
                  col.var = "black",
                  fill.var = "red", 
                  repel = TRUE,
                  pointshape = 21, 
                  pointsize = 4) +
  theme_bw() +
  labs(x = xrf.pc3lab,
       y = xrf.pc4lab) +
  geom_text(aes(label=xrf$sample_id, hjust=0.5,vjust=-1.0)) +
  geom_point(aes(fill = xrf$hue, 
                 shape = factor(xrf$material)), 
             size = 3) +
  scale_shape_manual(name = "Material", 
                     values=c(24,22,21),
                     guide = guide_legend(override.aes = list(fill = "black"),
                                          title.position="top",
                                          title.hjust = 0.5,
                                          order = 1)) +
  scale_color_manual(values = pca.colors,
                     guide = "none") +
  scale_fill_manual(name = "Hue",
                    values = pca.colors,
                    guide = guide_legend(override.aes = list(shape = 21,
                                                             fill = pca.colors,
                                                             color = "black",
                                                             size=3),
                                         title.position="top",
                                         title.hjust = 0.5,
                                         order = 2)) +
  ggnewscale::new_scale_fill() +
  stat_chull(aes(fill = xrf$cluster),
             alpha = 0.3, 
             geom = "polygon") + 
  scale_fill_manual(name = "Cluster", 
                    values = c( "#0072B2", "#D55E00", "#CC79A7"),
                    guide_legend(order = 1)) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

#Save score plot for PC1-PC2
ggsave("009-xrf-pca-aNOSR.png",
       fig2a,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)
#Save score plot for PC2-PC3
ggsave("009-xrf-pca-bNOSR.png",
       fig2b,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)

#Boxplot of XRF raw data by element and K-means cluster
fig3 <- 
  ggplot(xrf_long, aes(x=variable, y=value, fill = cluster)) +
  geom_boxplot() +  
  facet_wrap( ~ variable, scales = "free", labeller = labeller(Element = Elements)) +
  #facetted_pos_scales(y = scales) +
  theme_bw() +
  scale_fill_manual(name = "Cluster", 
                    values = c("#0072B2", "#D55E00", "#CC79A7")) +
  theme(
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(size=10, face = "bold", colour = "black", margin = margin(t = 10, b = -20), vjust=0.02, hjust=0.01),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())

#Save K-means boxplot
ggsave("009-xrf-kmeans-boxplotNOSR.png",
       fig3,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)

