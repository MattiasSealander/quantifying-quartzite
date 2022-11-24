suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(robCompositions))
suppressPackageStartupMessages(library(tidyverse))

#Import xrf data, set empty fields to NA
xrf.csv <-
  read.csv2("./analysis/data/raw_data/XRF/xrf_quantitative_data_20220407.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL"))

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2("./analysis/data/raw_data/NIR/asd_raw_data_20220407.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Import descriptive metadata
metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

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
         type == "Point" | type == "Point fragment" | type == "Preform", 
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") 

#fill the NA fields in the munsell hue column to mark them as colourless/translucent material
Points.xrf[8][is.na(Points.xrf[8])] <- "Colourless"

#filter XRF data on sample size (only include samples with smallest dimension >= 10mm), 
#select elements to include in PCA, as well as columns to group by
xrf <-
  Points.xrf %>%
  filter(max_length_mm >= 10 & max_width_mm >= 10) %>% 
  dplyr::select(reading_no, site_id, sample_id, hue, material, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`) %>% 
  `row.names<-`(., NULL) %>% 
  column_to_rownames(var = "reading_no")

#impute missing data using least trimmed squares regression
xrf.imp <- 
  impAll(xrf[,5:11])

#log-transform data due to it being percentages
xrf.log <- 
  log10(xrf.imp)

#if you want to transform data using centered log-ratio
#potential issues due to using geometric mean
#xrf.clr <- 
#  clr(xrf.imp, ifclose = TRUE)

#perform pca with mean centering and unit variance scaling
xrf.pca <- 
  prcomp(xrf.log, center = TRUE, scale.=TRUE)

#prepare labels for PCs
xrf.pc1var <- round(summary(xrf.pca)$importance[2,1]*100, digits=2)
xrf.pc2var <- round(summary(xrf.pca)$importance[2,2]*100, digits=2)
xrf.pc1lab <- paste0("PC1 (",as.character(xrf.pc1var),"%)")
xrf.pc2lab <- paste0("PC2 (",as.character(xrf.pc2var),"%)")

#prepare color/fill and symbols for score plots
pca.colors <- c("#9467BDFF", "#7F7F7FFF", "#FF7F0EFF", "#1F77B4FF")
pca.hue <- c("Colourless", "Dark", "Light", "White")

#Extract loadings from xrf.pca
loadings <- as.data.frame(xrf.pca$rotation[,1:2])
#As elements are stored as rownames, set them as first column so it can be used as a variable
loadings <- setDT(loadings, keep.rownames = "Element")[]
#Melt loadings from wide to long format
loadings.melted <- melt(loadings, id.vars="Element")
#Set names of columns
loadings <- setNames(loadings.melted, c("Element", "PC", "Value"))

#Show eigenvalue scree plot
#fviz_eig(xrf.pca, choice = "eigenvalue")

#keep PCs above eigenvalue 1
xrf.transform = as.data.frame(-xrf.pca$x[,1:2])

#Determine no. of clusters with elbow graph
#set.seed(100)
#fviz_nbclust(xrf.transform, kmeans, method = 'wss')

#Perform cluster analysis on the Principal Components
#In order to ensure that the cluster numbering is consistent between runs kmeans is calculated one time initially,
#center result is then ordered and passed to the kmeans in the second run. see: https://stackoverflow.com/questions/39906180/consistent-cluster-order-with-kmeans-in-r
centers <- kmeans(xrf.transform, centers = 3, nstart = 50)$centers

centers <- centers[order(-centers[,1], centers[,2]), ]

kmeans.xrf <- 
  kmeans(xrf.transform, centers = centers, nstart = 50)

#Prepare result of cluster analysis, so that they can be used in biplot
#In order to guard against potential mistakes, set rownames as column and use that for joining clusters to xrf dataframe
cluster <- 
  as.data.frame(factor(kmeans.xrf$cluster)) %>%
  rename("cluster" = "factor(kmeans.xrf$cluster)") %>% 
  rownames_to_column(var = "reading_no")

#Set reading_no as column again join with kmeans cluster output
#Then group by kmeans cluster and calculate Si mean, concatenating the mean value with cluster for display in PCA legend
xrf <- 
  xrf %>% 
  rownames_to_column(var = "reading_no") %>% 
  inner_join(cluster, by = "reading_no") %>%
  group_by(cluster) %>% 
  mutate(Si_mean = round(mean(`Si`), 1)) %>% 
  mutate(Si_mean = paste0("mean Si (", Si_mean, " %)")) %>% 
  unite(legend, sep = " - ", cluster, Si_mean, remove = FALSE)

#pivot data and add to new dataframe, so a boxplot can be generated for each cluster and element
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


ggsave("009-xrf-pca-load-scree.png",
       fig1,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=25, 
       height=20,
       units = "cm",
       dpi = 300)


#Biplot with ellipses showing the results of the k-means cluster analysis 
fig2 <- 
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
  #geom_text(aes(label=xrf$sample_id, hjust=0.5,vjust=-1.0)) +
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
  stat_chull(aes(fill = xrf$legend),
             alpha = 0.3, 
             geom = "polygon") + 
  scale_fill_manual(name = "Cluster", 
                    values = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF"),
                    guide_legend(order = 1)) +
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

#Save score plot for PC1-PC2
ggsave("009-xrf-pca.png",
       fig2,
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
                    values = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF")) +
  theme(
    legend.title = element_text(size = 12, face = "bold", colour = "black"),
    strip.text.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.x = element_blank(),
    plot.title = element_text(size=10, face = "bold", colour = "black", margin = margin(t = 10, b = -20), vjust=0.02, hjust=0.01),
    axis.title.y=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank())

#Save K-means boxplot
ggsave("009-xrf-kmeans-boxplot.png",
       fig3,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)

#Add XRF k-means cluster to NIR data to visualise in NIR PCA
#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.csv, by='sample_id'))

#Filter NIR data to focus on points and preforms made from quartz/quartzite material, exclude observations on dark material with mainly noise
#then summarise measurements by group (sample)
Points.nir <-
  nir.merged %>%
  filter(site_id == "Vilhelmina 1069" | site_id == "Vilhelmina 109" | site_id == "Vilhelmina 112" | site_id == "Vilhelmina 1124" | site_id == "Vilhelmina 1127" | site_id == "Vilhelmina 114" |
           site_id == "Vilhelmina 115" | site_id == "Vilhelmina 117" | site_id == "Vilhelmina 118" | site_id == "Vilhelmina 1254" | site_id == "Vilhelmina 216" | site_id == "Vilhelmina 235" |
           site_id == "Vilhelmina 240" | site_id == "Vilhelmina 245" | site_id == "Vilhelmina 252" | site_id == "Vilhelmina 263" | site_id == "Vilhelmina 335" | site_id == "Vilhelmina 356" |
           site_id == "Vilhelmina 399" | site_id == "Vilhelmina 411" | site_id == "Vilhelmina 419" | site_id == "Vilhelmina 439" | site_id == "Vilhelmina 444" | site_id == "Vilhelmina 450" |
           site_id == "Vilhelmina 458" | site_id == "Vilhelmina 539" | site_id == "Vilhelmina 542" | site_id == "Vilhelmina 611" | site_id == "Vilhelmina 619" | site_id == "Vilhelmina 636" | 
           site_id == "Vilhelmina 637" | site_id == "Vilhelmina 643" | site_id == "Vilhelmina 769" | site_id == "Vilhelmina 949" | site_id == "Vilhelmina 95" | site_id == "Åsele 101" | 
           site_id == "Åsele 107" | site_id == "Åsele 115" | site_id == "Åsele 117" | site_id == "Åsele 119" | site_id == "Åsele 129" | site_id == "Åsele 182" | site_id == "Åsele 188" |
           site_id == "Åsele 393" | site_id == "Åsele 56" | site_id == "Åsele 91" | site_id == "Åsele 92" | site_id == "Åsele 99", 
         type == "Point" | type == "Point fragment" | type == "Preform", 
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>% 
  dplyr::filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                                  "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                                  "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) %>% 
  replace_na(list(munsell_hue = "Colourless")) %>% 
  group_by(across(sample_id:river)) %>% 
  dplyr::summarise(across(`350.0`:`2500.0`, mean), .groups = "drop") %>% 
  left_join(xrf[,c("sample_id", "cluster")], by = "sample_id")

#perform PCA with SNV normalization and mean-center
nir.pca <-
  prcomp(Points.nir[,c(681:2180)], center = TRUE, scale = FALSE)

nir.pc1var <- round(summary(nir.pca)$importance[2,1]*100, digits=2)

#Prepare axis labels with variance in %
nir.pc1lab <- as.data.frame(paste0("PC1 (",as.character(round(summary(nir.pca)$importance[2,1]*100, digits=2)),"%)"))
nir.pc2lab <- as.data.frame(paste0("PC2 (",as.character(round(summary(nir.pca)$importance[2,2]*100, digits=2)),"%)"))
nir.pc3lab <- as.data.frame(paste0("PC3 (",as.character(round(summary(nir.pca)$importance[2,3]*100, digits=2)),"%)"))

#prepare a basic score plot with fviz_pca_ind using PC1 and PC2
basic_plot1 <-
  fviz_pca_ind(nir.pca, axes = c(1,2), label="none")

nir <- cbind(basic_plot1$data, Points.nir[, c(2181,10,1)])

#bind the basic fviz plot for PC 1 and 2, and use as basis for a more customizeable plot in ggpplot
fig4 <- 
  nir %>% 
  ggplot(aes(x=x, y=y, shape = material, fill = cluster)) + 
  geom_text(data=subset(nir, sample_id == 231 | sample_id == 404 | sample_id == 391), 
            aes(x,y,label=sample_id, hjust=-0.5,vjust=0.5)) +
  geom_point(size=4) +
  theme_bw() +
  ggtitle("Near infrared 1 000 - 2 500 nm") +
  labs(x = nir.pc1lab,
       y = nir.pc2lab) +
  scale_shape_manual(name = "Material",
                     values=c(24,22,21)) +
  scale_fill_manual(name = "XRF Cluster",
                    values=c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF"),
                    na.value="#D62728FF") +
  guides(shape = guide_legend(override.aes = list(fill = "black"), 
                              ncol = 3,
                              title.position="top", 
                              title.hjust = 0.5,
                              order = 2),
         fill = guide_legend(override.aes = list(shape = 22,
                                                 fill = c("#1F77B4FF", "#FF7F0EFF", "#2CA02CFF"), na.value="#D62728FF",
                                                 ncol = 4,
                                                 color = "black",
                                                 size=3,
                                                 order = 1),
                             title.position="top", title.hjust = 0.5))+
  theme(plot.title = element_text(size = 12, face = "bold", colour = "black", vjust = - 10, hjust = 0.02),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "bottom")

#Save NIR plot with kmeans cluster
ggsave("010-nir-pca-xrf-kmeans.png",
       fig4,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)