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
         type == "Point" | type == "Point fragment" | type == "Preform", 
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") 

#fill the NA fields in the munsell hue column to mark them as colourless/translucent material
Points.xrf[8][is.na(Points.xrf[8])] <- "Colourless"

#filter XRF data on sample size (only include samples with smallest dimension >= 10mm), 
#select elements to include in PCA, as well as columns to group by
xrf <-
  Points.xrf %>%
  filter(max_length_mm >= 10 & max_width_mm >= 10) %>% 
  dplyr::select(reading_no, site_id, sample_id, hue, material, `Mg`, `Al`, `Si`, `K`, `Ca`, `Fe`, `Sr`, `Zr`, `Ti`)


#impute missing data using least trimmed squares regression
xrf.imp <- 
  impAll(xrf[,6:14])

#transform data using centered log-ratio 
xrf.clr <- 
  clr(xrf.imp, ifclose = TRUE)

#perform pca with mean centering and unit variance scaling
xrf.pca <- 
  prcomp(xrf.clr, center = TRUE, scale.=TRUE)

#Determine no. of clusters with elbow graph
fig.1 <- 
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

#Extract loadings from xrf.pca
loadings <- as.data.frame(xrf.pca$rotation[,1:3])
#As elements are stored as rownames, set them as first column so it can be used as a variable
loadings <- setDT(loadings, keep.rownames = "Element")[]
#Melt loadings from wide to long format
loadings.melted <- melt(loadings, id.vars="Element")
#Set names of columns
loadings <- setNames(loadings.melted, c("Element", "PC", "Value"))

#bar graph of the loadings (P1-P2) for the XRF data
fig.2 <- 
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
    #coord_flip()

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(fig.1,fig.2, 
                    ncol = 1, 
                    nrow = 2,
                    labels = c("A", "B"))


ggsave("009-xrf-pca-load-scree.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=25, 
       height=20,
       units = "cm",
       dpi = 300)