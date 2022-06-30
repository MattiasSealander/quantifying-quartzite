suppressPackageStartupMessages(library(factoextra))
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
fig <- 
  fviz_eig(xrf.pca, 
           choice = "eigenvalue",
           barfill = "#999999",
           barcolor = "black",
           linecolor = "black",
           ncp = 10) +
  theme_classic() +
  theme(plot.title = element_blank(),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"),
        legend.background = element_rect(linetype = "solid", color = "black"),
        legend.position = "right")

ggsave("009-xrf-pca-scree.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       scale = 1, 
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)