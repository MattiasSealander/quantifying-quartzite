# Loading package
library(gridBase)
library(spectrolab)
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata_20220510.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2("./analysis/data/raw_data/NIR/asd_raw_data_20220127.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
nir.averaged <- 
  aggregate(nir.csv[, 4:2154], list(sample_id = nir.csv$sample_id), mean)

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.averaged, by='sample_id'))

#Fill the NA fields in the munsell_hue column and mark them as colourless samples
nir.merged[8][is.na(nir.merged[8])] <- "Colourless"
Points.nir[,1] <- as.numeric(Points.nir[,1]) 

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
Points.nir <-
  nir.merged %>%
  dplyr::filter(site_id == "Vilhelmina 1069" | site_id == "Vilhelmina 109" | site_id == "Vilhelmina 112" | site_id == "Vilhelmina 1124" | site_id == "Vilhelmina 1127" | site_id == "Vilhelmina 114" |
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
                           "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56"))

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
Points.c <- Points.nir %>%
  dplyr::filter(hue == "Colourless") %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
Points.c <- 
  suppressWarnings(melt(setDT(Points.c), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with dark hues 
#and select the NIR range 1 000 - 2 500 nm
Points.d <- Points.nir %>%
  dplyr::filter(hue == "Dark") %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
Points.d <- 
  suppressWarnings(melt(setDT(Points.d), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with light hues 
#and select the NIR range 1 000 - 2 500 nm
Points.l <- Points.nir %>%
  dplyr::filter(hue == "Light") %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
Points.l <- 
  suppressWarnings(melt(setDT(Points.l), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on material with white hues 
#and select the NIR range 1 000 - 2 500 nm
Points.w <- Points.nir %>%
  dplyr::filter(hue == "White") %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
Points.w <- 
  suppressWarnings(melt(setDT(Points.w), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the colourless spectra
p.c <-
  ggplot(Points.c,aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Colourless",
                     values = "purple") + 
  theme_classic() +
  theme(legend.position = c(.9,.95),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the dark spectra
p.d <-
  ggplot(Points.d, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Dark",
                     values = "black") + 
  theme_classic() +
  theme(legend.position = c(.1,.95),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the light spectra
p.l <-
  ggplot(Points.l, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Light",
                     values = "blue") + 
  theme_classic() +
  theme(legend.position = c(.1,.95),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the white spectra
p.w <-
  ggplot(Points.w, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "White",
                     values = "red") + 
  theme_classic() +
  theme(legend.position = c(.1,.95),
        axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.c, p.d, p.l, p.w, 
                    ncol = 2, 
                    nrow = 2)

#Save figure
ggsave("002-nir-dark-filtered-summary.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)