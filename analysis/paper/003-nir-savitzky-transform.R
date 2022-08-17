#Load packages
suppressPackageStartupMessages(library(prospectr))
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

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
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
  filter(sample_id %in% c("54","58","152","170","171","173","185","186","187","188","197","208","221","227","241","249","253","254","258","259","269","280","378","384","385","386","391","401","403","404","408",
                          "432"))


#Filter NIR data to focus on material with dark hues 
#and select the NIR range 1 000 - 2 500 nm
Points.d <- Points.nir %>%
  filter(hue == "Dark") %>% 
  select(1, c(681:2180))

#Melt into long format
Points.d <- 
  suppressWarnings(melt(setDT(Points.d), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the dark spectra
p.d <-
  ggplot(Points.d, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Dark",
                     values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = c(.1,.95),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Filter NIR data to focus on material with dark hues 
#and select the NIR range 1 000 - 2 500 nm
Points.d258 <- Points.nir %>%
  filter(sample_id == "258") %>% 
  select(1, c(681:2180))

#Melt into long format
Points.d258.long <- 
  suppressWarnings(melt(setDT(Points.d258), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the dark spectra
p.d258 <-
  ggplot(Points.d258.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = "", group = sample_id), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Sample 258",
                     values = "black") + 
  #coord_cartesian(x=c(1000,2500), clip = "off") +
  #annotate("label", x = 890, y = 1.20, label = "B", label.size=NA, size= 5, fontface = "bold", fill = "white") +
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = c(.15,.95),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Perform Savitzky-Golay 2nd derivative on sample 258
sg.d <- as.data.table(gapDer(X = Points.d258[,2:1501], m = 2, w = 11, s = 5))

#Transpose the data for ggplot
sg.d.long <-  suppressWarnings(melt(setDT(sg.d), variable.name = "Wavelength", value.name = "Absorbance", variable.factor = FALSE))

#plot the transformed dark spectra
p.sg.d <-
  ggplot(sg.d.long) + 
  geom_line(aes(x = as.numeric(Wavelength), y = as.numeric(Absorbance), colour = ""), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  scale_color_manual(name = "Sample 258",
                     values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 10)) +
  theme_classic() +
  theme(legend.position = c(.1,.95),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(
    ggpubr::ggarrange(p.d, p.d258, labels = c("A", "B")),
    p.sg.d, labels = c("", "C"), nrow = 2)

#Save figure
ggsave("003-nir-savitzky-transform.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)