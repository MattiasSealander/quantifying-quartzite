# Loading package
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(flextable))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2("./analysis/data/raw_data/NIR/asd_raw_data_20220407.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Import xrf data, set empty fields to NA
xrf.csv <-
  read.csv2("./analysis/data/raw_data/XRF/xrf_quantitative_data_20220407.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL"))

#aggregate observations by group(sample) and calculate average of wavelength measurements
nir.averaged <- 
  aggregate(nir.csv[, 4:2154], list(sample_id = nir.csv$sample_id), mean)

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.averaged, by='sample_id'))

#merge XRF data with metadata
xrf.merged <- 
  as.data.frame(merge(metadata.csv, xrf.csv, by='sample_id'))

#Filter nir data to focus on points and preforms made from quartz/quartzite material 
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
                material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") 

#Filter xrf data to focus on points and preforms made from quartz/quartzite material 
Points.xrf <-
  xrf.merged %>%
  dplyr::select(sample_id, `Al`, `Si`, `K`, `Ca`, `Fe`, `Zr`, `Ti`) %>% 
  dplyr::filter(sample_id == 231 | sample_id == 249| sample_id == 391)

ft <- 
  flextable(Points.xrf) %>% 
  set_header_labels(sample_id = "Sample") %>% 
  theme_vanilla() %>% 
  set_caption(caption = "XRF - Elemental content (%)") %>% 
  as_raster()

x.ft <- 
  ggplot() + 
  theme_void() + 
  annotation_custom(rasterGrob(ft), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
n.231 <- Points.nir %>%
  dplyr::filter(sample_id == 231) %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
n.231.long <- 
  suppressWarnings(melt(setDT(n.231), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
n.249 <- Points.nir %>%
  dplyr::filter(sample_id == 249) %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
n.249.long <- 
  suppressWarnings(melt(setDT(n.249), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#Filter NIR data to focus on colourless 
#and select the NIR range 1 000 - 2 500 nm
n.391 <- Points.nir %>%
  dplyr::filter(sample_id == 391) %>% 
  dplyr::select(1, c(681:2180))

#Melt into long format
n.391.long <- 
  suppressWarnings(melt(setDT(n.391), variable.name = "Wavelength", variable.factor = FALSE, value.name = "Absorbance"))

#plot the colourless spectra
p.231 <-
  ggplot(n.231.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = ""), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  labs(title="231") +
  scale_color_manual(name = "231",
                     values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = 0.02),        
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the colourless spectra
p.249 <-
  ggplot(n.249.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = ""), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  labs(title="249") +
  scale_color_manual(name = "249",
                     values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = 0.02),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))


#plot the colourless spectra
p.391 <-
  ggplot(n.391.long, aes(x = as.numeric(Wavelength))) + 
  geom_line(aes(y = Absorbance, colour = ""), size = 1, stat = "identity") +
  xlab("Wavelength (nm)") +
  ylab("Absorbance") +
  labs(title="391") +
  scale_color_manual(name = "391",
                     values = "black") + 
  scale_x_continuous(limits = c(1000, 2500), breaks = scales::pretty_breaks(n = 5)) +
  theme_classic() +
  theme(legend.position = "none",
        plot.title = element_text(size=16, 
                                  hjust = 0.02),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
        legend.title = element_text(size = 12, face = "bold", colour = "black"))



top <- 
  plot_grid(x.ft)
bottom <- 
  plot_grid(p.231, p.249, p.391,nrow=1)

fig <- 
  cowplot::plot_grid(top,bottom, nrow = 2)

#Save figure
ggsave("011-nir-sample-comparison.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       width=25, 
       height=10,
       units = "cm",
       dpi = 300)
