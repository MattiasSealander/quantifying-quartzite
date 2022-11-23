suppressPackageStartupMessages(library(ChemoSpec))
suppressPackageStartupMessages(library(spectrolab))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import raman data, set empty fields to NA
raman.csv <-
  read.csv2("./analysis/data/raw_data/RAMAN/raman_samples_data_non_treated_20220407.csv", sep = ";", dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
raman.averaged <- 
  aggregate(raman.csv[, 5:469], list(sample_id = raman.csv$sample_id), mean)

#merge raman de-trended data with metadata
raman.baseline.merged <- 
  as.data.frame(merge(metadata.csv, raman.averaged, by='sample_id'))

Points.raman <-
  raman.baseline.merged %>%
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
  replace_na(list(munsell_hue = "Colourless"))

#Filter raman data to focus on dark samples 
r.dark <- Points.raman %>%
  filter(hue == "Dark") %>% 
  select(1, c(30:494))

#Melt into long format
r.dark.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.dark), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data to focus on light samples
r.light <- Points.raman %>%
  filter(hue == "Light") %>% 
  select(1, c(30:494))

#Melt into long format
r.light.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.light), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data to focus on colourless samples
r.colourless <- Points.raman %>%
  filter(hue == "Colourless") %>% 
  select(1, c(30:494))

#Melt into long format
r.colourless.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.colourless), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#Filter raman data to focus on white samples 
r.white <- Points.raman %>%
  filter(hue == "White") %>% 
  select(1, c(30:494))

#Melt into long format
r.white.long <- 
  suppressWarnings(reshape2::melt(data.table::setDT(r.white), variable.name = "Wavenumber", variable.factor = FALSE, value.name = "Intensity"))

#plot the colourless spectra
p.c <-
  ggplot(r.colourless.long,aes(x = as.numeric(as.character(Wavenumber)))) + 
    geom_line(aes(y = Intensity, colour = "", group = sample_id), size = 1, stat = "identity") +
    xlab("Wavenumber (cm-1)") +
    ylab("Intensity") +
    scale_color_manual(name = "Colourless",
                       values = "#9467BDFF") + 
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the dark spectra
p.d <-
  ggplot(r.dark.long, aes(x = as.numeric(as.character(Wavenumber)))) + 
    geom_line(aes(y = Intensity, colour = "", group = sample_id), size = 1, stat = "identity") +
    xlab("Wavenumber (cm-1)") +
    ylab("Intensity") +
    scale_color_manual(name = "Dark",
                      values = "black") + 
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the light spectra
p.l <-
  ggplot(r.light.long, aes(x = as.numeric(as.character(Wavenumber)))) + 
    geom_line(aes(y = Intensity, colour = "", group = sample_id), size = 1, stat = "identity") +
    xlab("Wavenumber (cm-1)") +
    ylab("Intensity") +
    scale_color_manual(name = "Light",
                       values = "#FF7F0EFF") + 
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_text(size = 12, face = "bold", colour = "black"),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

#plot the white spectra
p.w <-
  ggplot(r.white.long, aes(x = as.numeric(as.character(Wavenumber)))) + 
    geom_line(aes(y = Intensity, colour = "", group = sample_id), size = 1, stat = "identity") +
    xlab("Wavenumber (cm-1)") +
    ylab("Intensity") +
    scale_color_manual(name = "White",
                       values = "#1F77B4FF") + 
    theme_classic() +
    theme(legend.position = c(.9,.95),
          axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
          axis.title.y = element_blank(),
          legend.title = element_text(size = 12, face = "bold", colour = "black"))

#Layout the plots in one figure
fig <-
  ggpubr::ggarrange(p.c, p.d, p.l, p.w, 
                    ncol = 2, 
                    nrow = 2)

#Save figure
ggsave("005-raman-spectra-summary.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)