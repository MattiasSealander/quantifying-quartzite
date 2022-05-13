#Load packages
suppressPackageStartupMessages(library(OmicsPLS))
suppressPackageStartupMessages(library(rgr))
suppressPackageStartupMessages(library(robCompositions))
suppressPackageStartupMessages(library(tidyverse))

#Import descriptive metadata
metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata_20220510.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#Import nir data, set empty fields to NA
xrf.csv <-
  read.csv2("./analysis/data/raw_data/XRF/xrf_quantitative_data_20220407.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Import nir data, set empty fields to NA
nir.csv <-
  read.csv2("./analysis/data/raw_data/NIR/asd_raw_data_20220127.csv", sep = ";", dec = ".", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#aggregate observations by group(sample) and calculate average of wavelength measurements
nir.averaged <- 
  aggregate(nir.csv[,4:2154], list(sample_id = nir.csv$sample_id), mean)

#merge NIR data with metadata
nir.merged <- 
  as.data.frame(merge(metadata.csv, nir.averaged, by='sample_id'))

#Fill the NA fields in the munsell_hue column and mark them as colourless samples
nir.merged[8][is.na(nir.merged[8])] <- "Colourless"

#merge NIR data with metadata
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
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite") %>% 
  filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                           "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                           "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) %>% 
  dplyr::select(sample_id, `Mg`, `Al`, `Si`, `K`, `Ca`, `Fe`, `Sr`, `Zr`, `Ti`)

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
  filter(!sample_id %in% c("153","167","168","169","172","174","175","177","182","183","190","191","193","194","196","198","200","204","207","210","213","214",
                           "215","216","229","234","235","237","238","251","262","265","268","269","272","278","281","282","359","377","385","392","393","397","405",
                           "406","410","411","413","414","415","416","417","424","425","426","428","430","432","55","56")) 

#impute missing data using least trimmed squares regression
Points.xrf[,2:10] <- 
  impAll(Points.xrf[,2:10])

#transform data using centered log-ratio 
Points.xrf[,2:10] <- 
  clr(Points.xrf[,2:10], ifclose = TRUE)

Points.xrf[,2:10] <- 
  scale(Points.xrf[,2:10])

Points.nir[,30:2180] <- 
  scale(Points.nir[,30:2180], scale = FALSE)

#merge NIR data with metadata
xrf.merged <- 
  as.data.frame(merge(Points.xrf, nir.averaged, by='sample_id'))

fit0 = o2m(Points.nir[,681:2180], Points.xrf[,2:10], 2, 3, 2)


p.x <- 
  #score plot of joint X scores 
  qplot(x=fit0$Tt[,1],
        y=fit0$Tt[,2],
        label = Points.nir$sample_id) +
  theme_bw() +
  labs(x = "Joint X - Component 1",
       y = "Joint X - Component 2") +
  geom_point(shape = 21, size = 4, aes(fill = Points.nir$hue)) +
  geom_text(hjust=0,
            vjust=-1,
            size=3) +
  scale_fill_manual(name = "Hue",
                    values = c("purple", "black","red","blue"))

#score plot of joint Y scores
p.y <-   
  qplot(x=fit0$U[,1], 
      y=fit0$U[,2],
      label = Points.nir$sample_id) +
  theme_bw() +
  labs(x = "Joint Y - Component 1",
       y = "Joint Y - Component 2") +
  geom_point(shape = 21, size = 4, aes(fill = Points.nir$hue)) +
  geom_text(hjust=0, 
            vjust=-1, 
            size=3) +
  scale_fill_manual(name = "Hue",
                    values = c("purple", "black","red","blue"))
fig <- 
  ggpubr::ggarrange(p.x, p.y,
                    nrow = 2,
                    common.legend = TRUE,
                    legend = "right")

#Save figure
ggsave("010-O2PLS.png",
       fig,
       device = "png",
       here::here("analysis/figures/"),
       width=20, 
       height=20,
       units = "cm",
       dpi = 300)
