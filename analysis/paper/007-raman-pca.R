#Load packages
suppressPackageStartupMessages(library(ChemoSpec))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(prospectr))
suppressPackageStartupMessages(library(tidyverse))

#Import raman data, set empty fields to NA
raman.csv <-
  read.csv2("./analysis/data/raw_data/RAMAN/raman_samples_data_non_treated_20220407.csv", sep = ";", dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

metadata.csv <-
  read.csv2("./analysis/data/raw_data/metadata_20220510.csv", sep = ";", header = TRUE, na = c("", "NA", "NULL"), encoding = "UTF-8")

#aggregate observations by group(sample) and calculate average of wavelength measurements
raman.averaged <- 
  aggregate(raman.csv[, 5:469], list(sample_id = raman.csv$sample_id), mean)

#transpose in prep for baseline de-trend
raman.transposed <- 
  t(raman.averaged)

#Write csv with transposed raman data for baseline de-trending in ChemoSpec package
#ChemoSpec requires data to be read from file
write.table(raman.transposed, file = "./analysis/data/derived_data/raman_transposed.csv", col.names = FALSE, sep = ";", dec = )

#prepare vector with sample_id for reading raman data into spectra object
sample_id <- 
  raman.averaged$sample_id

#read transposed raman data into spectra object
raman <-suppressWarnings(matrix2SpectraObject(
  gr.crit = sample_id,
  gr.cols = c("auto"),
  freq.unit ="Wavelength cm-1",
  int.unit ="Intensity",
  descrip ="Bifacial points measurements",
  in.file = "./analysis/data/derived_data/raman_transposed.csv",
  sep = ";",
  dec = ".",
  chk = TRUE,
  out.file = "./analysis/data/derived_data/raman_spec_object"))

#call on relevant baseline method (modified polynomial fitting) from baseline package and 
#return corrected spectra to spectra object
raman.baseline <- 
  baselineSpectra(raman,
                  int = FALSE,
                  method = "modpolyfit",
                  retC = TRUE)

#prevent baseline from plotting a figure
dev.off()

#spectra are stored in "Spectra" object, not a data frame
#need to extract samples and frequencies from that object 
raman.tmp <- 
  cbind(raman.baseline$freq, t(raman.baseline$data))
#transpose data to get frequencies as columns
raman.tmp <- 
  as.data.frame(t(raman.tmp))
#Add sample_id column and remove the initial "X" character introduced to sample names by R
raman.tmp$sample_id <- 
  c("sample_id", sub("^X*", "", raman.baseline$names))

#Make first row headers and remove the first row with the names
names(raman.tmp) <-
  raman.tmp[1,]
raman.tmp <- 
  raman.tmp[-1,]

#merge raman de-trended data with metadata
raman.baseline.merged <- 
  as.data.frame(merge(metadata.csv, raman.tmp, by='sample_id'))

#select and filter raman data to focus on relevant material and artefact types
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
         material == "Brecciated quartz" | material == "Quartz" | material == "Quartzite")

#fill the NA fields in the munsell hue column to mark them as colourless/translucent material
Points.raman[8][is.na(Points.raman[8])] <- "Colourless"

#perform PCA with SNV normalization and mean-center
raman.pca <-
  prcomp(standardNormalVariate(X = Points.raman[,c(30:321)]), center = TRUE, scale = FALSE)

#prepare labels for PCs
raman.pc1var <- round(summary(raman.pca)$importance[2,1]*100, digits=1)
raman.pc2var <- round(summary(raman.pca)$importance[2,2]*100, digits=1)
raman.pc3var <- round(summary(raman.pca)$importance[2,3]*100, digits=1)
raman.pc1lab <- paste0("PC1 (",as.character(raman.pc1var),"%)")
raman.pc2lab <- paste0("PC2 (",as.character(raman.pc2var),"%)")
raman.pc3lab <- paste0("PC3 (",as.character(raman.pc3var),"%)")

#extract components from pca and prepare for 3d visualization
components <- raman.pca[["x"]]

components <- data.frame(components)

components$PC2 <- -components$PC2

components$PC3 <- -components$PC3

components = cbind(components, Points.raman$hue)


#Plot using plotly with hue as color group
fig <- 
  plot_ly() %>%
  add_trace(data = components, 
            type = 'scatter3d', 
            mode = "markers",
            x = ~PC1, 
            y = ~PC3, 
            z = ~PC2,
            color = ~factor(Points.raman$hue),
            colors = "Set1",
            marker = list(opacity = 0.8,
                          line = list(
                            color = 'rgba(0,0,0)',
                            width = 0.2))) %>% 
  layout(
    showlegend=T,
    scene = list(
      xaxis=list(title=raman.pc1lab),
      yaxis=list(title=raman.pc3lab),
      zaxis=list(title=raman.pc2lab),
      camera = list(eye = list(x = -1.5, y = 1.0, z = 1.5))
    )
  )

save_image(fig,
           here::here("analysis/figures/007-raman-pca.png"),
           scale = 1, 
           width=1500, 
           height=1300,
           dpi = 300)