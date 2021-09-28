library(ChemoSpec)
library(spectrolab)
library(tidyverse)

#Read and create new file for chemospec
nir <-matrix2SpectraObject(
  gr.crit = c("Dark", "Light", "Translucent", "White"),
  gr.cols = c("black", "blue", "red","purple"),
  freq.unit ="Wavelength (nm)",
  int.unit ="Absorbance",
  descrip ="Bifacial points measurements",
  in.file = "NIR_master_averaged_20210914_transposed.csv",
  sep = ";",
  dec = ",",
  chk = TRUE,
  out.file = "nir")

#calculate Savitsky-Golay 2nd derivative, m = derivative order, n = filter length
temp <- sgfSpectra(nir, m=2, n=9)

#Bind the data and save as a csv
tmp <- cbind(temp$freq, t(temp$data))
colnames(tmp) <- c("freq", temp$names)
tmp <- as.data.frame(tmp)
write.csv2(tmp, 'temp.csv')

#Import file after manually transposing it
nir.csv <-
  read.csv2("temp.csv", sep = ";",
            dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

# Set up a centering function with 'scale()'
center_scale <- function(x) {
  scale(x, scale = FALSE)
}

#center data
nir.csv[,3:1499] <- center_scale(nir.csv[,3:1499])

#write new file
write.csv2(nir.csv, 'temp.csv', row.names = FALSE)

#Import file again 
nir.csv <-
  read.csv2("temp.csv", sep = ";",
            dec = ",", header = TRUE, check.names = FALSE, na = c("","NA","NULL",NULL))

#Filter and set up spectra groups then plot result
dark <- nir.csv %>% 
  filter(hue == "Dark")

spec.d <- 
  as_spectra(dark, name_idx = 1, meta_idxs = 2)


light <- nir.csv %>% 
  filter(hue == "Light")

spec.l <- 
  as_spectra(light, name_idx = 1, meta_idxs = 2)

white <- nir.csv %>% 
  filter(hue == "White")

spec.w <- 
  as_spectra(white, name_idx = 1, meta_idxs = 2)

translucent <- nir.csv %>% 
  filter(hue == "Translucent")

spec.t <- 
  as_spectra(translucent, name_idx = 1, meta_idxs = 2)

par(mfrow=c(2,2))

plot(mean(spec.d), lwd = 1.2, col = "red", ylim=c(0.00005,-0.00005), cex.lab = 1.5, cex.main = 1.5, xaxp=c(1000,2500,10),
     xlab = "Wavelength (nm)", main = "Dark")

plot(mean(spec.l), lwd = 1.2, col = "red", ylim=c(0.00005,-0.00005), cex.lab = 1.5, cex.main = 1.5, xaxp=c(1000,2500,10),
     xlab = "Wavelength (nm)", main = "Light")

plot(mean(spec.w), lwd = 1.2, col = "red", ylim=c(0.00005,-0.00005), cex.lab = 1.5, cex.main = 1.5, xaxp=c(1000,2500,10),
     xlab = "Wavelength (nm)", main = "White")

plot(mean(spec.t), lwd = 1.2, col = "red", ylim=c(0.00005,-0.00005), cex.lab = 1.5, cex.main = 1.5, xaxp=c(1000,2500,10),
     xlab = "Wavelength (nm)", main = "Translucent")