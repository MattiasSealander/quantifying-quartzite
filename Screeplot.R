library(robCompositions)
library(tidyverse)
library(factoextra)
library(FactoMineR)
library(missMDA)
library(VIM)
library(ggbiplot)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(scales)
library(plyr)


#Import eigenvalues from XRF Evince
XEig <-
  as.data.frame(read.csv2("analysis/data/raw_data//XRF_Eig_values.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Import eigenvalues from NIR Evince
NEig <-
  as.data.frame(read.csv2("analysis/data/raw_data//XRF_Eig_values.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Import eigenvalues from visNIR Evince
vNEig <-
  as.data.frame(read.csv2("analysis/data/raw_data//XRF_Eig_values.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Import eigenvalues from Raman Evince
REig <-
  as.data.frame(read.csv2("analysis/data/raw_data//XRF_Eig_values.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#XRF Screeplot of eigenvalues
ggplot(Eig, aes(x=factor(Vector), y=Eigenvalue, group=1)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10),"Eigenvalue") +
  theme_minimal() +
  xlab("Components")

#NIRScreeplot of eigenvalues
ggplot(Eig, aes(x=factor(Vector), y=Eigenvalue, group=1)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10),"Eigenvalue") +
  theme_minimal() +
  xlab("Components")

#visNIR Screeplot of eigenvalues
ggplot(Eig, aes(x=factor(Vector), y=Eigenvalue, group=1)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10),"Eigenvalue") +
  theme_minimal() +
  xlab("Components")

#Raman Screeplot of eigenvalues
ggplot(Eig, aes(x=factor(Vector), y=Eigenvalue, group=1)) +
  geom_line() +
  geom_point() +
  scale_y_continuous(limits = c(0,10), breaks = seq(0,10),"Eigenvalue") +
  theme_minimal() +
  xlab("Components")
