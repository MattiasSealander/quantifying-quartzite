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

#Import xrf data
XRF.csv <-
read.csv2("analysis/data/raw_data//XRF//XRF_quantitative.csv", sep = ";", dec = ".", header = TRUE)

#Import eigenvalues from Evince
Eig <-
    as.data.frame(read.csv2("analysis/data/raw_data//XRF_Eig_values.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Import loadings from Evince
P <-
  as.data.frame(read.csv2("analysis/data/raw_data//XRF_P_values.csv", sep = ";", dec = ".", header = TRUE, encoding = "UTF-8"))

#Import sample list, set empty fields to NA
samples <-
  as.data.frame(read.csv2("analysis/data/raw_data//bifacial_points.csv", sep = ";", dec = ".", header = TRUE, na.strings=c("", "NA")))

#Set reading no as row number
XRF <-
  csv2<- XRF.csv %>%
  remove_rownames %>%
  column_to_rownames(var="reading_no") %>%
  as.data.frame()

#Prepare filters for data
artefact <- c("Point", "Point fragment", "Preform")
material <- c("Quartzite", "Breccie quartz", "Breccie quartzite", "Felsitic porphyric quartz", "Quartz")

#Filter data to focus on quartzite points, point fragments and preforms made from quartz/quartzite materials
Points <-
  XRF.csv %>%
  filter(material == "Quartz" | material == "Quartzite" | material == "Breccie quartz" | material == "Breccie quartzite" | material  == "Felsitic porphyric quartz",
         type == "Point" | type == "Point fragment" | type == "Preform")

#Filter data to focus on quartzite points, point fragments and preforms made from quartz/quartzite materials
Points2 <-
  samples %>%
  filter(material == "Quartz" | material == "Quartzite" | material == "Breccie quartz" | material == "Breccie quartzite" | material  == "Felsitic porphyric quartz",
         type == "Point" | type == "Point fragment" | type == "Preform")


#Select only element columns
Points_xrf <-  Points %>%
  select(X12Mg:X92U)


colors <- c(na.value = "purple", "N1" = "grey4", "N2" = "grey7", "N3" = "grey9", "N3" = "grey11", "N4" = "grey13", "N5" = "grey15", "N6" = "grey17",
            "N8" = "grey38", "N9" = "gray68", "2.5PB" = "slateblue1", "2.5Y" = "cornsilk1","5Y" = "cornsilk2","7.5Y" = "cornsilk3",
            "10Y" = "bisque3","2.5GY" = "darkolivegreen2","5GY" = "darkolivegreen4","5YR" ="tomato1", "10B" = "skyblue4")


Points2 <-
  Points %>%
    select(X12Mg,X13Al, X14Si, X15P, X16S, X19K, X20Ca, X23V, X25Mn, X26Fe, X30Zn, X38Sr, X39Y, X40Zr, X56Ba, X81Ti)
var <- c("Mg", "Al", "Si", "P", "S", "K", "Ca", "V", "Mn", "Fe", "Zn", "Sr", "Y", "Zr", "Ba", "Ti")

#Boxplot of Mg content
Mg <- ggplot(Points,aes(y=X12Mg)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Mg") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X12Mg)$stats[c(1, 5)]
# scale y limits based on ylim1
Mg1 = Mg + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Al content
Al <- ggplot(Points,aes(y=X13Al)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Al") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X13Al)$stats[c(1, 5)]
# scale y limits based on ylim1
Al1 = Al + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Si content
Si <- ggplot(Points,aes(y=X14Si)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Si") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X14Si)$stats[c(1, 5)]
# scale y limits based on ylim1
Si1 = Si + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of P content
P <- ggplot(Points,aes(y=X15P)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("P") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X15P)$stats[c(1, 5)]
# scale y limits based on ylim1
P1 = P + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of S content
S <- ggplot(Points,aes(y=X16S)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("S") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X16S)$stats[c(1, 5)]
# scale y limits based on ylim1
S1 = S + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of K content
K <- ggplot(Points,aes(y=X19K)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("K") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X19K)$stats[c(1, 5)]
# scale y limits based on ylim1
K1 = K + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Ca content
Ca <- ggplot(Points,aes(y=X20Ca)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Ca") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X20Ca)$stats[c(1, 5)]
# scale y limits based on ylim1
Ca1 = Ca + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of V content
V <- ggplot(Points,aes(y=X23V)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("V") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X23V)$stats[c(1, 5)]
# scale y limits based on ylim1
V1 = V + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Mn content
Mn <- ggplot(Points,aes(y=X25Mn)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Mn") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X25Mn)$stats[c(1, 5)]
# scale y limits based on ylim1
Mn1 = Mn + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Fe content
Fe <- ggplot(Points,aes(y=X26Fe)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Fe") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X26Fe)$stats[c(1, 5)]
# scale y limits based on ylim1
Fe1 = Fe + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Zn content
Zn <- ggplot(Points,aes(y=X30Zn)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Zn") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X30Zn)$stats[c(1, 5)]
# scale y limits based on ylim1
Zn1 = Zn + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Sr content
Sr <- ggplot(Points,aes(y=X38Sr)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Sr") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X38Sr)$stats[c(1, 5)]
# scale y limits based on ylim1
Sr1 = Sr + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Y content
Y <- ggplot(Points,aes(y=X39Y)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Y") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X39Y)$stats[c(1, 5)]
# scale y limits based on ylim1
Y1 = Y + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Zr content
Zr <- ggplot(Points,aes(y=X40Zr)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Zr") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X40Zr)$stats[c(1, 5)]
# scale y limits based on ylim1
Zr1 = Zr + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Ba content
Ba <- ggplot(Points,aes(y=X56Ba)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Ba") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X56Ba)$stats[c(1, 5)]
# scale y limits based on ylim1
Ba1 = Ba + coord_cartesian(ylim = ylim1*1.05)

#Boxplot of Ti content
Ti <- ggplot(Points,aes(y=X81Ti)) +
  geom_boxplot(width=0.5) +
  theme_bw() +
  scale_y_continuous() +
  scale_x_discrete() +
  ylab("Ti") +
  theme(
    axis.title.y = element_text(size = 10, face = "bold", colour = "black"))
# compute lower and upper whiskers
ylim1 = boxplot.stats(Points$X81Ti)$stats[c(1, 5)]
# scale y limits based on ylim1
Ti1 = Ti + coord_cartesian(ylim = ylim1*1.05)

# 6 figures arranged in 3 rows and 4 columns
grid.arrange(Mg1 ,Al1 ,Si1 ,P1 ,S1 ,K1 ,Ca1 ,V1 ,Mn1 ,Fe1 ,Zn1 ,Sr1 ,Y1 ,Zr1 ,Ba1 ,Ti1 , ncol=4)

#Si content bar graph
Points %>%
  drop_na(X14Si) %>%
  ggplot(aes(x= reorder(reading_no, X14Si), y = X14Si)) +
  geom_bar(stat = "identity", position = "identity", fill = "black") +
  theme_minimal() +
  xlab("Readings") +
  ylab("SiO -%") +
  theme(
    axis.text.x=element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#calculate mean of a column
print(summary(Points$Bal))

#no missing rows MgO
dplyr::count(Points, X12Mg)
#no missing rows Al2O3
dplyr::count(Points, X13Al)
#no missing rows SiO2
dplyr::count(Points, X14Si)
#no missing rows K20
dplyr::count(Points, X19K)
#no missing rows MnO
dplyr::count(Points, X25Mn)
#no missing rows Fe2O3
dplyr::count(Points, X26Fe)

#Barplot of MgO calculated estimates for individual readings, 287 rows NA
b1 <- Points %>%
  drop_na(X12Mg) %>%
  ggplot(aes(x= reorder(reading_no, X12Mg), y = X12Mg)) +
  geom_bar(stat = "identity", position = "identity", fill = "black") +
  theme_minimal() +
  ggtitle("Mg") +
  xlab("Readings") +
  ylab("%") +
  theme(
    plot.title = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 20, b = -20)),
    axis.text.x=element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#Barplot of Al2O3 calculated estimates for individual readings, 6 NA
b2 <- Points %>%
  drop_na(X13Al) %>%
    ggplot(aes(x= reorder(reading_no, X13Al), y = X13Al)) +
    geom_bar(stat = "identity", position = "identity", fill = "black") +
    theme_minimal() +
    ggtitle("Al") +
    xlab("Readings > 0: 466") +
  ylab("%") +
    theme(
      plot.title = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 20, b = -20)),
      axis.text.x=element_blank(),
      axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#Barplot of SiO2 calculated estimates for individual readings
b3 <- Points %>%
  drop_na(X14Si) %>%
    ggplot(aes(x= reorder(reading_no, X14Si), y = X14Si)) +
    geom_bar(stat = "identity", position = "identity", fill = "black") +
    theme_minimal() +
    ggtitle("Si") +
    xlab("Readings") +
    ylab("%") +
    theme(
      plot.title = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 20, b = -20)),
      axis.text.x=element_blank(),
      axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
      axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#Barplot of K2O calculated estimates for individual readings
b4 <- Points %>%
  drop_na(X19K) %>%
  ggplot(aes(x= reorder(reading_no, X19K), y = X19K)) +
  geom_bar(stat = "identity", position = "identity", fill = "black") +
  theme_minimal() +
  ggtitle("K") +
  xlab("Readings") +
  ylab("%") +
  theme(
    plot.title = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 20, b = -20)),
    axis.text.x=element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#Barplot of MnO calculated estimates for individual readings, 318 NA
b5 <- Points %>%
  drop_na(X25Mn) %>%
  ggplot(aes(x= reorder(reading_no, X25Mn), y = X25Mn)) +
  geom_bar(stat = "identity", position = "identity", fill = "black") +
  theme_minimal() +
  ggtitle("Mn") +
  xlab("Readings") +
  ylab("%") +
  theme(
    plot.title = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 20, b = -20)),
    axis.text.x=element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"))

#Barplot of Fe2O3 calculated estimates for individual readings, 29 NA
b6 <- Points %>%
  drop_na(X26Fe) %>%
  ggplot(aes(x= reorder(reading_no, X26Fe), y = X26Fe)) +
  geom_bar(stat = "identity", position = "identity", fill = "black") +
  theme_minimal() +
  ggtitle("Fe") +
  xlab("Readings") +
  ylab("%") +
  theme(
    plot.title = element_text(size = 14, face = "bold", colour = "black", margin = margin(t = 20, b = -20)),
    axis.text.x=element_blank(),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold", colour = "black"))


# 6 figures arranged in 3 rows and 2 columns
grid.arrange(b1, b2, b3, b4, b5, b6, ncol=2)



#stacked bar chart with counts of samples per parish and type
#if geom_bar() is left empty it defaults to counts
ggplot(Points2, aes(x= fct_rev(parish), fill = type)) +
  geom_bar() +
  theme_minimal() +
  xlab("Count") +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(face = "bold"),
    axis.title.y=element_blank(),
    axis.text.y = element_text(size = 10, face = "bold", colour = "black"),
    axis.title.x = element_text(size = 12, face = "bold", colour = "black")) +
  scale_fill_jco() +
  coord_flip() +
  geom_text(stat ='count', aes(group = parish, label = ..count..), hjust = -0.5)

#Drop NA rows for site_id using tidyverse, make stacked bar chart with counts per site
#if geom_bar() is left empty it defaults to counts
Points2 %>%
  drop_na(site_id) %>%
      group_by(site_id) %>%
      filter(n() > 5) %>%
        ggplot(aes(x= site_id, fill = type)) +
          geom_bar() +
          theme_minimal() +
          xlab("Count") +
          theme(
            legend.position = "bottom",
            legend.title = element_text(size = 16, face = "bold"),
            legend.text = element_text(face = "bold"),
            axis.title.y=element_blank(),
            axis.text.x = element_text(angle = 90, size = 12, face = "bold", colour = "black"),
            axis.title.x = element_text(size = 12, face = "bold", colour = "black")) +
          scale_fill_jco() +
          geom_text(stat ='count', aes(group = site_id, label = ..count..), vjust = -0.5)


#Barplot of Fe2O3 calculated estimates for individual readings
#Not sure how null munsell_hue values are displayed, the legend seems to be the color of background
filter(Points, !is.na(X26Fe)) %>%
  ggplot() +
    geom_bar(mapping = aes(x= reorder(reading_no, X26Fe), y = X26Fe,  fill = munsell_hue), stat = "identity", position = "identity") +
    scale_fill_manual(values = colors, name = "Munsell") +
    theme_minimal() +
    theme(axis.text.x=element_blank()) +
    xlab("Readings") +
    ylab(bquote(~Fe[2]~O[3]~"-%"))


#Barplot of Al2O3 calculated estimates for individual readings
#Not sure how null munsell_hue values are displayed, the legend seems to be the color of background
filter(Points, !is.na(X13Al)) %>%
  ggplot() +
  geom_bar(mapping = aes(x= reorder(reading_no, X13Al), y = X13Al,  fill = munsell_hue), stat = "identity", position = "identity") +
  scale_fill_manual(values = colors, name = "Munsell") +
  theme_minimal() +
  theme(axis.text.x=element_blank()) +
  xlab("Readings") +
  ylab(bquote(~Al[2]~O[3]~"-%"))


#Scatterplot Al - Si, 10 missing values
ggplot(Points, aes(x=X14Si, y=X13Al, color=munsell_hue)) +
  geom_point() +
  scale_color_manual(values = colors, name = "Munsell") +
  theme_dark() +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~Al[2]~O[3]~"-%"))


#Scatterplot Fe - Si, 31 missing values
ggplot(Points, aes(x=X14Si, y=X26Fe, color=munsell_hue)) +
  geom_point() +
  scale_color_manual(values = colors, name = "Munsell") +
  theme_dark() +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~Fe[2]~O[3]~"-%"))

#Scatterplot Mg - Si, 287 missing values
ggplot(Points, aes(x=X14Si, y=X12Mg, color=munsell_hue)) +
  geom_point() +
  scale_color_manual(values = colors, name = "Munsell") +
  theme_dark() +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~Fe[2]~O[3]~"-%"))

#Scatterplot Bal - Si, 36 missing values
ggplot(Points, aes(x=X14Si, y=Bal, color=munsell_hue)) +
  geom_point() +
  scale_color_manual(values = colors, name = "Munsell") +
  theme_dark() +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~"Bal-%"))



##Prepare scatterplots for 4-plot figure
#Scatterplot Al - Si, 6 missing values
Al <- ggplot(Points, aes(x=X14Si, y=X13Al, fill="black")) +
  geom_point() +
  theme_minimal() +
  ggtitle("B") +
  theme(legend.position = "none") +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~Al[2]~O[3]~"-%"))

#Scatterplot Bal - Si, 36 missing values
Bal <- ggplot(Points, aes(x=X14Si, y=Bal, fill="black")) +
  geom_point() +
  theme_minimal() +
  ggtitle("A") +
  theme(legend.position = "none") +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~"Bal-%"))

#Scatterplot Fe - Si, 31 missing values
Fe <- ggplot(Points, aes(x=X14Si, y=X26Fe, fill="black")) +
  geom_point() +
  theme_minimal() +
  ggtitle("C") +
  theme(legend.position = "none") +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab(bquote(~Fe[2]~O[3]~"-%"))

#Scatterplot Mg - Si, 287 missing values
Mg <- ggplot(Points, aes(x=X14Si, y=X12Mg, fill="black")) +
  geom_point() +
  theme_minimal() +
  ggtitle("D") +
  theme(legend.position = "none") +
  xlab(bquote(~SiO[3]~"-%")) +
  ylab("MgO-%")


# 4 figures arranged in 2 rows and 2 columns
grid.arrange(Bal, Al, Fe, Mg, ncol=2)


##Prepare plots for 3-plot bargraph figure of loadings
#PC1 loadings
P1 <-
  ggplot(P, aes(x = fct_rev(Element), y = P1)) +
    geom_bar(
      stat = "identity", position = "identity",
      color = "black", fill = "darkgrey") +
  theme_minimal() +
  theme(plot.margin = margin(10,10,10,10),
        axis.title.y=element_blank(),
        axis.text.y = element_text(size=14, face="bold", colour = "black"),
        axis.title.x = element_text(size = 14, face = "bold", colour = "black")) +
  coord_flip()
#PC2 loadings
P2 <-
  ggplot(P, aes(x = fct_rev(Element), y = P2)) +
  geom_bar(
    stat = "identity", position = "identity",
    color = "black", fill = "darkgrey") +
  theme_minimal() +
  theme(plot.margin = margin(10,10,10,10),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", colour = "black")) +
  coord_flip()
#PC3 loadings
P3 <-
  ggplot(P, aes(x = fct_rev(Element), y = P3)) +
  geom_bar(
    stat = "identity", position = "identity",
    color = "black", fill = "darkgrey") +
  theme_minimal() +
  theme(plot.margin = margin(10,10,10,10),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 14, face = "bold", colour = "black")) +
  coord_flip()

#Composite figure of above bar graphs
grid.arrange(P1, P2, P3, ncol=3)






#Remove columns with only na
Points_na <- Filter(function(x)!all(is.na(x)), Points_xrf)

#Drop columns containing between 1-5 values as they may interfer with imputation
Points_na = subset(Points_na, select = -c(X92U))

Points_na = subset(Points_na, select = c(X12Mg, X13Al, X14Si, X15P, X16S,	X17Cl,	X19K,	X20Ca,	X23V, X25Mn,	X26Fe,	X27Co, X30Zn, X37Rb,	X38Sr,	X39Y,	X40Zr, X50Sn, X56Ba,	X81Ti, X90Th
))

#Impute missing values using PCA
imp <- imputePCA(Points_na,
                 ncp = 2,
                 scale = FALSE,
                 method = "Regularized",
                 maxiter = 1)

test <- pcaCoDa(imp$completeObs, method = "classical", solve = "eigen")

p1 <- pcaCoDa(test)


biplot(test, col = c(Points$hue, "red"))





