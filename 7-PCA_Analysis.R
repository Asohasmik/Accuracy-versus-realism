# Title:    PCA analysis for herpetufauna of Romania and background points
# Author:   Laura Karolin Steib and Dennis Rödder
# Date:     17.09.2025

library(sp)
library (sf)
library(rSDM)
library(terra)
library(raster)
library(stringr)

################################################################################
# PCA Analysis
################################################################################

# merge the files with the records

wd <- "D:/.../Records"
setwd(wd)

occ <- read.csv("D:/.../species.training.csv")
sp <- names(table(occ$Species))

sp.1 <- read.csv("records_Ablepharus_kitaibelii.csv")
rec.1 <- as.data.frame(nrow(sp.1))
colnames (rec.1)[1] <- "Records"
rec.1[,2] <- sp[1]
colnames(rec.1)[2] <- "Species"

for (j in 2:length(sp)){
  if (j == 37) next 
  sp.j <- read.csv(paste(wd, "/records_", sp[j], ".csv", sep=""))
  rec.j <- as.data.frame(nrow(sp.j))
  colnames (rec.j)[1] <- "Records"
  rec.j[,2] <- sp[j]
  colnames(rec.j)[2] <- "Species"
  rec.1 <- rbind(rec.1, rec.j)
  write.csv(rec.1, "NumberOfRecords.csv", row.names=F)
  sp.1 <- rbind(sp.1, sp.j)
  write.csv(sp.1, "AllRecords.csv", row.names=F)
}

# file with number of records, slope, normalized Average 
rec <- read.csv("D:/.../Records/NumberOfRecords.csv")
lm_ges <- read.csv("D:/.../MaxEntAllYears/_ClimaticallySuitCond/RasterWithBuffer/AreaSize/EachYear/lm_results_02_23.csv")
lm_ges <- as.data.frame(lm_ges[,2:3])
colnames(lm_ges)[2] <- "SlopeGesamt"
lm_ges$Species <- gsub(" ", "_", lm_ges$Species)
lm_ges$Species <- gsub("\\*", "", lm_ges$Species)
lm_rm <- read.csv("D:/.../MaxEntAllYears/_ClimSuit_Mask500Comb/RasterWithBuffer/AreaSize/EachYear/lm_results_02_23.csv")
lm_rm <- as.data.frame(lm_rm[,2:3])
colnames(lm_rm)[2] <- "SlopeMaskrm"
lm_rm$Species <- gsub(" ", "_", lm_rm$Species)
lm_rm$Species <- gsub("\\*", "", lm_rm$Species)
lm_whole <- read.csv("D:/.../MaxEntAllYears/_ClimSuitMaskWhole_Comb/RasterWithBuffer/AreaSize/EachYear/lm_results_02_23.csv")
lm_whole <- as.data.frame(lm_whole[,2:3])
colnames(lm_whole)[2] <- "SlopeMaskwhole"
lm_whole$Species <- gsub(" ", "_", lm_whole$Species)
lm_whole$Species <- gsub("\\*", "", lm_whole$Species)

comb <- merge(rec,lm_ges, by= "Species")
comb <- merge(comb, lm_rm, by = "Species")
comb <- merge(comb, lm_whole, by = "Species")
write.csv(comb, "D:/.../PCA/forGLM.csv")

hab_ges <- read.csv("D:/.../MaxEntAllYears/_ClimSuitMaskWhole_Comb/RasterWithBuffer/AreaSize/EachYear/normHabSizewithAverage.csv", row.names = "Years")
hab_ges <- t(hab_ges)
av <- hab_ges[,"Average "]
av <- as.data.frame(av)
colnames(av)[1] <- "AverageMaskwhole"
sp <- rownames(av)
av[2] <- c(sp)
colnames (av)[2] <- "Species"

comb <- merge (comb, av, by ="Species")

# extract PCA for coordinates 
PC1 <- rast("D:/.../PCA/PC1.asc")
PC2 <- rast("D:/.../PCA/PC2.asc")
records <- read.csv("D:/.../PCA/AllRecords.csv")

ext <- extract(PC2, records[,2:3], ID=FALSE)

add <- cbind(add,ext)

write.csv(add, "D:/.../PCA/AllRecordsWithPCA.csv")

# get the boxplot stats for PC1/PC2
wd <- "D:/.../PCA"
setwd(wd)

data <- read.csv( "D:/.../PCA/AllRecordsWithPCA.csv")
data[,1] <- rm()
data <- data [, c(1:3,13:14)]
sp <- names(table(data$Species))

out2 <- c()
for (j in 1:length(sp)){
  rec.j <- subset(data, data$Species == sp[j])
  rec.j <- as.data.frame(rec.j [, "PC2"])
  colnames (rec.j) [1] <- "PC2"
  out <- c()
  Species <- sp[j]
  out <- cbind(out,Species)
  med_PC2 <- median(rec.j$PC2)
  mean_PC2 <- mean(rec.j$PC2)
  out <- cbind(out,mean_PC2, med_PC2)
  sd_PC2 <- sd(rec.j$PC2)
  out <- cbind(out,sd_PC2)
  t_value <- qt(0.95, nrow(rec.j)-1)              # Nutzung der t-Verteilung 
  SE <- sd_PC2/sqrt(nrow(rec.j))
  CI_lower_PC2 <- mean (rec.j$PC2) - t_value *SE
  CI_upper_PC2 <- mean (rec.j$PC2) + t_value *SE
  out <- cbind(out, CI_lower_PC2,CI_upper_PC2)
  Q1 <- quantile(rec.j$PC2, 0.25)
  Q3 <- quantile(rec.j$PC2, 0.75)
  IQR <- Q3-Q1
  multiplier <- qnorm(0.95)
  whisker_low_PC2 <- max(min(rec.j$PC2), Q1 - multiplier * IQR)
  whisker_high_PC2 <- min(max(rec.j$PC2), Q3 + multiplier * IQR)
  out <- cbind (out, whisker_low_PC2, whisker_high_PC2)
  out <- as.data.frame(out)
  out2 <- rbind(out2,out)
}
write.csv(out2, "Stats.PC2.csv")

PC1 <- read.csv("Stats.PC1.csv")
PC1 [,1] <- rm()
PC2 <- read.csv("Stats.PC2.csv")
PC2[,1] <- rm()

comb <- merge(PC1, PC2, by = "Species")
write.csv(comb, "PCAStats.csv")

# get rid of 10% of the data and get the the min and max value than 
out <- c()
rec.j <- subset(data, data$Species == sp[1])
PC1_lower_q <- quantile(rec.j$PC1, 0.05)
PC1_upper_q <- quantile(rec.j$PC1, 0.95)
PC2_lower_q <- quantile(rec.j$PC2, 0.05)
PC2_upper_q <- quantile(rec.j$PC2, 0.95)
out <- as.data.frame(cbind(PC1_lower_q,PC1_upper_q,PC2_lower_q, PC2_upper_q))
rownames(out) <- sp[1]
for (j in 2:length(sp)){
  rec.j <- subset(data, data$Species == sp[j])
  PC1_lower_q <- quantile(rec.j$PC1, 0.05)
  PC1_upper_q <- quantile(rec.j$PC1, 0.95)
  PC2_lower_q <- quantile(rec.j$PC2, 0.05)
  PC2_upper_q <- quantile(rec.j$PC2, 0.95)
  out2 <- as.data.frame(cbind(PC1_lower_q,PC1_upper_q,PC2_lower_q, PC2_upper_q))
  rownames(out2) <- sp[j]
  out <- rbind(out,out2)
}

# merge with other species data
GLM <- read.csv("forGLM.csv")
GLM[,1] <- rm()
PCA <- read.csv("PCAStats.csv")
PCA[,1] <- rm()
comb <- merge(GLM,PCA, by = "Species")
write.csv(comb, "GLM_PCA.csv", row.names = FALSE)

# plot PCA
library(ggplot2)
library(ggrepel)
data <- read.csv("GLM_PCA.csv")			  

SlopeGesamt <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = SlopeGesamt), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = SlopeGesamt), size = 0.01) +  # Vertical whisker for PC2
  geom_point(aes(color = SlopeGesamt, shape = Category, size = Records),show.legend=FALSE) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "blue", high = "red") +
  #geom_text_repel(data = data, 
  #              aes(label = Species), 
  #             box.padding = 0.5, 
  #            fontface = c('italic'),
  #           size = 3, 
  #          max.overlaps = Inf,  # Ensures all labels are drawn
  #   		segment.color = 'black',  # Color of the connecting line for geom_text_repel
  #       segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Slope of Change in Habitat Size without Habitat Mask")	

plot(SlopeGesamt)

SizeGesamt <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = AverageSizeTotal), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = AverageSizeTotal), size = 0.01) +  # Vertical whisker for PC2
  geom_point(aes(color = AverageSizeTotal, shape = Category, size = Records),show.legend=FALSE) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "#EE7942", high = "#008B45",limits = c(0,1)) +
  geom_text_repel(data = data, 
                  aes(label = Species), 
                  box.padding = 0.5, 
                  fontface = c('italic'),
                  size = 3, 
                  max.overlaps = Inf,  # Ensures all labels are drawn
                  segment.color = 'black',  # Color of the connecting line for geom_text_repel
                  segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Standardized Average Habitat Size without Habitat Mask")	


SlopeRemoved <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = SlopeMaskrm), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = SlopeMaskrm), size = 0.01) +  # Vertical whisker for PC2
  geom_point(aes(color = SlopeMaskrm , shape = Category, size = Records),show.legend=FALSE) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "blue", high = "red") +
  geom_text_repel(data = data, 
                  aes(label = Species), 
                  box.padding = 0.5, 
                  fontface = c('italic'),
                  size = 3, 
                  max.overlaps = Inf,  # Ensures all labels are drawn
                  segment.color = 'black',  # Color of the connecting line for geom_text_repel
                  segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Slope of Change in Habitat Size with Habitat Mask (Marginal Habitat Types Removed)")	


SizeRemoved <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = AverageSizeMaskrm), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = AverageSizeMaskrm), size = 0.01) +  # Vertical whisker for PC2
  geom_point(aes(color = AverageSizeMaskrm, shape = Category, size = Records),show.legend=FALSE) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "#EE7942", high = "#008B45",limits = c(0,1)) +
  geom_text_repel(data = data, 
                  aes(label = Species), 
                  box.padding = 0.5, 
                  fontface = c('italic'),
                  size = 3, 
                  max.overlaps = Inf,  # Ensures all labels are drawn
                  segment.color = 'black',  # Color of the connecting line for geom_text_repel
                  segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Standardized Average Habitat Size with Habitat Mask(Marginal Habitat Types Removed)")	


SlopeWhole <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = SlopeMaskwhole), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = SlopeMaskwhole), size = 0.01) +  # Vertical whisker for PC2
  geom_point(aes(color = SlopeMaskwhole , shape = Category, size = Records),show.legend=FALSE) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "blue", high = "red") +
  geom_text_repel(data = data, 
                  aes(label = Species), 
                  box.padding = 0.5, 
                  fontface = c('italic'),
                  size = 3, 
                  max.overlaps = Inf,  # Ensures all labels are drawn
                  segment.color = 'black',  # Color of the connecting line for geom_text_repel
                  segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Slope of Change in Habitat Size with Habitat Mask (All Habitat Types)")	


SizeWhole <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = AverageSizeMaskwhole), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = AverageSizeMaskwhole), size = 0.01) +  # Vertical whisker for PC2
  geom_point(aes(color = AverageSizeMaskwhole, shape = Category, size = Records),show.legend=T) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "#EE7942", high = "#008B45",limits = c(0,1)) +
  geom_text_repel(data = data, 
                  aes(label = Species), 
                  box.padding = 0.5, 
                  fontface = c('italic'),
                  size = 3, 
                  max.overlaps = Inf,  # Ensures all labels are drawn
                  segment.color = 'black',  # Color of the connecting line for geom_text_repel
                  segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Standardized Average Habitat Size with Habitat Mask(All Habitat Types)")	
plot(SizeWhole)	 

par(mfrow=c(2,3))  
grid_arrange <- grid.arrange(
  SlopeGesamt, SizeGesamt, SlopeRemoved,
  SizeRemoved, SlopeWhole, SizeWhole,
  ncol = 2  # Number of columns in the grid
)
pdf("PDA_withNames.pdf")
plot(SlopeGesamt)
plot(SizeGesamt)
plot(SlopeRemoved)
plot(SizeRemoved)
plot(SlopeWhole)
plot(SizeWhole)
dev.off()

# GLM
library(MASS)

GLM <- read.csv("GLM_PCA.csv")

glm_size_ges <- glm(formula =  AverageSizeTotal ~ SlopeGesamt , data =GLM, family = gaussian()) # AIC -16.81
glm_size_rm <- glm(formula =  AverageSizeMaskrm ~ med_PC2   , data =GLM, family = gaussian()) #AIC - 10.67
glm_size_whole <- glm(formula =  AverageSizeMaskwhole ~  med_PC2  , data =GLM, family = gaussian()) # AIC - 16.47

glm_slope_ges <- glm(formula =  SlopeGesamt ~ med_PC1 + med_PC2  , data =GLM, family = gaussian())
glm_slope_rm <- glm(formula =  SlopeMaskrm ~  med_PC1 + med_PC2 , data =GLM, family = gaussian())
glm_slope_whole <- glm(formula =  SlopeWhole ~  med_PC1  + med_PC2   , data =GLM, family = gaussian())

summary(glm_slope_ges)
summary(glm_slope_rm)
summary(glm_slope_whole)

dropterm(glm_size_ges, test = "Chisq")

stepAIC(glm_slope_ges)
stepAIC(glm_slope_rm)
stepAIC(glm_slope_whole)

sink("glm_slope_whole.txt")
print(summary(glm_slope_whole))
sink()

################################################################################
# # PCA Density Plotimag
################################################################################

wd <- "D:/.../PCA"kde <- kde2d(values12[,1], values12[,2], n = 1000)
setwd(wd)

library(terra)
library(MASS)
library(raster)
library(fields)

pc1 <- rast("PC1.asc")
pc2 <- rast("PC2.asc")

values1 <- values(pc1)
values2 <- values(pc2)
values12 <- cbind(values1, values2)
values12 <- na.omit(values12)
values12 <- values12[sample(nrow(values12),replace=F,size=100000),]


kde <- kde2d(values12[,1], values12[,2], n = 1000)
image(kde, col=tim.colors(100), ylab = "PC2 (17.33%)", xlab= "PC1 (71.0%)")
rast.kde <- raster(kde)
kde.values <- extract(rast.kde, values12[,1:2])

rast.1 <- rast.kde
rast.1[rast.1 < min(kde.values)] <- NA
plot(rast.1, col=tim.colors(100))

writeRaster(rast.1, "PCA.density.tif")
image(kde)


save(kde, file = "kde_output.RData")
load("kde_output.RData")

kde_df <- data.frame(expand.grid(x = kde$x, y = kde$y), density = as.vector(kde$z))
kde_df$density[kde_df$density < min(kde.values)] <- NA

tim_colors <- tim.colors(100)
tim_colors_scale <- scale_fill_gradientn(colours = alpha(c("blue", "green", "yellow", "red"), 0.6),na.value = "transparent")

# Visualize KDE with ggplot2
p <- ggplot(kde_df, aes(x = x, y = y, fill = density)) +
  geom_raster(na.rm = TRUE, interpolate = TRUE) +
  tim_colors_scale +
  theme_minimal() +
  labs(title = "KDE Visualization with tim.colors",
       x = "X-axis",
       y = "Y-axis",
       fill = "Density")
#########	   
combined_plot <- ggplot() +
  # KDE raster layer
  geom_raster(data = kde_df, aes(x = x, y = y, fill = density), interpolate = TRUE) +
  tim_colors_scale +  # Adjust color scale as needed
  
  # Data points layer
  geom_point(data = data, aes(x = med_PC1, y = med_PC2, color = SlopeGesamt, shape = Category, size = Records), show.legend = FALSE) +  # Data points
  scale_color_gradient(low = "blue", high = "red") + # Color gradient for data points
  
  # Data segments layer (whiskers)
  geom_segment(data = data, aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2,color = SlopeGesamt), size = 0.5) +  # Horizontal whisker for PC1
  geom_segment(data = data, aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = SlopeGesamt), size = 0.5) +  # Vertical whisker for PC2
  
  #geom_point(data = rec, aes(x = PC1, y = PC2), color = "black", fill = "transparent", size = 3) 
  # Text labels layer
  #geom_text_repel(data = data, 
  #              aes(label = Species, x = med_PC1, y = med_PC2),  # Specify x and y aesthetics
  #	box.padding = 0.5, 
  #            fontface = 'italic',
  #           size = 3, 
  #	   max.overlaps = Inf,  # Ensures all labels are drawn
  #          segment.color = 'black',  # Color of the connecting line for geom_text_repel
  #         segment.size = 0.5) +  # Size of the connecting line for geom_text_repel
  
  # Points layer from 'rec' dataframe
  #geom_point(data = rec.j, aes(x = PC1, y = PC2), color = "black", fill = "transparent", size = 3) + 
  
  # Labels and titles
  xlab("PC1 (71.0%)") +
  ylab("PC2 (17.33%)") +
  ggtitle("Slope of Change in Habitat Size without Habitat Mask") +
  
  # Ensure a proper theme
  theme_minimal()

# Display the combined plot
print(combined_plot)

##########

SlopeGesamt <- ggplot(data, aes(x = med_PC1, y = med_PC2 )) +
  geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = AverageSizeTotal), size = 0.01) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = AverageSizeTotal), size = 0.01)  # Vertical whisker for PC2
geom_point(aes(color = SlopeGesamt, shape = Category, size = Records),show.legend=FALSE) +  # Plotting the points with gradient color and increased size
  scale_color_gradient(low = "blue", high = "red") +
  #geom_text_repel(data = data, 
  #              aes(label = Species), 
  #             box.padding = 0.5, 
  #            fontface = c('italic'),
  #	  size = 3, 
  #          max.overlaps = Inf,  # Ensures all labels are drawn
  #   		segment.color = 'black',  # Color of the connecting line for geom_text_repel
  #        segment.size = 0.) + # Size of the connecting line for geom_text_repel
  xlab("PC1 (71.0%)")+
  ylab("PC2 (17.33%)")+
  ggtitle ("Slope of Change in Habitat Size without Habitat Mask")

geom_segment(aes(x = whisker_low_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = SlopeGesamt), size = 0.5) +  # Horizontal whisker for PC1
  geom_segment(aes(x = med_PC1, xend = med_PC1, y = whisker_low_PC2, yend = whisker_high_PC2, color = SlopeGesamt), size = 0.5) +  # Vertical whisker for PC2
  geom_segment(aes(x = whisker_low_PC1, xend = med_PC1, y = med_PC2, yend = med_PC2, color = SlopeGesamt), size = 0.5) +  # Connecting line from whisker_low_PC1 to med_PC1
  geom_segment(aes(x = med_PC1, xend = whisker_high_PC1, y = med_PC2, yend = med_PC2, color = SlopeGesamt), size = 0.5)+  # Connecting line from med_PC1 to whisker_high_PC1


################################################################################
# 3D PCA
################################################################################
wd <- "D:/.../Clim/bio4km"
setwd(wd)
#####
years <- 2002:2023

for (i in 1:length(years)){
  wd.i <- paste(wd, "/", years[i], sep="")
  setwd(wd.i)
  pc1 <- rast("PC1.asc")
  pc2 <- rast("PC2.asc")
  pc3 <- rast("PC3.asc")
  
  values1 <- values(pc1)
  values2 <- values(pc2)
  values3 <- values(pc3)
  
  values123 <- cbind(values1, values2, values3)
  values123 <- na.omit(values123)
  
  set.seed(123)  # Für Reproduzierbarkeit
  sampled_values <- values123[sample(nrow(values123), size = 10000, replace = FALSE), ]
  
  write.csv(sampled_values, "Random_PCA_Values.csv", row.names=F)
}
##### plot one year

my.files <- list.files(wd, pattern="Random_PCA_Values.csv", full.names = TRUE, recursive=T)
years <- 2002:2023

sampled_values <- read.csv(my.files[1])

kde <- kde3d(sampled_values$PC1, sampled_values$PC2, sampled_values$PC3, n = 100)
kde$d[kde$d <= min(kde$d)] <- 0
kde$d[kde$d < 1e-10] <- 0  
open3d()
plot3d(xlim <- c(-11, 5),xlab = "PC1 - 47.16%", ylim <- c(-8, 0),ylab = "PC2 - 19.43%", zlim <- c(-6, 4), zlab = "PC3 - 12.41%")
contour3d(kde$d, x= kde$x, y= kde$y, z= kde$z, 
          level = seq(min(kde$d), max(kde$d), length.out = 10), 
          color = tim.colors(15), alpha = 0.1, add=T)

axes3d(edges = "bbox", labels = TRUE, tick = TRUE, col = "black")


#####

# Starte eine neue 3D-Sitzung
my.files <- list.files(wd, pattern="Random_PCA_Values.csv", full.names = TRUE, recursive=T)
data <- read.csv("D:/.../PCA_Hypervolume/PCA_Stats/All_Stats_Merged_prep_GLM.csv")
years <- 2002:2023
for (k in 22:length(my.files)) {
  # kde_data <- read.csv(my.files[k])
  kde_data <- read.csv("D:/.../PCA_Hypervolume/Clim/bio4km/Average/Random_PCA_Values.csv") #for average 
  year <- years[k]
  kde <- kde3d(kde_data[,1], kde_data[,2], kde_data[,3], n = 100) 
  kde$d[kde$d < 1e-10] <- 0                          # set very small values 0
  kde$d[kde$d < 0] <- 0         
  
  open3d()
  par3d(windowRect = c(0, 0, 700, 700))
  plot3d() #xlim <- c(-11, 5),xlab = "PC1 - 47.16%", ylim <- c(-8, 0),ylab = "PC2 - 19.43%", zlim <- c(-6, 4), zlab = "PC3 - 12.41%")
  #title3d(main = paste(years[k]), line = 5, cex = 1.0)
  title3d(main = "Average", line = 5, cex = 1.0) #for average
  contour3d(kde$d, kde$x, kde$y, kde$z, 
            level = seq(min(kde$d), max(kde$d), length.out = 10), 
            color = tim.colors(15), alpha = 0.1, add=T)
  #points3d(kde_data, color="grey")	
  
  
  ########## plot PC and quantiles ######
  
  max_size <- 0.8
  min_size <- 0.2
  
  # calculate size
  data$scaled_size <- min_size + (data$Records / max(data$Records)) * (max_size - min_size)
  
  # add colors, based on SlopeMaskrm
  color_map <- colorRampPalette(c("blue", "red"))(length(data$SlopeMaskrm))
  data$color <- color_map[rank(data$SlopeMaskrm)]
  
  # plot without data with correct size 
  plot3d(data$med_PC1_year, data$med_PC2_year, data$med_PC3_year,
         col = data$color, xlab = "PC1", ylab = "PC2", zlab = "PC3",
         main = "3D Plot with Quantiles",add=F, size=8)
  
  # add median with size weighted by number of records 
  spheres3d(data$med_PC1_year, data$med_PC2_year, data$med_PC3_year,
            radius = data$scaled_size / 2,
            col = data$color)
  
  for (i in 1:38) {
    
    # add Species 
    segments3d(rbind(c(data$med_PC1_year[i], data$med_PC2_year[i], data$med_PC3_year[i]),
                     c(data$med_PC1_year[i], data$med_PC2_year[i], data$med_PC3_year[i] - 0.1)), lwd = 1.3, col = "black")
    
    text3d(data$med_PC1_year[i], data$med_PC2_year[i], data$med_PC3_year[i] - 0.12,  
           texts = data$Species[i], 
           col = "black",  # Farbe der Beschriftungen
           cex = 0.8, 
           font= 2)
  }
  for (i in 19:38) {
    
    # add Species 
    segments3d(rbind(c(data$med_PC1_year[i], data$med_PC2_year[i], data$med_PC3_year[i]),
                     c(data$med_PC1_year[i], data$med_PC2_year[i], data$med_PC3_year[i] + 0.1)), lwd = 1.3, col = "black")
    
    text3d(data$med_PC1_year[i], data$med_PC2_year[i], data$med_PC3_year[i] + 0.12,  
           texts = data$Species[i], 
           col = "black",  # Farbe der Beschriftungen
           cex = 0.8, 
           font= 2)
  }
}		 
#Quantiles PC1
segments3d(x = c(data$PC1_lower_q_year[i], data$PC1_upper_q_year[i]),
           y = rep(data$med_PC2_year[i], 2),
           z = rep(data$med_PC3_year[i], 2),
           col = data$color[i], lwd = 0.9, alpha = 0.5)

# Quantiles for PC2
segments3d(x = rep(data$med_PC1_year[i], 2),
           y = c(data$PC2_lower_q_year[i], data$PC2_upper_q_year[i]),
           z = rep(data$med_PC3_year[i], 2),
           col = data$color[i], lwd = 0.9, alpha = 0.5)

# Quantiles for PC3
segments3d(x = rep(data$med_PC1_year[i], 2),
           y = rep(data$med_PC2_year[i], 2),
           z = c(data$PC3_lower_q_year[i], data$PC3_upper_q_year[i]),
           col = data$color[i], lwd = 0.9, alpha = 0.5)
}
movie3d(spin3d(axis = c(1,1,1), rpm = 4), duration = 15, fps = 25, convert=NULL, movie= paste0(years[k], "_"), type=".gif", dir = "D:/.../PCA_Hypervolume/PCA_Stats/Test")
}


### plot PC Density and the PC_Medians and PC Quantiles for the Static Data ####

sampled_values <- read.csv("D:/.../PCA_Hypervolume/Clim/bio4km/Average/Random_PCA_Values.csv") 
kde <- kde3d(sampled_values$PC1, sampled_values$PC2, sampled_values$PC3, n = 100)
kde$d[kde$d <= min(kde$d)] <- 0
kde$d[kde$d < 1e-10] <- 0  

open3d()
plot3d(xlim <- c(-11, 5),xlab = "PC1 (47.16%)", ylim <- c(-8, 0),ylab = "PC2 (19.43%)", zlim <- c(-6, 4), zlab = "PC3 (12.41%)", line = 3, cex= 0.8)
title3d(main = "Average", line = 5, cex = 0.8)
contour3d(kde$d, kde$x, kde$y, kde$z, 
          level = seq(min(kde$d), max(kde$d), length.out = 10), 
          color = tim.colors(15), alpha = 0.1, add=T)


########## plot PC and whiskers ######

####alternative####
data <- read.csv("D:/.../PCA_Hypervolume/PCA_Stats/All_Stats_Merged_prep_GLM.csv")

max_size <- 0.8
min_size <- 0.2

# calculate size
data$scaled_size <- min_size + (data$Records / max(data$Records)) * (max_size - min_size)

# add colors, based on SlopeMaskrm
color_map <- colorRampPalette(c("blue", "red"))(length(data$SlopeMaskrm))
data$color <- color_map[rank(data$SlopeMaskrm)]

# plot without data with correct size 
plot3d(data$med_PC1_av, data$med_PC2_av, data$med_PC3_av,
       col = data$color, xlab = "PC1", ylab = "PC2", zlab = "PC3",
       main = "3D Plot with Quantiles",add=T)

# add median with size weighted by number of records 
spheres3d(data$med_PC1_av, data$med_PC2_av, data$med_PC3_av,
          radius = data$scaled_size / 2,
          col = data$color)

for (i in 1:nrow(data)) {
  
  # add Species 
  #segments3d(rbind(c(data$med_PC1_av[i], data$med_PC2_av[i], data$med_PC3_av[i]),
  #  c(data$med_PC1_av[i], data$med_PC2_av[i], data$med_PC3_av[i] + 0.3)), lwd = 1.3, col = "black")
  
  #text3d(data$med_PC1_av[i], data$med_PC2_av[i], data$med_PC3_av[i] + 0.35,  
  #        texts = data$Species[i], 
  #         col = "black",  # Farbe der Beschriftungen
  #        cex = 0.9)
  
  #Quantiles PC1
  segments3d(x = c(data$PC1_lower_q_av[i], data$PC1_upper_q_av[i]),
             y = rep(data$med_PC2_av[i], 2),
             z = rep(data$med_PC3_av[i], 2),
             col = data$color[i], lwd = 0.9, alpha = 0.5)
  
  # Quantiles for PC2
  segments3d(x = rep(data$med_PC1_av[i], 2),
             y = c(data$PC2_lower_q_av[i], data$PC2_upper_q_av[i]),
             z = rep(data$med_PC3_av[i], 2),
             col = data$color[i], lwd = 0.9, alpha = 0.5)
  
  # Quantiles for PC3
  segments3d(x = rep(data$med_PC1_av[i], 2),
             y = rep(data$med_PC2_av[i], 2),
             z = c(data$PC3_lower_q_av[i], data$PC3_upper_q_av[i]),
             col = data$color[i], lwd = 0.9, alpha = 0.5)
}
movie3d(spin3d(axis = c(1,1,1), rpm = 4), duration = 15, fps = 25, convert=NULL, movie= paste0("Average", "_"), type=".gif", dir = "D:/.../PCA_Hypervolume/PCA_Stats/Test")
}


### PCA Variable Plot with 3 PC
dat <- read.csv("D:/.../PCA_Hypervolume/summary.PCA.all.csv")
dat <- dat[-c(10:11),]

plot3d(dat$Comp.1, dat$Comp.2, dat$Comp.3,
       xlim = c(-1, 1), ylim = c(-1, 1), zlim = c(-1, 1),
       type = "p", col = "black", size = 4, 
       xlab = "PC1 (47.16%)", ylab = "PC2 (19.43%)", zlab = "PC3 (12.41%)")
text3d(dat$Comp.1, dat$Comp.2, dat$Comp.3,
       texts = dat$X, adj = 1.2, cex = 1.2)
for (i in 1:nrow(dat)) {
  arrow3d(p0 = c(0, 0, 0),  # Start at the origin
          p1 = c(dat$Comp.1[i], dat$Comp.2[i], dat$Comp.3[i]), 
          col = "darkgrey", lwd = 0.5, barblen = 0.02, type = "lines")
}
