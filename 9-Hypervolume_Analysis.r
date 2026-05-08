
# Title:   Hypervolume analysis for herpetufauna of Romania
# Author:   Laura Karolin Steib and Dennis Rödder
# Date:     09.01.2025

library(sp)
library (sf)
library(rSDM)
library(terra)
library(raster)



########################### second try ###############################
wd <- "D:/.../hypervolumes/zSec_Try"
setwd(wd)
  
################# get datasets containing background points (one species, as all the same) 
 
 ## for static data 
bg <- read.csv("D:/.../hypervolumes/average/background_Ablepharus_kitaibelii.csv")
bg$ID <- paste("bg_av.")
out <- bg 
for (j in 1:length(sp)){
  if (j==37) next 
  sp.j <- read.csv(paste("D:/.../hypervolumes", "/average/", "records_", sp[j], ".csv", sep=""))
  sp.j$ID <- paste("rec_av.",  sp[j], sep="")
  out <- rbind(out,sp.j)
}
write.csv(out, "static_data_sec.csv")
  
  ## for dynamic data  
bg <- read.csv("D:/.../hypervolumes/singleyears/background_Ablepharus_kitaibelii.csv")
bg$ID <- paste("bg_years.")
out <- bg 
for (j in 1:length(sp)){
  if (j==37) next 
  sp.j <- read.csv(paste("D:/.../hypervolumes", "/singleyears/", "records_", sp[j], ".csv", sep=""))
  sp.j$ID <- paste("rec_years.",  sp[j], sep="")
  out <- rbind(out,sp.j)
}
write.csv(out, "dynamic_data_sec.csv")

### conduct PCA for static and dynamic model ###

  ## for static data
static <- read.csv("static_data_sec.csv")
dynamic <- read.csv("dynamic_data_sec.csv")
combined <- rbind(static,dynamic)

write.csv(combined, "all.csv", row.names = F)
combined <- static
combined <- read.csv("all.csv", h = T) # read in data 
combined_filtered <- combined[,4:12]          # select only columns containig clim cariables 

write.csv(combined_filtered, "all_raw.pca.csv", row.names =F)

pca <- princomp(combined_filtered, cor=T)             #perform PCA 
eig <- pca$sdev^2                                   #compute eigenvalues
anz.PC=sum(as.vector(eig)>=1 )           #number of PCs with Eigenvalues > 1 and with value to get 3PC
                
exp.variance <- c()                                 # compute explained variance per PC
for (x in 1:length(eig)) {                        
    exp.variance.x <- eig[x]/sum(as.vector(eig))*100
    exp.variance <- cbind(exp.variance, exp.variance.x)
}

factor_loadings=c()
factor_loadings<-cbind(factor_loadings, colnames(combined_filtered))
yy <- data.frame(pca$scores[, 1:anz.PC])

for(i in 1:ncol(yy)){
  xPC=c()
  for(k in 1:ncol(combined_filtered)){
   xPC=c(xPC, cor.test(combined_filtered[,k],yy[,i])$estimate)    
  }
  factor_loadings=cbind(factor_loadings, xPC)
}
colnames(factor_loadings)<- c("", colnames(yy))

factor_loadings=rbind(factor_loadings, c("Eigenvalues" ,eig[1:anz.PC]), c("Explained Variance", exp.variance[1:anz.PC]))
write.csv(factor_loadings, "summary.PCA.all.csv", row.names = FALSE)

pc_scores <- pca$scores[,1:anz.PC]             #get the PC for each data point     
pc_output <- cbind(combined,pc_scores)
write.csv(pc_output, "all_with_PC.csv", row.names = FALSE)



############ calculate Hypervolume of Background & Records for static & dynamic and calculate overlaps
 wd <- "D:/.../hypervolumes/zSec_Try"
setwd(wd)
  
 library(hypervolume)
 library(dplyr)
 library(rgl)
 library(alphahull)
 
 data <- read.csv("all_with_PC.csv", h=T)
 background <- subset(data, data$Species == "Backgound")
 records <- subset(data, data$Species != "Backgound")
 static_bg <- subset(background, background$ID == "bg_av.")
 static_rec <- records %>% filter(grepl("rec_av.", ID))
 dynamic_bg <- subset(background, background$ID == "bg_years.")
 dynamic_rec <- records %>% filter(grepl("rec_years.", ID))
 
  #background static - 3PC and 2PC
static_bg <- na.omit(static_bg[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.static_bg <- hypervolume_svm(static_bg, name = "Static_Background") # use the PCA results of the background point of the static SDM to get the hypervolume (fundamental niche)
vol.static_bg <- hv.static_bg@Volume

  #records static - 3PC and 2PC
static_rec <- na.omit(static_rec[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.static_rec <- hypervolume_svm(static_rec, name = "Static_Records") # use the PCA results of the background point of the static SDM to get the hypervolume (fundamental niche)
vol.static_rec <- hv.static_rec@Volume
  
  #background dynamic - 3PC and 2PC
dynamic_bg <- na.omit(dynamic_bg[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.dynamic_bg <- hypervolume_svm(dynamic_bg, name = "Dynamic_Background") # use the PCA results of the background point of the dynamic SDM to get the hypervolume (fundamental niche)
vol.dynamic_bg <- hv.dynamic_bg@Volume

  #records dynamic - 3PC and 2PC
dynamic_rec <- na.omit(dynamic_rec[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.dynamic_rec <- hypervolume_svm(dynamic_rec, name = "Dynamic_Records") # use the PCA results of the background point of the dynamic SDM to get the hypervolume (fundamental niche)
vol.dynamic_rec <- hv.dynamic_rec@Volume
 
  # get all the niche volumes 
volumes <- rbind(vol.dynamic_bg,vol.static_bg, vol.dynamic_rec, vol.static_rec)
write.csv(volumes, "hypervolumes.csv")

  #overlap static and dynamic background
hv.set_bg <- hypervolume_set(hv.dynamic_bg, hv.static_bg, check.memory=FALSE)
hv.overlap_bg_stats <- hypervolume_overlap_statistics(hv.set_bg)            # get the niche overlap of the two models 

write.csv(hv.overlap_bg_stats, "hypervolumeStats_background.csv")           # save all the stuff 
pdf("hypervolume_overlap_background.pdf")
plot(hv.set_bg)
dev.off()

#overlap static and dynamic records
hv.set_rec <- hypervolume_set(hv.dynamic_rec, hv.static_rec, check.memory=FALSE)
hv.overlap_rec_stats <- hypervolume_overlap_statistics(hv.set_rec)            # get the niche overlap of the two models 

write.csv(hv.overlap_rec_stats, "hypervolumeStats_records.csv")           # save all the stuff 
pdf("hypervolume_overlap_records.pdf")
plot(hv.set_rec)
dev.off()

#overlap static - background and records 
hv.set_static <- hypervolume_set(hv.static_bg, hv.static_rec, check.memory=FALSE)
hv.overlap_static_stats <- hypervolume_overlap_statistics(hv.set_static)            # get the niche overlap of the two models 

write.csv(hv.overlap_static_stats, "hypervolumeStats_static.csv")           # save all the stuff 
pdf("hypervolume_overlap_static.pdf")
plot(hv.set_static)
dev.off()

#overlap dynamic - background and records 
hv.set_dynamic <- hypervolume_set(hv.dynamic_bg, hv.dynamic_rec, check.memory=FALSE)
hv.overlap_dynamic_stats <- hypervolume_overlap_statistics(hv.set_dynamic)            # get the niche overlap of the two models 

write.csv(hv.overlap_dynamic_stats, "hypervolumeStats_dynamic.csv")           # save all the stuff 
pdf("hypervolume_overlap_dynamic.pdf")
plot(hv.set_dynamic)
dev.off()

# plot 3d settings
hv <- hypervolume_join(hv.static_bg, hv.dynamic_bg, hv.static_rec, hv.dynamic_rec)
hv.back <- hypervolume_join(hv.static_bg, hv.dynamic_bg)
hv.rec <- hypervolume_join(hv.static_rec, hv.dynamic_rec)
hv.stat <- hypervolume_join(hv.static_bg, hv.static_rec)
hv.dyn <- hypervolume_join(hv.dynamic_bg, hv.dynamic_rec)

plot(hv.back, show.3d=F, show.random=F, show.density=T,show.data=T, colors= c("grey","black", "blue", "red"))

plot(hv.dyn, show.3d=TRUE, show.random=F, show.density=F,show.data=T, colors= c( "blue", "red"),num.points.max.data = 1000, num.points.max.random = 1000, cex.random=3,cex.data=3)
  hypervolume_save_animated_gif(image.size = 800, axis = c(1, 1, 1), duration = 30,file.name = "dynamic")
  rgl.close()


################################## get hypervolume of records for each species 
data <- read.csv("all_with_PC.csv", h=T)
records <- subset(data, data$Species != "Backgound")
 
dynamic <- records %>% filter(grepl("rec_years.", ID))
static <- records %>% filter(grepl("rec_av.", ID))

sp <- read.csv("D:/.../Species.csv")
 sp <- names(table(sp$Species))
 sp <- sp[-37]

selected <- c(1, 2, 6, 8, 15, 21, 27, 31, 33, 34, 37)

out.vol <- c()
out.stats <- c()

for (j in selected){
dynamic.j <- subset(dynamic, dynamic$Species == sp[j])
dynamic.j <- na.omit(dynamic.j[,c('Comp.1', 'Comp.2', 'Comp.3')])
hv.dynamic.j <- hypervolume_svm(dynamic.j, name = paste("Dynamic", sp[j]))
vol.dynamic.j <- hv.dynamic.j@Volume

static.j <- subset(static, static$Species == sp[j])
static.j <- na.omit(static.j[,c('Comp.1', 'Comp.2', 'Comp.3')])
hv.static.j <- hypervolume_svm(static.j, name = paste("Static", sp[j]))
vol.static.j <- hv.static.j@Volume

hv.set.j <- hypervolume_set(hv.dynamic.j,hv.static.j, check.memory=FALSE)
hv.overlap_stats.j <- hypervolume_overlap_statistics(hv.set.j)
hv.join.j <- hypervolume_join(hv.dynamic.j,hv.static.j)

plot(hv.join.j, show.3d=TRUE, show.random=F, show.density=F,show.data=T, colors= c( "blue", "red", "grey", "black"),num.points.max.data = 1000, num.points.max.random = 1000, cex.random=3,cex.data=3)
  hypervolume_save_animated_gif(image.size = 800, axis = c(1, 1, 1), duration = 30,file.name = sp[j])
  rgl.close()
}
vol.j <- as.data.frame(cbind(vol.dynamic.j, vol.static.j))
vol.j$Species <- paste(sp[j])

hv.overlap_stats.j <- as.list(hv.overlap_stats.j)
hv.overlap_stats.j <- as.data.frame(hv.overlap_stats.j)
hv.overlap_stats.j$Species <- paste(sp[j])

out.vol <- rbind(out.vol,vol.j)
out.stats <- rbind(out.stats, hv.overlap_stats.j)
}
write.csv(out.vol, "hypervolume_per_species.csv")
write.csv(out.stats, "hypervolume_Stats_per_species.csv")



###### redo PCA with following Hypervolume ########
v <- vect("D:/.../ROU_adm0.shp")
setwd(wd)
library(raster)
library(dismo)
library(shapefiles)
library(vegan)

###############################################################################
# Load grids from folder
###############################################################################
#grids <- list.files(path=wd, pattern='asc', full.names=T)
#predictors <- stack(grids)
#names(predictors)

############ prepare the clim layers ##################################


################ run PCA ###########################################
wd <- "D:/.../PCA_Hypervolume"
setwd(wd)
combined <- read.csv("all_sortiert.csv", h = T) # read in data 
combined_filtered <- combined[,4:12]          # select only columns containig clim cariables 

write.csv(combined_filtered, "all_raw.pca.csv", row.names =F)

pca <- princomp(combined_filtered, cor=T)             #perform PCA 
#pca <- prcomp(combined_filtered, scale = T)
eig <- pca$sdev^2                                   #compute eigenvalues
anz.PC=sum(as.vector(eig)>=1 )           #number of PCs with Eigenvalues > 1 and with value to get 3PC
                
exp.variance <- c()                                 # compute explained variance per PC
for (x in 1:length(eig)) {                        
    exp.variance.x <- eig[x]/sum(as.vector(eig))*100
    exp.variance <- cbind(exp.variance, exp.variance.x)
}

factor_loadings=c()
factor_loadings<-cbind(factor_loadings, colnames(combined_filtered))
yy <- data.frame(pca$scores[, 1:anz.PC])

for(i in 1:ncol(yy)){
  xPC=c()
  for(k in 1:ncol(combined_filtered)){
   xPC=c(xPC, cor.test(combined_filtered[,k],yy[,i])$estimate)    
  }
  factor_loadings=cbind(factor_loadings, xPC)
}
colnames(factor_loadings)<- c("", colnames(yy))

factor_loadings=rbind(factor_loadings, c("Eigenvalues" ,eig[1:anz.PC]), c("Explained Variance", exp.variance[1:anz.PC]))
write.csv(factor_loadings, "summary.PCA.all.csv", row.names = FALSE)

pc_scores <- pca$scores[,1:anz.PC]             #get the PC for each data point     
pc_output <- cbind(combined,pc_scores)
write.csv(pc_output, "all_with_PC.csv", row.names = FALSE)
saveRDS(pca, file = "pca_model.rds")

#pca <- readRDS("pca_model.rds")

# Project PCA and write output in asciis
###############################################################################
wd <- "D:/.../PCA_Hypervolume"
setwd(wd)

for (i in 6:length(years)){
setwd(wd)
pca <- readRDS("D:/.../PCA_Hypervolume/pca_model.rds")
eig <- pca$sdev^2
anz.PC=sum(as.vector(eig)>=1 ) 
years <- 2002:2023
years <- as.list(years)
PCs <- c()
# wd.i <- paste ("D:/.../PCA_Hypervolume/Clim/bio4km/", years[i], sep="") # for years
wd.i <- paste("D:/.../PCA_Hypervolume/Clim/bio4km/Average")
setwd(wd.i)
grids <- list.files(path=wd.i, pattern='.asc', full.names=T)
predictors <- stack(grids)
names(predictors)
PCs <- predict(predictors, pca, index=1:anz.PC)          
  for (k in 1:anz.PC) {
    myfilename <- paste("PC", k, ".asc", sep="")
    writeRaster(PCs[[k]], filename = myfilename, datatype='FLT4S', overwrite=TRUE)
}
options(scipen=99)
}

# calculate a density plot of the PCA

wd <- "D:/.../PCA_Hypervolume/Clim/bio4km/2002"
setwd(wd)

library(terra)
library(MASS)
library(raster)
library(fields)
library(rgl)
library(plotly)

pc1 <- rast("PC1.asc")
pc2 <- rast("PC2.asc")
pc3 <- rast("PC3.asc")

values1 <- values(pc1)
values2 <- values(pc2)
values3 <- values(pc3)
values123 <- cbind(values1, values2,values3)
values123 <- na.omit(values123)
set.seed(123)
values123_sampled <- values123[sample(nrow(values12),replace=F,size=100000),]


kde <- kde2d(values123[,1], values123[,2], values123[,3]n = 1000)
image(kde, col=tim.colors(100), ylab = "PC2 (17.33%)", xlab= "PC1 (71.0%)", zlab ="PC3")
rast.kde <- raster(kde)
kde.values <- extract(rast.kde, values123[,1:3])

rast.1 <- rast.kde
rast.1[rast.1 < min(kde.values)] <- NA
plotly(x=kde$x,y=kde$y,z=kde$z, col=tim.colors(100))

writeRaster(rast.1, "PCA.density.tif")


############ calculate Hypervolume of Background & Records for static & dynamic and calculate overlaps
 wd <- "D:/.../PCA_Hypervolume"
setwd(wd)
  
 library(hypervolume)
 library(dplyr)
 library(rgl)
 library(alphahull)
 
 data <- read.csv("all_with_PC.csv", h=T)
 background <- subset(data, data$Species == "Backgound")
 records <- subset(data, data$Species != "Backgound")
 static_bg <- subset(background, background$ID == "bg_av")
 static_rec <- records %>% filter(grepl("rec_av", ID))
 dynamic_bg <- subset(background, background$ID == "bg_years")
 dynamic_rec <- records %>% filter(grepl("rec_years", ID))
 
  #background static 
static_bg <- na.omit(static_bg[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.static_bg <- hypervolume_svm(static_bg, name = "Static_Background", scale.factor =1) # use the PCA results of the background point of the static SDM to get the hypervolume (fundamental niche)
vol.static_bg <- hv.static_bg@Volume

  #records static 
static_rec <- na.omit(static_rec[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.static_rec <- hypervolume_svm(static_rec, name = "Static_Records") # use the PCA results of the background point of the static SDM to get the hypervolume (fundamental niche)
vol.static_rec <- hv.static_rec@Volume
  
  #background dynamic 
dynamic_bg <- na.omit(dynamic_bg[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.dynamic_bg <- hypervolume_svm(dynamic_bg, name = "Dynamic_Background") # use the PCA results of the background point of the dynamic SDM to get the hypervolume (fundamental niche)
vol.dynamic_bg <- hv.dynamic_bg@Volume

  #records dynamic 
dynamic_rec <- na.omit(dynamic_rec[,c('Comp.1', 'Comp.2', 'Comp.3')])

hv.dynamic_rec <- hypervolume_svm(dynamic_rec, name = "Dynamic_Records") # use the PCA results of the background point of the dynamic SDM to get the hypervolume (fundamental niche)
vol.dynamic_rec <- hv.dynamic_rec@Volume
 
  # get all the niche volumes 
volumes <- rbind(vol.dynamic_bg,vol.static_bg, vol.dynamic_rec, vol.static_rec)
write.csv(volumes, paste(wd, "/hypervol_results/","hypervolumes.csv"))

  #overlap static and dynamic background
hv.set_bg <- hypervolume_set(hv.dynamic_bg, hv.static_bg, check.memory=FALSE)
hv.overlap_bg_stats <- hypervolume_overlap_statistics(hv.set_bg)            # get the niche overlap of the two models 

write.csv(hv.overlap_bg_stats, paste(wd, "/hypervol_results/", "hypervolumeStats_background.csv", sep=""))           # save all the stuff 
pdf(paste(wd, "/hypervol_results/","hypervolume_overlap_background.pdf", sep=""))
plot(hv.set_bg)
dev.off()

#overlap static and dynamic records
hv.set_rec <- hypervolume_set(hv.dynamic_rec, hv.static_rec, check.memory=FALSE)
hv.overlap_rec_stats <- hypervolume_overlap_statistics(hv.set_rec)            # get the niche overlap of the two models 

write.csv(hv.overlap_rec_stats, paste(wd, "/hypervol_results/", "hypervolumeStats_records.csv",sep=""))           # save all the stuff 
pdf(paste(wd, "/hypervol_results/","hypervolume_overlap_records.pdf",sep=""))
plot(hv.set_rec)
dev.off()

#overlap static - background and records 
hv.set_static <- hypervolume_set(hv.static_bg, hv.static_rec, check.memory=FALSE)
hv.overlap_static_stats <- hypervolume_overlap_statistics(hv.set_static)            # get the niche overlap of the two models 

write.csv(hv.overlap_static_stats,paste(wd, "/hypervol_results/", "hypervolumeStats_static.csv",sep=""))           # save all the stuff 
pdf(paste(wd, "/hypervol_results/","hypervolume_overlap_static.pdf",sep=""))
plot(hv.set_static)
dev.off()

#overlap dynamic - background and records 
hv.set_dynamic <- hypervolume_set(hv.dynamic_bg, hv.dynamic_rec, check.memory=FALSE)
hv.overlap_dynamic_stats <- hypervolume_overlap_statistics(hv.set_dynamic)            # get the niche overlap of the two models 

write.csv(hv.overlap_dynamic_stats, paste(wd, "/hypervol_results/","hypervolumeStats_dynamic.csv",sep=""))           # save all the stuff 
pdf(paste(wd, "/hypervol_results/","hypervolume_overlap_dynamic.pdf",sep=""))
plot(hv.set_dynamic)
dev.off()

# plot 3d settings
hv <- hypervolume_join(hv.static_bg,hv.dynamic_bg, hv.static_rec, hv.dynamic_rec)
hv.back <- hypervolume_join(hv.static_bg, hv.dynamic_bg)
hv.rec <- hypervolume_join(hv.static_rec, hv.dynamic_rec)
hv.stat <- hypervolume_join(hv.static_bg, hv.static_rec)
hv.dyn <- hypervolume_join(hv.dynamic_bg, hv.dynamic_rec)

wd.save <- "D:/.../PCA_Hypervolume/hypervol_results"

     #save as movie
plot(hv, show.3d=T, show.random=F, show.density=F,show.data=F, colors= c( "red", "blue", "green", "purple"),num.points.max.data = 1000, num.points.max.random = 1000, cex.random=3,cex.data=4)
  hypervolume_save_animated_gif(image.size = 800, axis = c(1, 1, 1), duration = 30,file.name = "all")
  rgl.close()
  
    #save as pdf
pdf(paste(wd.save,"/", "hypervolume_overlap_all.pdf", sep=""))
plot(hv, show.3d=F, show.random=T, show.density=F,show.data=F, colors= c("red", "blue", "green", "purple"),num.points.max.data = 1000, num.points.max.random = 1000, cex.random=0.7,cex.data=0.7)
dev.off()
	
	
	plot(hv, show.3d=F, show.random=F, show.density=F, show.data=F, 
     colors=c("red", "blue", "green", "purple"), 
     num.points.max.data = 2000, num.points.max.random = 2000, 
     cex.random=0.7, cex.data=0.7,
	 xlim=c(-5, 3), ylim=c(-3, 5))
	 
################################ get hypervolume of each species #####################
wd <- "D:/.../PCA_Hypervolume"
setwd(wd)
data <- read.csv("all_with_PC.csv", h=T)
records <- subset(data, data$Species != "Backgound")
sp <- unique(records$Species)

out <- c()
for (j in 1:length(sp)){
sp.j <- subset(records, records$Species == sp[j])
sp.j.year <- subset(sp.j, sp.j$ID == "rec_years")
sp.j.av <- subset(sp.j, sp.j$ID == "rec_av")
sp.j.year <- na.omit(sp.j.year[,c('Comp.1', 'Comp.2', 'Comp.3')])
sp.j.av <- na.omit(sp.j.av[,c('Comp.1', 'Comp.2', 'Comp.3')])
hv.sp.j.year <- hypervolume_svm(sp.j.year, name = paste(sp[j], "year", sep =""))
hv.sp.j.av <- hypervolume_svm(sp.j.av, name = paste(sp[j], "av", sep =""))
vol_year <- hv.sp.j.year@Volume
vol_av <- hv.sp.j.av@Volume 
set <- hypervolume_set(hv.sp.j.year, hv.sp.j.av, check.memory=F)
stats <- as.data.frame(t(hypervolume_overlap_statistics(set)))
out.j <- as.data.frame(cbind(vol_year,vol_av,stats))
out.j$Species <- sp[j]
out <- rbind(out,out.j)
}

## plot frequency distribution of climate ## 
wd <- "D:/.../PCA_Hypervolume"
setwd(wd)

library(sm)
library(dismo)
library(tidyr)
library(RColorBrewer)


records <- read.csv("all_with_PC.csv", h=T)

## add 3 and 5 year average ##
background_3y <- read.csv("D:/.../PCA_Hypervolume/Clim/bio4km/background_3y_avg_withPC.csv",h=T)
background_5y <- read.csv("D:/.../PCA_Hypervolume/Clim/bio4km/background_5y_avg_withPC.csv",h=T)
records_3y <- read.csv("D:/.../PCA_Hypervolume/Clim/bio4km/rec_3_year_avg_withPC.csv",h=T)
records_5y <- read.csv("D:/.../PCA_Hypervolume/Clim/bio4km/rec_5_year_avg_withPC.csv",h=T)

comb <- rbind(records, background_3y, background_5y, records_3y, records_5y)

write.csv(comb, "all_with_PC_and_3_5_y_avg.csv", row.names=F)
#grids <- list.files(path=getwd(), pattern='asc', full.names=T)
#predictors <- stack(grids)

#scores <- extract(predictors, records[,2:3])
#scores <- cbind(records, scores)
records <- read.csv("all_with_PC_and_3_5_y_avg.csv")
mydf_comb <- na.omit(records)

# density estmimates - frequency distribution climate 

#  for background 
mydf_comb1 <- subset(mydf_comb, mydf_comb$Species == "Backgound")

Species <- names(table(mydf_comb1$ID))
mycolnames <- colnames(mydf_comb1)
 
mydf_comb1$ID <- factor(mydf_comb1$ID, levels = c("bg_av", "bg_years", "bg_3y_avg", "bg_5y_avg"))
colors <- c(1:4)
new_names <- c("Static", "Dynamic - Yearly")# "Dynamic - 3 Years", "Dynamic - 5 Years")

pdf(paste(wd, "/Climate_dis/", "freqPlot_back.pdf", sep=""))
par(mfrow=c(5,3), mar=c(4,4,1,1))
for (i in 4:15) {
    sm.density.compare(mydf_comb1[,i], mydf_comb1$ID, col = colors, xlab=mycolnames[i], model="none", fill=colors, lwd=1.2)
    }
	  # Eindeutige Werte von 'ID'
	  plot.new()
	  mtext("Background", side=3, line=0, cex=0.8)
legend("center", legend=new_names, fill=colors,border=NULL)
dev.off()


# for records
mydf_comb2 <- subset(mydf_comb, mydf_comb$Species != "Backgound")

Species <- names(table(mydf_comb2_sep$ID))
mycolnames <- colnames(mydf_comb2_sep)

mydf_comb2$ID <- factor(mydf_comb2$ID, levels = c("rec_av", "rec_years", "3y_avg", "5y_avg"))
colors <- c(1:4)
new_names <- c("Static", "Dynamic - Yearly", "Dynamic - 3 Years", "Dynamic - 5 Years")

pdf(paste(wd, "/Climate_dis/", "freqPlot_recs.pdf", sep=""))
par(mfrow=c(5,3), mar=c(4,4,1,1))
for (i in 4:15) {
    sm.density.compare(mydf_comb2[,i], mydf_comb2$ID, col = colors , xlab=mycolnames[i], model="none",fill=colors, lwd=1.2)
    }
	plot.new()
	mtext("Records", side=3, line=0, cex=0.8)
legend("center", legend=new_names, fill=colors,border=NULL)
dev.off()


# for each species 

mydf_comb2 <- subset(mydf_comb, mydf_comb$Species != "Backgound")

mydf_comb2$ID <- factor(mydf_comb2$ID, levels = c("rec_av", "rec_years", "3y_avg", "5y_avg"))
colors <- c(1:4)
new_names <- c("Static", "Dynamic - Yearly", "Dynamic - 3 Years", "Dynamic - 5 Years")

sp <- names(table(mydf_comb2$Species))

for (j in 1:length(sp)){
sub.j <- subset(mydf_comb2, mydf_comb2$Species == sp[j])
pdf(paste(wd, "/Climate_dis/", "freqPlot_", sp[j], ".pdf", sep=""))
par(mfrow=c(5,3), mar=c(4,4,1,1))
for (i in 4:15) {
    sm.density.compare(sub.j[,i], sub.j$ID, col = colors, xlab=mycolnames[i], model="none",fill=colors, lwd=1.2)
    }
	plot.new()
	mtext(paste(sp[j]), side=3, line=0, cex=0.8)
legend("center", legend=new_names, fill=colors, border=NULL)
dev.off()
}


###### prepare PCA plot with PCs of species #####################
# get the boxplot stats for PC1/PC2/PC3
library(tidyr)
wd <- "D:/.../PCA_Hypervolume"
setwd(wd)

data <- read.csv("all_with_PC.csv")
data <- subset(data,data$Species != "Backgound")
data$ID2 <- data$ID                              # Datensatz nach dynamic und static aufteilen 
data <- data %>%                                 #
separate(ID2, into = c("lev", "sp_"), sep="\\.")  #
data <- subset(data,data$lev == "rec_years")        # rec_av wenn static data und rec_years wenn dynamic 
data <- data [, c(1:3,13:16)]
sp <- names(table(data$Species))

out2 <- c()
for (j in 1:length(sp)){
rec.j <- subset(data, data$Species == sp[j])
rec.j <- as.data.frame(rec.j [, "Comp.3"])
colnames (rec.j) [1] <- "Comp.3"
out <- c()
Species <- sp[j]
out <- cbind(out,Species)
med_PC3_year <- median(rec.j$Comp.3)
mean_PC3_year <- mean(rec.j$Comp.3)
out <- cbind(out,mean_PC3_year, med_PC3_year)
sd_PC3_year <- sd(rec.j$Comp.3)
out <- cbind(out,sd_PC3_year)
t_value <- qt(0.95, nrow(rec.j)-1)              # Nutzung der t-Verteilung 
SE <- sd_PC3_year/sqrt(nrow(rec.j))
CI_lower_PC3_year <- mean (rec.j$Comp.3) - t_value *SE
CI_upper_PC3_year <- mean (rec.j$Comp.3) + t_value *SE
out <- cbind(out, CI_lower_PC3_year,CI_upper_PC3_year)
Q1 <- quantile(rec.j$Comp.3, 0.25)
Q3 <- quantile(rec.j$Comp.3, 0.75)
IQR <- Q3-Q1
multiplier <- qnorm(0.95)
whisker_low_PC3_year <- max(min(rec.j$Comp.3), Q1 - multiplier * IQR)
whisker_high_PC3_year <- min(max(rec.j$Comp.3), Q3 + multiplier * IQR)
out <- cbind (out, whisker_low_PC3_year, whisker_high_PC3_year)
out <- as.data.frame(out)
out2 <- rbind(out2,out)
}
write.csv(out2, "Stats.PC3_year.csv")

# get rid of 10% of the data and get the the min and max value than - quantiles manually calculated and later merged in Excel 
out <- c()
out2 <- c()
for (j in 1:length(sp)){
rec.j <- subset(data, data$Species == sp[j])
PC1_lower_q <- quantile(rec.j$Comp.1, 0.05)
PC1_upper_q <- quantile(rec.j$Comp.1, 0.95)
PC2_lower_q <- quantile(rec.j$Comp.2, 0.05)
PC2_upper_q <- quantile(rec.j$Comp.2, 0.95)
PC3_lower_q <- quantile(rec.j$Comp.3, 0.05)
PC3_upper_q <- quantile(rec.j$Comp.3, 0.95)
out2 <- as.data.frame(cbind(PC1_lower_q,PC1_upper_q,PC2_lower_q, PC2_upper_q, PC3_lower_q, PC3_upper_q ))
rownames(out2) <- sp[j]
out <- rbind(out,out2)
}
write.csv(out, "Quantiles_year.csv")

## first data from all, av and year are merged and later all PCA_stats 
PC1 <- read.csv("PCAStats_all.csv")
PC2 <- read.csv("PCAStats_av.csv")
PC3 <- read.csv("PCAStats_year.csv")

comb12 <- merge(PC1,PC2, by = "Species")
comb <- merge(comb12, PC3, by = "Species")
write.csv(comb, "All_Stats_Merged.csv")



# merge with other species data
GLM <- read.csv("forGLM.csv")
GLM[,1] <- rm()
PCA <- read.csv("PCAStats.csv")
PCA[,1] <- rm()
comb <- merge(GLM,PCA, by = "Species")
write.csv(comb, "GLM_PCA.csv", row.names = FALSE)

# conduct the GLM Analysis 
library(MASS)
wd <- "D:/.../PCA_Hypervolume/PCA_Stats"
setwd(wd)

GLM <- read.csv("D:/.../PCA_Hypervolume/PCA_Stats/All_Stats_Merged_prep_GLM_onlySignSlopes.csv")
 
glm_size_tot <- glm(formula =  AverageSizeTotal ~ SlopeGesamt + med_PC1_all + med_PC2_all + med_PC3_all + Records, data =GLM, family = gaussian()) # AIC: -41.141
stepAIC(glm_size_tot)
glm_size_tot <- glm(formula =  AverageSizeTotal ~  med_PC2_all + med_PC3_all  , data =GLM, family = gaussian()) # AIC: -47.007
summary(glm_size_tot)
sink("glm_size_tot_all_onlySig.txt")
print(summary(glm_size_tot))
sink()

glm_size_whole <- glm(formula =  AverageSizeMaskwhole ~ SlopeMaskwhole + med_PC1_all + med_PC2_all + med_PC3_all + Records, data =GLM, family = gaussian()) # AIC: -33.449
stepAIC(glm_size_whole)
glm_size_whole <- glm(formula =  AverageSizeMaskwhole ~ SlopeMaskwhole +  med_PC2_all + med_PC3_all , data =GLM, family = gaussian()) # AIC: -38.609
summary(glm_size_whole)
sink("glm_size_whole_all_onlySig.txt")
print(summary(glm_size_whole))
sink()

glm_size_rm <- glm(formula =  AverageSizeMaskrm ~ SlopeMaskrm + med_PC1_all + med_PC2_all + med_PC3_all + Records, data =GLM, family = gaussian()) # AIC: -21.438
stepAIC(glm_size_rm)
glm_size_rm <- glm(formula =  AverageSizeMaskrm ~ med_PC2_all + med_PC3_all , data =GLM, family = gaussian()) # AIC: -26.25
summary(glm_size_rm)
sink("glm_size_rm_all_onlySig.txt")
print(summary(glm_size_rm))
sink()

###

glm_slope_tot <- glm(formula =    SlopeGesamt ~ AverageSizeTotal+ med_PC1_all + med_PC2_all + med_PC3_all + Records, data =GLM, family = gaussian()) # AIC: -245.03
stepAIC(glm_slope_tot)
glm_slope_tot <- glm(formula =  SlopeGesamt ~ AverageSizeTotal+ med_PC1_all + med_PC2_all + med_PC3_all , data =GLM, family = gaussian()) # AIC: -249.4
summary(glm_slope_tot)
sink("glm_slope_tot_all_onlySig.txt")
print(summary(glm_slope_tot))
sink()

glm_slope_whole <- glm(formula =    SlopeMaskwhole ~ AverageSizeMaskwhole + med_PC1_all + med_PC2_all + med_PC3_all + Records, data =GLM, family = gaussian()) # AIC: -252.33
stepAIC(glm_slope_whole)
glm_slope_whole <- glm(formula =  SlopeMaskwhole ~ AverageSizeMaskwhole + med_PC1_all + med_PC2_all + med_PC3_all , data =GLM, family = gaussian()) #  AIC: -256
summary(glm_slope_whole)
sink("glm_slope_whole_all_onlySig.txt")
print(summary(glm_slope_whole))
sink()

glm_slope_rm <- glm(formula =   SlopeMaskrm ~  AverageSizeMaskrm + med_PC1_all + med_PC2_all + med_PC3_all + Records, data =GLM, family = gaussian()) # AIC: -255.64
stepAIC(glm_slope_rm)
glm_slope_rm <- glm(formula =  SlopeMaskrm ~ med_PC1_all + med_PC2_all , data =GLM, family = gaussian()) # AIC: -259.2
summary(glm_slope_rm)
sink("glm_slope_rm_all_onlySig.txt")
print(summary(glm_slope_rm))
sink()