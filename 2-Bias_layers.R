
# Title:    create bias grids over all years
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     03.09.2025

library(sp)
library (sf)
library(rSDM)
library(MASS)
library(dismo)
library(terra)
library(raster)
library(magrittr)


# Set working directory
wdir <- "C:/.../data"
setwd(wdir)


################################################################################
# Uploading required data and datasets
################################################################################
occs <- read.csv("species.records.years.csv", h=T)
head(occs)
climdat <- rast("iucn.habitats.tif")

Romania <- vect("Romania.shp") # shapefile for Romanian border


################################################################################
# 1- Generate a bias layer based on the species occurrence for Romanian herpetofauna
################################################################################

locations <- as.matrix(locations[,1:2])
occur.ras <- rasterize(locations, climdat,   1)
plot(occur.ras)

occur.ras <- raster(occur.ras)
presences <- which(values(occur.ras) == 1)
pres.locs <- coordinates(occur.ras)[presences, ]

dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)))
dens.ras <- raster(dens)
plot(dens.ras)
dens.ras.res <- resample(dens.ras, climdat)
dens.crop <- crop(dens.ras.res, v, mask=T)
writeRaster(dens.crop, "bias.grid.tif", overwrite=T)

dens <- rast("bias.grid.tif")
plot(dens)

################################################################################
# 2- Generate a bias layer for each year
################################################################################

r <- raster(paste0(wdir, "/masks/Ablepharus_kitaibelii.mask.tif"))
plot(r)
sp.dat <- read.csv("species.records.years.all.in.habitats.cellnr.500m.csv",h=T)
year <- names(table(sp.dat$out))
year <- year[12:33]

for (i in 1:length(year)){
  i = 1
  sp.dat.i <- subset(sp.dat, sp.dat$out == year[i])
  occur.ras <- rasterize(as.matrix(sp.dat.i[,16:17]), r, 1)
  presences <- which(values(occur.ras) == 1)
  pres.locs <- coordinates(occur.ras)[presences, ]
  
  dens <- kde2d(pres.locs[,1], pres.locs[,2], n = c(nrow(occur.ras), ncol(occur.ras)))
  dens.ras <- raster(dens)
  plot(dens.ras)
  dens.ras <- resample(dens.ras, r)
  dens.ras <- rast(dens.ras)
  dens.ras <- mask(dens.ras, Romania)
  plot(dens.ras)
  writeRaster(dens.ras, paste0(wdir, "/Years_bias_layers/bias.file.500m.",year[i],".tif", sep=""))
}


################################################################################
# 3- Generate a set of random samples
################################################################################

sp.dat <- read.csv("species.records.years.all.in.habitats.cellnr.500m.csv", h=T)
head(sp.dat)
years <- c(2002:2023)

out.rec <- c()
out.back <- c()
for (i in 1:length(years)){
  setwd(paste0(wdir, "/bio4km/",years[i]))
  sp.i <- subset(sp.dat, sp.dat$out == years[i])
  r.i <- rast("bio.tif")
  sp.i.bio <- extract(r.i, sp.i[,16:17])
  sp.i.bio <- cbind(sp.i$NEW_Names, sp.i[,16:17], sp.i.bio[,2:ncol(sp.i.bio)])
  colnames(sp.i.bio)[1:3] <- c("Species","x","y")
  out.rec <- rbind(out.rec, sp.i.bio)
  
  samplesize <- nrow(sp.i) * 2 # get in total 10000 backround points
  b.i <- rast(paste("P:/Herpetologie/Laura/bias/bias.file.500m.", years[i],".tif", sep=""))
  b.i[is.na(b.i)] <- 0
  s <- xyFromCell(b.i, sample(1:ncell(b.i), samplesize, prob=b.i[]))
  s.r.i <- extract(r.i, as.matrix(s))
  s.r.i <- cbind("Backgound", s, s.r.i)
  colnames(s.r.i)[1] <- "Species"
  out.back <- rbind(out.back, s.r.i)
  image(b.i)
  points(s)
}

write.csv(out.rec, "species.training.csv", row.names=F)
write.csv(na.omit(out.back), "random.sample.csv", row.names=F)