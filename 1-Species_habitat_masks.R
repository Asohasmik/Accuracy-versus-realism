
# Title:    Creat habitat maks for herpetufauna of Romania
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     03.09.2025

library(sp)
library (sf)
library(rSDM)
library(terra)
library(raster)


# Set working directory
wdir <- "C:/..../data"
setwd(wdir)


################################################################################
# Uploading required data and datasets
################################################################################

# Assigned habitat types for species 
list <- read.csv("Species.inc.Habitats.Updated.csv", h=T)


# raster file for IUCN habitat type of Romania
hab <- rast("iucn.habitats.tif") 


################################################################################
# 1- Generating a make for each species using IUCN suitable habitats
################################################################################

sp <- names(table(list$sp.i.name))
dir.create(paste0(wdir, "/masks"))
setwd(paste0(wdir, "/masks"))


for (j in 1:length(sp)){                            #for each species from sp (list of all species we have)
  sp.i <- subset(list, list$sp.i.name == sp[j])   #create a subset (sp.i) from the liste (named sp.i.name) that are the same as the names in sp for each species [j]
  scores <- unique(hab)                             #loads the IUCN categories present in hab (iucn.habitat.tif) of Romania
  hab.i <- sp.i$NewCode                             #all habitats suitable in Romania for a certain species
  scores <- cbind(scores, 0)                        #create scores which represent if a habitat is suitable for a certain species (1) or not suitable (0)
  for (i in 1:nrow(scores)){
    if (is.element(scores[i,1], hab.i)) {
      scores[i,2] <- 1
    }
  }
  rcl <- scores                                     #create new dataset with the same data as in score
  sp.i.rcl <- classify(hab, scores)                 #create a raster with all Habitats in Romania that are suitable for a ceratain species; classifies the groups from hab(habitats in Romania) to the scores of the species (habitat suitable 1, unsuitable 0)
  sp.i.rcl[sp.i.rcl == 0] <- NA
  plot(sp.i.rcl)                                    #creates the raster with suitable and not suitable areas for a species
  writeRaster(sp.i.rcl, paste0(sp[j],".mask.tif"), NAflag = -9999, overwrite=T)   #writes rasters for each species and saves it in the foulder creates ar the beginning
}


################################################################################
# 2- check if the occurance points are in suitable habitat types and in which habiat type they occur
################################################################################

# Occurrence dataset
sp <- read.csv("species.records.years.csv", h=T)
head(sp)


sp.names <- names(table(sp$NEW_Names))


out <- c()

for (i in 1:length(sp.names)){
  r.i <- rast(paste0(wdir, "/masks/",sp.names[i], ".mask.tif"))
  sp.i <- subset(sp, sp$NEW_Names == sp.names[i])
  r.i <- c(r.i, hab)
  hab.sp.i <- extract(r.i, sp.i[,1:2], cells=T)     # cell= T <- also extract grid cell of the coordinate
  sp.i <- cbind(sp.i, hab.sp.i)
  out <- rbind(out, sp.i)
}

write.csv(out, "habitat.check.csv")

# Remove the NA values
out2 <- na.omit(out)
write.csv(out2, "habitat.check.reduced.csv")

# Species without habitat
out3 <- is.na(out)
write.csv(out3, "habitat.check.reduced.notinhabita.csv")
################################################################################
# 3- move records out of suitable habitat type into suitable habitat
################################################################################

# to use
sp.dat <- read.csv("habitat.check.reduced.csv",h=T)
dim(sp.dat)
sp <- names(table(sp.dat$NEW_Names))
setwd(paste0(wdir, "/masks"))

out4 <- c()
for (i in 1:length(sp)){
temp <- raster(paste0(sp[i],".mask.tif")) # habitat file, same projection as records
rec <- subset(sp.dat, sp.dat$NEW_Names == sp[i])
rec.na <- extract(temp, rec[,2:3])
rec.na <- cbind(rec, rec.na)
# colnames(rec.na)[14] <- "inhab"
out4 <- rbind(out4, rec.na)
}

setwd(wdir)
write.csv(out4, "habitat.check.reduced.inhabitat.csv")
head(out4)


sp.daoutsp.dat <- read.csv("habitat.check.reduced.notinhabita.csv",h=T)
sp <- names(table(sp.dat$NEW_Names))
sp <- sp[-29] # remove Rana dalmatina -> only one record
n.sp <- as.data.frame(table(sp.dat$NEW_Names))
setwd(paste0(wdir, "/masks"))
out5 <- c()

for (i in 1:length(sp)){
  temp <- raster(paste0(sp[i],".mask.tif")) # habitat file, same projection as records
  rec <- subset(sp.dat, sp.dat$NEW_Names == sp[i])
  rec <- SpatialPoints(rec[,2:3], proj4string=crs(temp))
  # convert to terra objects (required by points2nearestcell)
  temp_terra <- rast(temp)
  rec_vect   <- vect(rec)
  check.coords <- capture.output (points2nearestcell(rec, temp), file = "test.txt")
  tab <- read.table("test.txt", nrows = length(rec))
  rec.dat <- subset(sp.dat, sp.dat$NEW_Names == sp[i])
  tab <- cbind(tab, rec.dat)
  out5 <- rbind(out5, tab)
}
setwd(wdir)
write.csv(wdir, "records.moved.csv")



################################################################################
# 4- remove records that are more far away than 200 m  and remove records to close to each other
################################################################################

sp.dat <- read.csv("species.records.years.all.in.habitats.cellnr.csv",h=T)
sp <- names(table(sp.dat$NEW_Names))
sp.dat <- subset(sp.dat, sp.dat$distances < 200)
hist(sp.dat$distances, 50)
summary1 <- as.data.frame(table(sp.dat$NEW_Names))

r <- rast(paste0(wdir, "/masks/Ablepharus_kitaibelii.mask.tif"))
r <- aggregate(r, 5)
e <- extract(r, sp.dat[,16:17], cells=T)
sp.dat.filt <- cbind(sp.dat, e$cell)
colnames(sp.dat.filt)[20] <- "cell"
write.csv(sp.dat.filt, "species.records.years.all.in.habitats.cellnr.csv")

# dups removed in excel
sp.dat <- read.csv("species.records.years.all.in.habitats.cellnr.500m.csv",h=T)
sp <- names(table(sp.dat$NEW_Names))
hist(sp.dat$distances, 50)
summary2 <- as.data.frame(table(sp.dat$NEW_Names))
summary.12 <- cbind(summary1, summary2)












