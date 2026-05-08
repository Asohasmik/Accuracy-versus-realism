# Title:    .... for herpetufauna of Romania and background points
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     08.09.2025

install.packages("adehabitatHR")



library(sp)
library (sf)
#library(rSDM)
library(terra)
library(raster)
#library(stringr)
#library(foreach)
#library(rmaxent)
library(devtools)
library(tidyverse)
library(doParallel)
library(adehabitatHR)
library(rangeBuilder)
#install_github('johnbaums/rmaxent')


# Set working directory
wdir <- "C:/.../data/"
setwd(wdir)



################################################################################
# Uploading required data and datasets
################################################################################
# occurrnce dataset
occ <- read.csv("species.training.csv") [, 1:3] 

occ$x <- as.numeric(occ$x)
occ$y <- as.numeric(occ$y)
head(occ)

# Vector of Romania
Romania <- vect(paste0(wdir, "Romania.shp"))

# MaxEnt Results
thresholds <- read.csv(paste0(wdir, "/MaxEntAllYears/MaxEntResults/maxentResultsAverage.csv"))
head(thresholds)


sp <- names(table(thresholds$Species))
years <- 1989:2023


##calculate final parameter of maxent

# |> # I couldn't find this folder in the datast
final <- list.files(paste0(wdir, "final.parameters"))
head(final)
################################################################################
# Generating the final MaxEnt results
################################################################################ 

f <- c()
for (i in 1:length(final)){
  f.i <- read.csv(final[i], h=T)
  f <- rbind(f, f.i)
}
write.csv(f, "final.settings.csv")

# create one file with all MaxEnt Results

f <- list.files(getwd(), "maxentResults.csv", recursive=T, full.names=T)
out <- c()


for (j in 1:length(f)){
  csv.j <- read.csv(f[j])
  csv.j <- (csv.j[26:26,])
  out <- rbind(out, csv.j)
}
write.csv(out, "maxentResultsAverage.csv")

################################################################################
# Getting distribution range for each species/year
################################################################################ 
# MaxEnt results folder
wd <- "D:/Romanian_Herpetofauna/MaxEntAllYears/MaxEntResults"


error_df <- data.frame()

# create 0/1 maps
for (j in 1:length(sp)){
  for (i in 1:length(years)){
    file_path <- paste0(wd, "/", sp[j], "_", years[i], "_avg.mxe")
    
    # Skip missing files
    if (!file.exists(file_path)) {
      error_df <- rbind(error_df, data.frame(
        species = sp[j], year = years[i], stringsAsFactors = FALSE
      ))
      next
    }
    
    # Read the .mxe file safely
    dis.i <- tryCatch(
      read_mxe(file_path),
      error = function(e) NULL
    )
    
    # Skip if reading failed
    if (is.null(dis.i)) next
    
      wd.j <- "D:/Romanian_Herpetofauna/sp_year_rasters"
      writeRaster (dis.i, paste0(wd.j, "/", sp[j], "_", years[i], "withoutNA.tif"),overwrite=TRUE)
  }
}



################################################################################
# overlaying the suitable habitats to create a map over all years
################################################################################ 
# Results of SDM


files <- list.files(paste0(wdir, "/results/sp_year_rasters"), 
                    pattern = "withoutNA.tif", full.names = TRUE)



# Function for masking resulting raster by species habitat and range
mask_builder <- function(sp, raster, occ){
  
  mask.j <- rast(paste0(wdir, "/masks.margrem.incl500/", sp , ".rm2.mask.tif", sep=""))
  com.j <- raster * mask.j
  
  occ.j <- subset(occ, occ$Species == sp)
  sp.j <- occ.j[, c("Species", "x", "y")]
  coordinates(sp.j) <- c("x", "y")
  proj4string(sp.j) <- CRS( "EPSG:4326" )
  
  sp.vec <- getDynamicAlphaHull(occ.j, 
                                fraction = 1, 
                                buff = 1000, 
                                coordHeaders = c("x", "y"),
                                clipToCoast = "terrestrial")
  sp.vec.v <- vect(sp.vec[[1]])
  area.j <- mask(raster, sp.vec.v) 
}



for(j in 1:length(sp)){
  j = 39

  threshold.j <- (thresholds[j, "X10.percentile.training.presence.Cloglog.threshold"]) # 10% threshold for species
  sp_files <- files[str_detect(files, sp[j])]
  
  # Read all as SpatRaster objects
  rasters <- lapply(sp_files, rast)
  rstack <- rast(rasters)

  # Merge them spatially

  
  # weighted mean (if you want to weight recent years more)
  # create weights vector length = nlyr(r)
  n <- nlyr(rstack)
  weights <- seq(1, n) / sum(seq(1, n))  # example: linearly increasing weight (older -> lower)
  w_mean <- app(rstack, fun = function(v) sum(v * weights, na.rm = TRUE) / sum(weights))
  #plot(w_mean)

  # classify the mean raster
  w_mean [w_mean < threshold.j] <- 0
  w_mean [w_mean >= threshold.j] <- 1
  #plot(w_mean)
  
  # dir <- "D:/Romanian_Herpetofauna/dynamic_ranges/"
  # writeRaster (rmean, paste0(dir, sp[j], "_dy_range.tif"), overwrite = TRUE)
  
  
  # Mask by species range
  area <- mask_builder(sp[j], w_mean, occ)   
  #plot (area)
  
  
  writeRaster (area, paste0(wdir, "/results/dynamic_ranges_mask/", sp[j], "_dy_masked.tif"), overwrite = TRUE)
}
  
  





