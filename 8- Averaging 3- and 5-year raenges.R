
# Title:    3- and 5-years average SDMs 
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     25.03.2026

library(sp)
library (sf)
library(rSDM)
library(MASS)
library(dismo)
library(terra)
library(raster)
library(magrittr)
library(tidyverse)


# Set working directory
wdir <- "C:/.../data"
setwd(wdir)


################################################################################
# Uploading required data and datasets
################################################################################
# species occurrences data
occs <- read.csv("species.records.years.csv", h=T) %>% 
  select(species = NEW_Names, 
         x, y, year = out) %>% 
  filter(!is.na(species))
head(occs)

################################################################################
# Generate average over 3- and 5-year across the study time
################################################################################

# Read files and extract metadata from the dynamic SDMs
wd <- "C:/.../data/results/sp_year_rasters"

files <- list.files(wd, pattern = "withoutNA.*\\.tif$", full.names = TRUE)

# Extract species name
species <- gsub(".*/|_\\d{4}withoutNA\\.tif", "", files)

# Extract year
years <- as.numeric(gsub(".*_(\\d{4})withoutNA\\.tif", "\\1", files))

df <- data.frame(file = files, species = species, year = years)
head(df)


# Split by species
species_list <- split(df, df$species)


# Function to process ONE species
process_species <- function(df_species, window = 3, out_dir = wd) {
  
  # Sort by year
  df_species <- df_species[order(df_species$year), ]
  
  files <- df_species$file
  years <- df_species$year
  sp <- unique(df_species$species)
  
  # Create time groups (robust to missing years)
  group_id <- floor((years - min(years)) / window)
  groups <- split(files, group_id)
  
  # Loop through groups
  for(i in seq_along(groups)) {
    
    file_group <- groups[[i]]
    
    r <- rast(file_group)
    r_mean <- mean(r, na.rm = TRUE)
    
    # Get years for naming
    yrs <- years[group_id == as.numeric(names(groups)[i])]
    
    start_year <- min(yrs)
    end_year <- max(yrs)
    
    fname <- paste0(
      out_dir, "/", sp,
      "_mean_", window, "yr_",
      start_year, "_", end_year, ".tif"
    )
    
    writeRaster(r_mean, fname, overwrite = TRUE)
  }
}

# Average across 3 years
out_3 <- "C:/.../data/results/3_ranges_years"
dir.create(out_3)
lapply(species_list, process_species, window = 3, out_dir = out_3)

# Average across 5 years
out_5 <- "C:/.../data/results/5_ranges_years"
dir.create(out_5)
lapply(species_list, process_species, window = 5, out_dir = out_5)


library(terra)
library(stringr)


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
  j = 1
  
  threshold.j <- (thresholds[j, "X10.percentile.training.presence.Cloglog.threshold"]) # 10% threshold for species
  sp_files <- files[str_detect(files, sp[j])]
  
  # Read all as SpatRaster objects
  rasters <- lapply(sp_files, rast)
  rstack <- rast(rasters)
  
  # Merge them spatially (weighted mean)
  
  # Define window size (3 or 5)

  window = 3
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







