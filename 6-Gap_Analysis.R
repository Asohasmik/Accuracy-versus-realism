# Title:    Gap analysis for herpetufauna of Romania and background points
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     11.09.2025


library(sp)
library (sf)
library(terra)
library(raster)
library(ggtext)
library(stringr)
library(tidyverse)
library(shapefiles)



# Set working directory
wdir <- "C:/.../data"
setwd(wdir)


################################################################################
# Uploading required data and datasets
################################################################################

# occurrnce dataset
occ <- read.csv("species.training.csv")
head(occ)

sp_list <- sort(unique(occ$Species))

PAs <- vect (paste0(wdir, "/PAs/cleanedandtransformed_no_overlap.shp"))
plot(PAs)



# Dynamic range per year
files <- list.files(paste0("D:/Romanian_Herpetofauna/dynamic_ranges/"), pattern = ".tif", full.names = TRUE)

# Species ranges dataset
sp_ranges_df <- read.csv(paste0(wdir, "/maxent/Area_size/AreasSize_all_species.csv"))
head(sp_ranges_df)

################################################################################
# Gap analysis
################################################################################

# Protected areas of each studied species using average Dynamic SDMs
out <- data.frame() # A dataframe for collecting the results
for (j in 1:length(sp_list)){
  # Dynamic species range
  sp_r.i <- rast(files[str_detect(files, sp_list[j])])
  #plot(sp_r.i)
  tmpfile <- tempfile(fileext = ".tif")
  size.j <- cellSize(sp_r.i, mask=TRUE, unit="km", filename=tmpfile, overwrite=TRUE)
  rast.pres <- mask(size.j, sp_r.i, maskvalues=0)
  res.j <- round(sum(values(rast.pres), na.rm=TRUE), 2)
  res.j
  
  # Maske the species range with PAs
  PAs_mask <- mask(sp_r.i, PAs)
  #plot(PAs_mask)
  tmpfile <- tempfile(fileext = ".tif")
  size.p <- cellSize(PAs_mask, mask=TRUE, unit="km", filename=tmpfile, overwrite=TRUE)
  rast.pres2 <- mask(size.p, sp_r.i, maskvalues=0)
  res.p <- round(sum(values(rast.pres2), na.rm=TRUE), 2)
  res.p
  
  
  out.i <- as.data.frame(cbind(species = sp_list [j], 
                               species_areas_dy = res.j,
                               percent_protected = round(res.p/res.j*100, 2)))
  
  out <- rbind(out, out.i) 
}
write.csv(out, paste0(wdir, "/maxent/Area_size/Protected_sp_ranges.csv"))


# Number of protected species within each PA



# Dynamic species range

# Initialize vector to store species counts per polygon
poly_PAs <- st_read(paste0(wdir, "/PAs/cleanedAndTransformed.shp")) %>% 
  select(ID, 
         name = SitName)
poly_PAs <- st_make_valid(poly_PAs)
poly_PAs <- poly_PAs[!st_is_empty(poly_PAs), ]
head(poly_PAs)
plot(poly_PAs)

all_occ <- data.frame()
ricnhess_df <- data.frame()

for(f in 1:length(sp_list)){
  r <- rast(files[str_detect(files, sp_list[f])])
  agg_r <- aggregate(r, fact=10)
  # plot(agg_r)
  
  vals <- as.data.frame(agg_r, xy=TRUE) %>% 
    filter(mean == 1) %>% 
    mutate(species = sp_list[f]) %>% 
    rename(value = mean)
  
  all_occ <- bind_rows(all_occ, vals)
  
  occ_sf <- st_as_sf(all_occ, coords = c("x","y"), crs = crs(poly_PAs))
  occ_pa <- st_join(occ_sf, poly_PAs, join = st_within)
  head(occ_pa)
  
  richness <- occ_pa %>% 
    group_by( ID, name) %>% 
    summarize(
      richness = n_distinct(species)
    ) %>% 
    select(ID, richness) %>% 
    st_drop_geometry() %>% 
    mutate(species = sp_list[j])
  ricnhess_df <- rbind(richness, ricnhess_df)
}
  
write.csv(ricnhess_df, paste0(wdir, "/maxent/PAs_richness/Richness_in_PAs.csv"))


# Add the richness to the PAs polygons
PAs_richness <- poly_PAs %>% 
  left_join(ricnhess_df, by = "ID")

plot(PAs_richness["richness"])
st_write(PAs_richness, paste0(wdir, "/maxent/PAs_richness/Richness_in_PAs.shp"))
