
# Title:    Extract climate values for herpetufauna of Romania and background points
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     08.09.2025

library(sp)
library (sf)
library(rSDM)
library(terra)
library(raster)


# Set working directory
wdir <- "C:/.../data"
setwd(wdir)


################################################################################
# Uploading required data and datasets
################################################################################
wdirclim <- paste0(wdir, "/bio4km")
clim <- list.files (wdirclim, ".tif", recursive=T, full.names=T)


speciesdata <- read.csv("species.records.years.csv", h=T)
head(speciesdata)


################################################################################
# Extracting the climate values for each year
################################################################################

years <- names(table(speciesdata$out))

out <- c()

for (i in 1: length(years)){
  sp.clim <- rast (paste (wdirclim ,"/", years[i] ,"/bio.tif", sep=""))
  sp.year <- subset(speciesdata, speciesdata$out == years[i])
  clim.year <- extract(sp.clim, sp.year[,c("x", "y")], cells=T)
  sp.year <- cbind (sp.year,clim.year)
  out <- rbind (out,sp.year)
}
write.csv(out, "results/species.records.with.clim.data.csv")


################################################################################
# Checking for the correlation among the climate variables
################################################################################

# just load in the variables you want to correlate
# take comma delimited CSV files
corr_matrix=function(file_name, method="pearson"){
  corr_data=read.table(file=file_name,h=T,sep=",")
  res=c()
  res=matrix(nrow=(ncol(corr_data)),ncol=(ncol(corr_data)))
  
  for(i in 1:(ncol(corr_data))){
    for(k in 1:(ncol(corr_data))){
      res[i,k]=round(as.numeric(cor.test(corr_data[,i],corr_data[,k],method=method)$estimate^2),3)
    }
  }
  colnames(res)<-colnames(corr_data)
  rownames(res)<-colnames(corr_data)
  write.table(res,file=paste("results_",method,"_squared.csv",sep=""),sep=",")
  return(res)
}


## apply function to your data (type file name in like in the example below
## pearson rank correlation coefficient is set as default!
corr_matrix(file="random.sample.csv")

## for spearman correlation type methods="spearman"
corr_matrix(file="random.sample.csv", method="spearman")
