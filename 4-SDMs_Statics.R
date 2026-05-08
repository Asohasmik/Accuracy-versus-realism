# Title:    SDMs for herpetufauna of Romania
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     08.09.2025



library(sp)
library (sf)
library(rSDM)
library(terra)
library(raster)
library(stringr)
library(foreach)
library(rmaxent)
library(devtools)
library(doParallel)
library(adehabitatHR)
library(rangeBuilder)
#install_github('johnbaums/rmaxent')


# Set working directory
wdir <- "D:/.../data"
setwd(wdir)


source("AICc.functions.r")


ncluster = 6
cl<-makeCluster(ncluster)
doParallel::registerDoParallel(cl)


reps<- c(5) # Number of reapiting the script

clim <- list.files(paste(wd,"/climate", sep=""), ".asc", full.names=TRUE) # Environmental variables
clim <- stack(clim)

# proj paths
scenarios <- c(paste(wd,"/proj/current",sep="")) 

bias <- raster("bias_layer.asc")

rec <- read.csv("records.csv", h=TRUE) # Species occurrences
nrec <- table(rec$species)
nrec <- as.data.frame(nrec)
nrec <- subset(nrec, nrec[,2] > 9)

sp <- names(table(rec$species)) # List of the species in the dataset

foreach (s = 1:length(sp), .packages = c("adehabitatHR","terra", "raster", "usdm", "dismo","ENMeval", "rSDM")) %dopar% {
  # create background files
  sp.i <- rec[rec$species==sp[s],]
  
  dir.create(paste(wd,"/results_", sp[s],"/", sep=""))
  clim.sp.s <- raster::extract(clim, rec[,2:3])
  sp.s <- na.omit(cbind(rec, clim.sp.s))
  
  # test scheme
  reg <- seq(0.5, 2.5, 0.25)
  reg <- c(reg, 5, 10)
  
  # create all possible feature combinations
  classesToTest <- c("-p", "-q", "-h", "threshold=true")
  classGrid <- expand.grid(rep(list(c(1, 0)), length(classesToTest)))
  names(classGrid) <- classesToTest
  
  
  feat <- c()
  for (f in 1:nrow(classGrid)){
    feat.f <- colnames(classGrid[f,])[which(classGrid[f,] == 1)]
    feat.f.i <- c()
    for (fi in 1:length(feat.f)){
      feat.f.i <- paste(feat.f.i, feat.f[fi])
    }
    feat <- c(feat, feat.f.i)
  }
  feat <- feat[-length(feat)]
  
  AICc <- c()
  MResults<-c()
  for (b in 1:length(reg)){
    for (j in 1:length(feat)){
      test.i.j <- paste("java -mx2000m -jar maxent.jar -e"," ",  
                        "nowarnings"," ", 
                        "noprefixes"," ", 
                        "nopictures"," ", 
                        "outputdirectory=", paste(wd,"/results_", sp[s],"/", sep="")," ", 
                        "samplesfile=", paste(wd, "/records.csv",sep="")," ", 
                        "environmentallayers=", paste(wd, "/climate", sep="")," ", 
                        "biasfile=", paste(wd, "/bias_layer.asc", sep="")," ",
                        "biastype=3"," ",
                        "randomseed"," ", 
                        "noaskoverwrite"," ", 
                        "nowriteclampgrid"," ", 
                        "nooutputgrids", " ",
                        "randomtestpoints=20"," ", 
                        "replicates=",reps," ", 
                        "replicatetype=bootstrap"," ", 
                        "outputformat=raw"," ",
                        "writebackgroundpredictions"," ",
                        "betamultiplier=", reg[b]," ",
                        feat[j]," ",
                        "-z autorun",
                        sep="")
      
      system(test.i.j)
      MResults.i.j <- read.csv(paste(wd,"/results_", sp[s],"/maxentResults.csv", sep=""), h=TRUE)
      MResults <- rbind(MResults, MResults.i.j)
      write.csv(MResults, paste(sp[s],"_maxentResults.summary.csv", sep=""))
      
      my.lambdas <- list.files(paste(wd,"/results_", sp[s],"/", sep=""), ".lambdas",full.names = TRUE)
      #		my.ascii <- list.files(paste(wd, "/results_", sp[s], sep=""), ".asc", full.names = T)
      my.occ <- list.files(paste(wd,"/results_", sp[s],"/", sep=""), "samplePredictions", full.names = TRUE)
      my.bg <- list.files(paste(wd,"/results_", sp[s],"/", sep=""), "backgroundPredictions", full.names = TRUE)
      
      for (k in 1:length(my.lambdas)){
        x <- readLines(my.lambdas[k])
        nparam <- get.params(x)
        my.occ.k <- read.csv(my.occ[k], h=TRUE)
        my.occ.k <- my.occ.k[my.occ.k$Test.or.train == "train",]
        #				my.ascii.k <- raster(my.ascii[k])
        my.ascii.k <- read.csv(my.bg[k])
        AIC.c.k <- calc.aicc(nparam, my.occ.k$Raw.prediction, my.ascii.k$raw)
        AICc.k <- c(my.bg[k],reg[b], feat[j], nparam, AIC.c.k$AICc)
        AICc <- rbind(AICc, AICc.k)
      }	
      write.csv(AICc, paste(sp[s],"_AICc.csv", sep=""))
      trash <- list.files(paste(wd,"/results_", sp[s],"/", sep=""),full.names = TRUE)
      file.remove(trash)
    }
  }
  colnames(AICc) <- c("Model","Regularization","Features","nParameters","AICc")
  write.csv(AICc, paste(sp[s],"_AICc.csv", sep=""))
}
stopCluster(cl)


############################# be Happy!!!#######################
## Summarize AICc results

for (i in 1:length(sp)){
  if (file.exists(paste(sp[i], "_AICc.csv",sep=""))==TRUE){
    AICc.sp.i <- read.csv(paste(sp[i], "_AICc.csv",sep=""), h=TRUE)
    #reps
    maxentResults.sp.i <- read.csv(paste(sp[i], "_maxentResults.summary.csv",sep=""), h=TRUE)
    
    for (n in 1:nrow(maxentResults.sp.i)) {
      maxentResults.sp.i[n,1] <- paste(wd,"/results_", sp[i],"/", maxentResults.sp.i[n,2], ".asc", sep="")
    }
    maxentResults.sp.i <- maxentResults.sp.i[maxentResults.sp.i[,2] != paste(sp[i]," (average)",sep=""),]
    ntests <- nrow(maxentResults.sp.i)/reps
    starts <- seq(1, (ntests)*reps, reps)
    starts <- starts[-length(starts)]
    res <- c()
    res.median <- c()
    for (j in 1:length(starts)){
      starts.j <- starts[j]
      stops <- starts[j] + reps -1 
      maxentResults.sp.i.j <- maxentResults.sp.i[starts.j:stops,]
      colnames(maxentResults.sp.i.j)[1] <- c("Model")
      AICc.sp.i.j <- AICc.sp.i[starts.j:stops,]
      out <- cbind(AICc.sp.i.j, maxentResults.sp.i.j)
      if (!any(is.na(out$nParameters)) && min(out$nParameters) != 0) {
        out.median <- c(sp[i], out$Regularization[1], as.character(out$Features[1]), median(na.omit(out$nParameters)), median(na.omit(out$AICc)), median(na.omit(out$Training.AUC)), median(na.omit(out$Test.AUC)))
        res.median <- rbind(res.median, out.median)
        res <- rbind(res, out)
      }
    }
    write.csv(res, paste(sp[i], "_maxentResults.AICc.csv", sep=""), row.names=F)
    colnames(res.median) <- c("Species", "Regularization", "Features", "nParameters", "AICc", "Training.AUC", "Test.AUC")
    write.csv(res.median, paste(sp[i], "_median.maxentResults.AICc.csv", sep=""), row.names=F)
  }
}

# Compute final models
#############################
#ncluster = 2
cl<-makeCluster(ncluster)
doParallel::registerDoParallel(cl) 

foreach (s = 1:length(sp), .packages = c("raster", "usdm", "dismo","ENMeval")) %dopar% {
  #for (s in 1:length(sp)){
  AICc.select <- read.csv(paste(wd, "/", sp[s], "_median.maxentResults.AICc.csv", sep=""), h=TRUE)
  AICc.select <- subset(AICc.select, AICc.select$Training.AUC > 0.7)
  if (nrow(AICc.select)>1){
    AICc.select <- AICc.select[which.min(AICc.select$AICc),]
    write.csv(AICc.select, paste(sp[s], "_final.parameters.csv"))
    
    test.i.j <- paste("java -mx2000m -jar maxent.jar -e"," ",  
                      "nowarnings"," ", 
                      "noprefixes"," ", 
                      "nopictures"," ", 
                      "outputdirectory=", paste(wd, "/results_", sp[s], sep="")," ", 
                      "samplesfile=", paste(wd, "/records.csv",sep="")," ", 
                      "environmentallayers=", paste(wd, "/climate", sep="")," ", 
                      "biasfile=", paste(wd, "/bias_layer.asc", sep="")," ",
                      "biastype=3"," ",
                      "projectionlayers=",scenarios," ",
                      "randomseed"," ", 
                      "noaskoverwrite"," ", 
                      "nowriteclampgrid"," ", 
                      "randomtestpoints=20"," ", 
                      "replicates=100"," ", 
                      "replicatetype=bootstrap"," ", 
                      "writebackgroundpredictions"," ",
                      "nooutputgrids", " ",
                      "betamultiplier=", AICc.select$Regularization," ",
                      AICc.select$Features," ",
                      "-z autorun",
                      sep="")
    system(test.i.j)
    
    # MESS
    proj.mess <- unlist(strsplit(scenarios, "[,]"))
    for (m in 1:length(proj.mess)){
      mess.m <- paste("java -cp maxent.jar density.tools.Novel ", paste(wd, "/climate", sep=""), proj.mess[m], " ",  paste(wd, "/results_", sp[s], "/", sp[s], "_mess_",  unlist(strsplit(proj.mess[m], "[/]"))[4], ".asc", sep=""))
      system(mess.m)
    }
  }
}
stopCluster(cl)

final <- list.files(wd, "final.parameters")
f <- c()
for (i in 1:length(final)){
  f.i <- read.csv(final[i], h=TRUE)
  f <- rbind(f, f.i)
}
write.csv(f, "final.settings.csv")

################################################################################
# Uploading required data and datasets
################################################################################
# occurrnce dataset
occ <- read.csv("species.training.csv") [, 1:3] 

occ$x <- as.numeric(occ$x)
occ$y <- as.numeric(occ$y)
head(occ)

# Vector of Romania
Romania <- vect("Romania.shp")

# MaxEnt Results
thresholds <- read.csv(paste0("D:/Romanian_Herpetofauna/static/maxent2/maxentResultsAverage.csv"))
thresholds$Species <- gsub(" \\(average\\)", "", thresholds$Species)
head(thresholds)


sp <- names(table(thresholds$Species))


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
wd <- "D:/Romanian_Herpetofauna/static/maxent2"


# create 0/1 maps
for (j in 1:length(sp)){

    file_path <- paste0(wd, "/", sp[j], "_current_avg.mxe")
    
    # Read the .mxe file safely
    dis.i <- read_mxe(file_path)

      wd.j <- "D:/Romanian_Herpetofauna/static/SDM_statics"
      writeRaster (dis.i, paste0(wd.j, "/", sp[j], "_statics.tif"),overwrite=TRUE)
  }




################################################################################
# Masking the statics raster by habiat and range
################################################################################ 
# Results of SDM
wd <- "D:/Romanian_Herpetofauna/static/SDM_statics" # directory
files <- list.files(paste0(wd), pattern = "_statics.tif", full.names = TRUE)



# Function for masking resulting raster by species habitat and range
mask_builder <- function(sp, raster, occ){
 
  mask.j <- rast(paste0(wdir, "/masks.margrem.incl500/", sp , ".rm2.mask.tif"))
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


# Dynamic binary distribution per year
for (j in 1:length (sp)){
j = 37
  sp_files <- files[str_detect(files, sp[j])]
  threshold.j <- (thresholds[j, "X10.percentile.training.presence.Cloglog.threshold"]) # 10% threshold for species
  
  rast.i <- rast(sp_files)

    # classify the mean raster
    rast.i [rast.i < threshold.j] <- 0
    rast.i [rast.i >= threshold.j] <- 1
    plot(rast.i)
    area.i <- mask_builder(sp[j], rast.i, occ)     
    # plot(area.i)
    wd.i <- "D:/Romanian_Herpetofauna/static_ranges/"
    writeRaster (area.i, paste0(wd.i, sp[j], "_st_range.tif"), overwrite = TRUE)
}
  

