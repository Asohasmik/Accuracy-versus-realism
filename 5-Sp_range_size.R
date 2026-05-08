# Title:    Species range size for herpetufauna of Romania and background points
# Author:   Sajad Noori, Laura Karolin Steib and Dennis Rödder
# Date:     10.09.2025



library(sp)
library (sf)
library(rSDM)
library(terra)
library(raster)
library(stringr)
library(rmaxent)
library(ggtext)
library(devtools)
library(tidyverse)
library(adehabitatHR)
library(rangeBuilder)
#install_github('johnbaums/rmaxent')


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

# Dynamic range per year
files <- list.files(paste0(wdir, "/maxent/dynamic_ranges_years/"), pattern = ".tif", full.names = TRUE)

# Study period
years <-1989:2023

################################################################################
# Calculate the size of the species ranges
################################################################################


out <- data.frame() # A dataframe for collecting the results
for (j in 1:length(sp_list)){
  j = 1
  sp_dy_ranges <- files[str_detect(files, sp_list[j])]
  
  for (i in 1:length(sp_dy_ranges)){
    rast.i <- rast (sp_dy_ranges[i])
    
    tmpfile <- tempfile(fileext = ".tif")
    size.i <- cellSize(rast.i, mask=TRUE, unit="km", filename=tmpfile, overwrite=TRUE)
    rast.pres <- mask(size.i, rast.i, maskvalues=0)
    res.i <- round(sum(values(rast.pres), na.rm=TRUE), 2)
    res.i
    out.i <- as.data.frame(cbind(species = sp_list [j], 
                                 year = years[i], 
                                 size = res.i))
     
    out <- rbind(out, out.i) 
  }
  
}

# Normalize the species ranges
out_df <- out %>% 
  group_by(species) %>%
  mutate(norm_size = as.numeric(size) / max(as.numeric(size), na.rm=TRUE)) %>%
  ungroup() %>% 
  mutate(average_size = mean(as.numeric(norm_size)))
out_df

write.csv(out_df, paste0(wdir, "/maxent/Area_size/AreasSize_all_species.csv"))



################################################################################
# get the regression lines of the different species
# for the entire time span and the two different time spans (1989:2001; 2002:2023)
################################################################################
# Periods
years_1 <- years[1:13]
years_2 <- years[14:35]

sp_ranges_df <- read.csv(paste0(wdir, "/maxent/Area_size/AreasSize_all_species.csv"))
head(sp_ranges_df)

# A function to apply regression on the species ranges over study period
get_regression_stats <- function(data) {
  model <- lm(size ~ year, data =data)
  summary_model <- summary(model)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  r_squared <- summary_model$r.squared
  p_value <- summary_model$coefficients[2, 4]
  label_x <- min(data$year)
  label_y <- max(data$size) 
  
  equation <- sprintf("y = %.3f  + %.3f * x", intercept, slope)
  r_squared_text <- sprintf("R2 = %.3f", r_squared)
  
  return(data.frame(
    species = unique(data$species),
    intercept = intercept,
    slope = slope,
    r_squared = r_squared,
    p_value = p_value,
    equation = equation,
    r_squared_text = r_squared_text,
    label_x = label_x,
    label_y = label_y
  ))
}

lm_results <- get_regression_stats(sp_ranges_df)

ggplot(sp_ranges_df, aes(x = year, y = norm_size, color = species)) +
  geom_col(fill = "lightgrey", color = "lightgrey", position = "dodge", alpha = 0.6)+  
  geom_smooth(method = "lm", se = FALSE, color = "#404040")  + 
  geom_text(data = lm_results,
            aes(label = paste(equation, r_squared_text, sep = "\n")),
            x = label_x, y = label_y,
            size = 2.8, show.legend = FALSE, hjust = 0, vjust = 1, color ="black") +
  geom_hline(aes(yintercept = average_size, color = species),
             color = "black", linetype = "dashed", size = 1) +
  labs(title = paste0("Variation in Habitat Sizes of", sp_list[j], " over 2002-2023"),
       x = "Year",
       y = "Normalized Habitat Size") +
  theme_minimal() +
  theme(legend.position = "none",
        strip.text = element_markdown(size = 8),  # Use element_markdown for italic text
        axis.text.x = element_text(angle = 45, hjust = 1),
        #axis.text.y = element_text(face = "italic")
        ) +
  facet_wrap(~ species, ncol = 4)


