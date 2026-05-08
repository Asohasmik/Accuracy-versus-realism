# Accuracy-versus-realism
This repository provides R code for comparing Species Distribution Models (SDMs) based on two different climate modeling approaches: (1) static models using averaged environmental variables and (2) dynamic models incorporating annual climatic variability in species habitat suitability predictions.

Repository Contents
1. Species Habitat Masks

Scripts for generating habitat masks and defining species-specific accessible areas based on IUCN habitat masks (Jung et al. 2020), occurrence records, and spatial boundaries.

2. Bias Layers

Tools for creating sampling bias layers and generating background pints  to reduce spatial sampling bias in species distribution models and improve model reliability.

3. Climate Data

Workflows for extracting environmnetal variables for occurrences and background points across environmental raster datasets for ecological analyses and SDMs.

4. SDMs Statistics

Scripts for evaluating species distribution models using multiple statistical metrics, model performance assessments, threshold selection, and predictive analyses.

5. Species Range Size

Methods for estimating and comparing species range sizes using occurrence data and predicted suitable habitats.

6. Gap Analysis

Spatial conservation gap analysis workflows for assessing the overlap between species distributions and protected areas to identify conservation priorities.

7. PCA Analysis

Principal Component Analysis (PCA) workflows for exploring environmental gradients, climatic differentiation, and ecological patterns among species occurrences and background data.

8. Averaging 3- and 5-Year Ranges

Scripts for generating temporal averages of species ranges and environmental conditions across 3-year and 5-year periods to reduce temporal variability.

9. Hypervolume Analysis

Methods for estimating ecological niches and multidimensional environmental spaces using hypervolume approaches for niche comparison and characterization.

10. AICc.function
Functions code for assessing the models AICc for comparison to find the optimal models in SDMs.
