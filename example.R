### EXAMPLE: running cluster analysis on Visium sample

## Load libraries and functions
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(hdf5r)
library(raster)

source("ClusterTendency.R") 

## Load data and preprocess

folder_path <- "" 
samp1 <- Load10X_Spatial(data.dir = folder_path, "filtered_feature_bc_matrix.h5")

## Check counts

plot1 <- VlnPlot(samp1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(samp1, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

## Normalize counts

samp1 <- SCTransform(samp1, assay = "Spatial", verbose = FALSE)

## Plot expression of genes of interest

SpatialFeaturePlot(samp1, features = c("PDGFRA")) 

## Run cluster analysis 

clust_tend(Visium_sample = samp1, gene = "PDGFRA", threshold = 1, num_sims = 10000)
