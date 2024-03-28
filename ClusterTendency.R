# FUNCTION: clust_tend
# 
# Description: evaluating tendency for spots to cluster or disperse
#
# @param Visium_sample Seurat object with a Visium sample dataset
# @param gene Name of gene of interest
# @param threshold Cutoff value for gene of interest, i.e. spots with expression values.
# higher than the threshold will be considered positive. 1 by default
# @param num_sims Number of simulation executions. Use histogram output to determine 
# if more simulations are needed, and decrease if runtime too long. 10000 by default
#
# @return numpos Number of positive spots
# @return numspots Number of total spots on Visium sample
# @return x Number of spots in clusters observed on Visium sample
# @return percentile Percentile of observed number of clusters compared to simulated distribution of clusters
# @return p p-value of observed number of clusters compared to median of simulated distribution of clusters
#
# plot output: Histogram of simulated numbers of clusters

clust_tend <- function(
    Visium_sample, 
    gene, 
    threshold = 1, 
    num_sims = 10000
    ){
  
  # function for evaluating number of clusters
  source("NumClusters.R") 
  
  ## Get all coordinates in Visium sample 
  allcoord <- cbind(Visium_sample@images[["slice1"]]@coordinates[["row"]], Visium_sample@images[["slice1"]]@coordinates[["col"]])
  #plot(allcoord) 

  ## Get number of positive spots
  
  # Pull the raw expression matrix from the original Seurat object containing only the genes of interest
  subset.matrix <- FetchData(Visium_sample, vars = gene) 
  subset.matrix$cells <- rownames(subset.matrix) 
  
  # Identify cells that meet threshold
  subset.matrix <- subset.matrix[subset.matrix[,gene] > threshold,] 
  
  # Subset on cells that meet threshold
  senc <- subset(Visium_sample, cells = subset.matrix$cells) 
  
  #plot(senc@images[["slice1"]]@coordinates[["row"]],senc@images[["slice1"]]@coordinates[["col"]])
  senc_xy <- cbind(senc@images[["slice1"]]@coordinates[["row"]],senc@images[["slice1"]]@coordinates[["col"]])
  
  ## Simulation: get distribution for number of clusters given coordinates and number of positive
  numpos <- nrow(senc_xy)
  numspots <- nrow(allcoord)
  
  clust_distr <- NULL
  
  for (i in 1:num_sims){ # recommend 10000 simulations
    
    randspots <- sample(1:numspots, numpos, replace = FALSE)
    randcoords <- allcoord[randspots,]
    
    #plot(randcoords) 
    #num_clusters(randcoords)
    
    clust_distr <- c(clust_distr, num_clusters(randcoords))
    
  }
  
  hist(clust_distr, main = "Histogram of simulated distribution", xlab = "Number of clusters") 
  
  # Compare number of spots in clusters
  x <- num_clusters(senc_xy)
  
  # Directly calculate percentile of simulated distribution (non-parametric)
  percentile <- ecdf(clust_distr)
  
  # Use conservative estimate of p value
  if(percentile(x) > 0.5){
    larger <- Filter(function(a) a >= x, clust_distr)
    p <- 2*length(larger)/num_sims
  } else{
    smaller <- Filter(function(a) a <= x, clust_distr)  
    p <- 2*length(smaller)/num_sims
  }
  
  if(p > 1) p = 1
  
  result <- list("numpos" = numpos, "numspots" = numspots, "numclusters" = x, "percentile" = percentile(x), "p" = p)
  
  return(result)
}
