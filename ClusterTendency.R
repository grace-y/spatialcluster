### Function for evaluating tendency for spots to cluster or disperse

clust_tend <- function(
    Visium_sample = temple2b, 
    gene = "CDKN1A", 
    threshold = 0.5, 
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
  
  hist(clust_distr) 
  
  ## Compare number of clusters observed to distribution
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
  
  result <- list("numpos" = numpos, "numspots" = numspots, "percentile" = percentile(x), "p" = p)
  
  return(result)
}




