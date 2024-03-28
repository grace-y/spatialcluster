# FUNCTION: num_clusters
#
# Description: calculate number of clusters in a Visium slide. 
# Cluster is defined as positive spots adjacent to each other.
#
# @param coords Coordinates of all positive spots
# 
# @return Number of spots in clusters

num_clusters <- function(coords){
  
  d <- pointDistance(coords, lonlat=FALSE)
  
  d[d == 0] <- Inf
  
  r <- apply(d, 1, min)
  
  r2 <- r/sqrt(2)
  
  return(sum(r2 == 1)) 
}
