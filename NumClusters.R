## Function for evaluating number of clusters
num_clusters <- function(coords){
  
  d <- pointDistance(coords, lonlat=FALSE)
  
  d[d == 0] <- Inf
  
  r <- apply(d, 1, min)
  
  r2 <- r/sqrt(2)
  
  return(sum(r2 == 1)) 
}
