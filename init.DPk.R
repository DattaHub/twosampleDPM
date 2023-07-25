## Step 0. Choosing initial $\lambda[1:n]$ wisely.
## Using hierarchical clustering
init.DPk <- function(y){#01
  
  # intial EDA estimate of th[1..n]
  # cluster data, and cut height H=10, to get 10 clusters
  hc <- hclust(dist(y)^2, "cen")
  s <- cutree(hc,k=2) # cluster membership indicators
  lams <- sapply(split(y,s),mean) # cluster specific means
  lam <- lams[s] # return th[i] = ths[s_i])
  
  return(lam)
}#01