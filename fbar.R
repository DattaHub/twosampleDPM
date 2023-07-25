
## Step 4.
## The macro fbar() implements the draw from f ~ p(f|\theta^*,\sigma).

fbar <- function(xgrid,lam,n,M,a,b){#41
  ## conditional draw F ~ p(F| th,sig,y) (approx -- will talk about this)
  ##
  nj <- table(lam)  # counts
  lams <- as.numeric(names(lam)) # unique values
  k <- length(nj)
  fx <- M/(M+n) * dnbinom(xgrid, size = a, prob = 1/(1+b))
  for (j in 1:k){#42
    fx <- fx + nj[j]/(n+M) * dpois(xgrid, lambda = lams[j])
  }#42
  return(fx)
}#41
