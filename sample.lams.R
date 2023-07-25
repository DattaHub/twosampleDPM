## Step 2.
## Next we update $\theta^*_j$.
sample.lams <- function(y,lam,a,b){#21
  ## sample lams[j] ~ p(lams[j]|...)
  ## unique values and counts
  nj <- table(lam) #counts
  lams <- sort(unique(lam))  # unique values
  ## use sort(.) to match table counts
  k <- length(nj)
#   cat(format(length(nj)),format(k),"\n")
  for (j in 1:k){#22
    ## find Sj = { i: s[i]=j} and compute sample average over Sj
    idx <- which(lam == lams[j])
    ysumj <- sum(y[idx])
    ## posterior moments for p(lams[j]|...)
    aj <- ysumj+a
    bj <- nj[j] + b
#     cat("sample.lam",anyNA(nj[j]),"\n")
    lamsj <- rgamma(1,shape=aj,rate=bj)
    #     cat(format(is.nan(bj)),"\n")
    ## record the new ths[j] by replacing all the th[i], i in Sj.
    lam[idx] <- lamsj
  }#22
  return(lam)
}#21
