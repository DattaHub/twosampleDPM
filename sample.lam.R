## Step 1. The transition probability generates from
## $p(\lambda_i | \lambda_(-i),y)

sample.lam <- function(y,lam,M,a,b){#11
  ## sample
  ## lam[i] ~ p(lam[i] | lam[-i], y)
  ## returns updated lam vector
  n <- length(y)
  svec <- rep(0,n)
  for (i in 1:n){#12
    ## unique values and counts
    nj <- table(lam[-i])  # counts
    lams <- as.numeric(names(nj)) # unique values
    k <- length(nj)
    ## likelihood
    fj <- dpois(y[i],lambda=lams)  # p(y|lam*_j), j=1,..k
    fj[ is.nan(fj) ] <- 0
    f0 <- dnbinom(y[i],size=a,prob=1/(1+b)) #q0
    pj <- c(fj*nj,f0*M) # p(s[i]=j|...),j=1..k, k+1
#     cat("fj",format(fj),"\n")
#     cat(format(pj),"\n")
    s <- sample(1:(k+1),1,prob=pj)  # sample s[i]
    svec[i] <- s
    if (s==(k+1)){#13
      ## generate new th[i] value
      a.new <- y[i]+a
      b.new <- b+1
      lamb.new <- rgamma(1,shape=a.new, rate=b.new)
      lams <- c(lams,lamb.new)
      lam[i] <- lamb.new  # record new th[i]
    }#13
    else lam[i] <- lams[s]
  }#12
  return(list(lam=lam, svec = svec))
}#11
