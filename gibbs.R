## Steps
##
## Step 0. Initializing: Determine s[i] deterministically
## and plot the kernel density estimate

## Step 1. Update s[i] and theta[i] using (3.13)[Incorrect
## uses (3.12)
## Question: If theta[i]'s are updated here, then what is
## the need to update theta*[i] in step 2?

## Step 2. Update theta*[i] using (3.14)

## Step 3. Update sig (not applicable for poisson - one parameter)

## Step 4. Generate f ~ p(f|...) using a poor man's method

## Step 5. Put all this together and iterate n.iter times
## and plot the estimated f with changing colors

##################################################
## The macro gibbs() implements n.iter steps of the MCMC,
## calling in turn macros for each of the transition probabilities.
## Before the actual for loop, the macro initializes several lists that
## will accumulate the states in each transition,
## and starts a plot of a simple kernel density estimate,
## to which the draes f ~ p(f|y) will be added.

gibbs <- function(y,n.iter=1000,M,a,b,ymax,n.burn=500,n.thin=5){##51
#   printf <- function(...) invisible(print(sprintf(...)))
  cat(b,"\n")
  n <- length(y)
  lam <- init.DPk(y)  ## initializes th[1..n]
#   plot(lam,type="l")
  ## set up data structures to record imputed posterior draws .. l
  xgrid <- seq(from=0,to=ymax,length=ymax+1)
  fgrid <- NULL ## we will record imputed draws of f
  njlist <- NULL ## record sizes of 8 largest clusters
  klist <- NULL
#   lamssum <- rep(0,length(lam)) ## record lams values
  lamssum <- NULL ## record lams values
  ## start with a plot of a kernel density of the data
  #   plot(density(y), xlab="X", ylab="Y", bty="l", type="l",
  #        xlim=c(0,200), ylim=c(0,0.6), main="")
  ### Gibbs sampler
  for (iter in 1:n.iter){##52
    lam.s <- sample.lam(y,lam,M,a,b)  ## 1. [th_i |...]
    svec <- lam.s$svec
    lam <- lam.s$lam
#     cat("sample.lam",anyNA(lam),"\n")
    lam <- sample.lams(y,lam,a,b) ## 3. [ths_j | ...]
#     cat("sample.lams",anyNA(lam),"\n")
    ## update running summaries ###############
    f <- fbar(xgrid,lam,n,M,a,b)
    #     lines(xgrid,f,col=iter,lty=3)
    fgrid <- rbind(fgrid,f)
    nj <- table(lam)    ## counts
    njlist <- rbind(njlist,sort(nj,decr=T)[1:8])
    klist <- c(klist,length(nj))
    if ( iter > n.burn && iter%%n.thin==0)
    {
      lamssum <- cbind(lamssum,lam)
    }
#     nclust <- as.numeric(names(which.max(table(klist))))
#     cat(nclust,"\n")
#     M <- length(nj)
#     cat(length(lamssum),"\n")
#     lines(lamssum/iter)
  }##52
  ## report summaries #######################
#   fbar <- apply(fgrid,2,mean)
#   lines(xgrid,fbar,lwd=3,col=2)
  njbar <- apply(njlist,2,mean,na.rm=T) ## why na.rm=T??
  cat("Average cluster sizes:\n",format(njbar),"\n")
  pk <- table(klist)/length(klist)
  cat("Posterior probs p(k): (row1 = k,row2 =p(k) \n")
  print(pk/sum(pk))
  return(list(fgrid=fgrid,klist=klist,njlist=njlist,lamssum=lamssum,svec = svec))
}##51
