twosample_gibbs_2<-function(y,s,w,M,a,b){

  nrow(y) -> n
  s.1 = s[,1];s.2=s[,2]
  stopifnot(length(s.1)==length(w))
  #s.2[(w==0)]<-s.1[(w==0)]
  for(i in 1:n)
  {
    y.til.j <- sapply(split(y,s),sum) #
    #   N.j.1 = as.vector(table(s.1))
    #   N.j.2 = as.vector(table(s.2))
    M.j = as.vector(table(s.2))
    W.j <- sapply(split(w,s.2),length)
    W = sum(W.j)
    y.til.minus.j = sapply(split(y[-i,1],s.1[-i]),sum) #
    N.minus.j.1 = as.vector(table(s.2[-i]))
    #N.minus.j.2 = as.vector(table(s.2[-i]))
    M.minus.j =  N.minus.j.1 #M_j- = M_j - N_jx- 
    L.i2 = length(y.til.j)
    #cat("length",length(y.til.minus.j),length(y.til.j),"L.i2",L.i2,"\n")
    fj = log(dnbinom(y[i,2],size=a+y.til.minus.j,prob=(M.minus.j)/(b+M.j)))
    f0 = log(dnbinom(y[i,2],size=a,prob=1/(1+b)))
    fj[ is.nan(fj) ] <- 0
    #browser()
    #cat("Why?",length(fj),length(W.j),"\n")
    pj <- c(exp(fj)*W.j,exp(f0)*M)
    #pj[L.i2+1]= M/(M+W)*prod(gamma(a+y.til.j)/(b+M.j)^(a+y.til.j))
    pj[is.infinite(pj)] = 0
    pj <- pj/sum(pj,na.rm=TRUE)
    # cat("element",i,"prob",pj,"\n")
    # cat(pj)
    if (w[i]==0){s.2[i]<-s.1[i]}
    else{
    # cat("first loop\n")
    s.2[i] = sample(1:length(pj),1,prob=pj)
    }
    s = cbind(s.1,s.2)
    if(s.2[i]!=s.1[i]){w[i]<-1}
    else{
      # cat("second loop \n")
      k = s.2[i]
      w[i] = rbinom(1,1,pj[k]*0.1)
      # cat("sampled w",w[i],"\n")
    }
  }
  return(list(s=s,w=w,pj=pj))
}