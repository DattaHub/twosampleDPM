generate_nhpp <- function(Lambda,Tmax){

  n=rpois(1,Lambda(Tmax))
  Ft=function(x) Lambda(x)/Lambda(Tmax)
  Ftinv=function(u){
    a=0
    b=Tmax
    for(j in 1:20){
      if(Ft((a+b)/2)<=u){binf=(a+b)/2;bsup=b}
      if(Ft((a+b)/2)>=u){bsup=(a+b)/2;binf=a}
      a=binf
      b=bsup
    }
    return((a+b)/2)
  }
  X0=rep(NA,n)
  for(i in 1:n){
    X0[i]=Ftinv(runif(1))
  }
  X=sort(X0)
  return(X)
}