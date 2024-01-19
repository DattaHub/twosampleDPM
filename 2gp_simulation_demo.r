##################################################################
### R code for two-groups test with DPM-Poisson-Gamma     
### Simulation example as shown in Datta, Banerjee & Dunson (2023). 
##################################################################

setwd("~/GitHub/twosampleDPM")
library(here)
library(pracma)
Tmax = 250
set.seed(2914)

#### Generate data from non-homogeneous Poisson process 
source("generate_nhpp.R")

### Pre-defined intensity function for NHPP 

lambda3=function(x) 2*exp(-((x-50)/10)^2)+20*exp(-((x-10)/10)^2)
Lambda3=function(t) integrate(f=lambda3,lower=0,upper=t)$value

X3<- generate_nhpp(Lambda3,Tmax)

lambda4=function(x) 20*exp(-((x-50)/10)^2)+2*exp(-((x-10)/10)^2)
Lambda4=function(t) integrate(f=lambda4,lower=0,upper=t)$value

X4<- generate_nhpp(Lambda4,Tmax)

### Define observation windows and calculating number of points 

obswindow <- seq(0,70,by = 2)

npts1 <- rep(0,length(obswindow)-1)
npts2 <- rep(0,length(obswindow)-1)
for (i in 2:length(obswindow)){
  npts1[i] <- sum((X3<=obswindow[i]))-sum((X3<=obswindow[i-1]))
  npts2[i] <- sum((X4<=obswindow[i]))-sum((X4<=obswindow[i-1]))
}
y <- cbind(npts1,npts2)

## Initialize group means by running a k-means
fit <- kmeans(y, centers = 2)

hc.1 <- hclust(dist(y[,1])^2, "cen")
init.s.1 <- cutree(hc.1,k=2)
hc.2 <- hclust(dist(y[,2])^2, "cen")
init.s.2 <- cutree(hc.2,k=2)

w = as.numeric(init.s.1!=init.s.2)
s = cbind(init.s.1,init.s.2)

### Source codes for fitting DPM

source("init.DPk.R");source("sample.lam.R");source("sample.lams.R")
source("fbar.R");source("gibbs.R");source("Mode.R")

source("twosample_gibbs_2.R") ## Gibbs for two-sample tests

printf <- function(...) invisible(print(sprintf(...)))

n.iter = 1000
burn.in = n.iter/2
w.all <- NULL
for (iter in 1:n.iter){
  if(iter%%(n.iter/10)==0){printf("Iteration = %d",iter)}
  fit <- twosample_gibbs_2(y,s,w,M=1,a=0.5,b=0.5)
  w <- fit$w; s <- fit$s
  w.all <- rbind(w.all,w)
}
cbind(y,colMeans(w.all))
post.prob = colMeans(w.all)
post.prob = colMeans(w.all[(burn.in+1):n.iter,])

## Plot the posterior probabilities and means

plot.data = rbind(data.frame(class="Group 1",type="Parameters",value =lambda3(obswindow), x=obswindow),
                  data.frame(class="Group 2",type="Parameters",value =lambda4(obswindow), x=obswindow),
                  data.frame(class="Group 1",type="Observations",value =y[,1], x=obswindow),
                  data.frame(class="Group 2",type="Observations",value =y[,2], x=obswindow),
                  data.frame(class="Shrinkage",type="Posterior Probabilty",value = post.prob, x = obswindow)
                  )

library(ggplot2)
var.plot <- ggplot(data=plot.data,aes(x=x,y=value,group=class,colour=class))+
  geom_line(aes(colour=class),position="identity",alpha=0.8,size=1,show.legend = TRUE)+
  facet_wrap(~type,ncol=1,scale="free")+theme_bw()+
  xlab("index")#+
   
var.plot <- var.plot + theme(axis.title.y = element_text(size = rel(1), angle = 90))+
  theme(axis.title.x = element_blank()) #theme(axis.title.x = element_text(size = rel(1)))
var.plot<- var.plot+ theme(axis.text = element_text(size = rel(1))) +
  theme(legend.title=element_text(size=12, face="bold"),legend.text=element_text(size=12),
        legend.position="bottom")
var.plot <- var.plot+theme(strip.text.x = element_text(size=10, face="bold"),strip.text.y = element_text(size=10, face="bold"))
print(var.plot)
