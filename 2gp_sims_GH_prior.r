##################################################################
### R code for two-groups test with Gauss-hypergeometric prior  
### Simualtion example as shown in Datta, Banerjee & Dunson (2023). 
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

obswindow <- seq(0,70, length.out = 50)

npts1 <- rep(0,length(obswindow)-1)
npts2 <- rep(0,length(obswindow)-1)
for (i in 2:length(obswindow)){
  npts1[i] <- sum((X3<=obswindow[i]))-sum((X3<=obswindow[i-1]))
  npts2[i] <- sum((X4<=obswindow[i]))-sum((X4<=obswindow[i-1]))
}
y <- cbind(npts1,npts2)


### Stan codes ----
library(plyr)
# library(rstan)
library(parallel)
# library(rbenchmark)
library(reshape2)
library(here)
# rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

## Use CMDSTAN 

library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE, quiet = TRUE)
library(posterior)
library(bayesplot)
color_scheme_set("brightblue")


library(here)

seed.val = 100
stan.iters = 5000
gh.data = list('n'=50, 'tau2'=0.1, 'alpha'=1/2, 'phi' = 0.1,
               'y'= npts1,
               'x' = npts2)

cauchy.code <- file.path(here::here("twosample-cauchy.stan"))
cauchy.mod <- cmdstan_model(cauchy.code)


cauchy.fit <- cauchy.mod$sample(
  data = gh.data, seed = 123,  chains = 1,
  parallel_chains = 1,  refresh = 10,
  iter_sampling = stan.iters, iter_warmup = stan.iters/2,
  init = 0)

library(posterior)

# draws <- theta.array <- cauchy.fit$draws()
# # Convert to data frame using posterior::as_draws_df
# cauchy.draws.df <- as_draws_df(draws)

# c("theta1","theta2","gamma","lambda")

theta1.array <- cauchy.fit$draws("theta1")
theta1.u.df <- as_draws_df(theta1.array) # as_draws_matrix() for matrix
theta1.u.summary <- summarise_draws(theta1.u.df) ## you can also use custome set of summary like ("mean", "median", Mode, "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")

theta2.array <- cauchy.fit$draws("theta2")
theta2.u.df <- as_draws_df(theta2.array) # as_draws_matrix() for matrix
theta2.u.summary <- summarise_draws(theta2.u.df) ## you can also use custome set of summary like ("mean", "median", Mode, "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")

gamma.array <- cauchy.fit$draws("gamma")
gamma.u.df <- as_draws_df(gamma.array) # as_draws_matrix() for matrix
gamma.u.summary <- summarise_draws(gamma.u.df) ## you can also use custome set of summary like ("mean", "median", Mode, "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")

lambda.array <- cauchy.fit$draws("lambda")
lambda.u.df <- as_draws_df(lambda.array) # as_draws_matrix() for matrix
lambda.u.summary <- summarise_draws(lambda.u.df) ## you can also use custome set of summary like ("mean", "median", Mode, "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")

gh.mean.theta1 = theta1.u.summary$mean
gh.mean.theta2 = theta2.u.summary$mean
gh.mean.gamma = gamma.u.summary$mean
gh.mean.lambda = lambda.u.summary$mean

## Separation Probaility
kappa3 <- gh.mean.lambda*gh.mean.gamma/(1+(1+gh.mean.gamma)*gh.mean.lambda)
fit.kappa <- kmeans(kappa3,centers=2)
plot(kappa3,col=fit.kappa$cluster,pch=17)

## Inclusion Probability
kappa2 <- gh.mean.lambda/(1+gh.mean.lambda)
fit.kappa <- kmeans(kappa2,centers=2)
plot(kappa2,col=fit.kappa$cluster,pch=17)

## Inclusion Probability
kappa1 <- gh.mean.gamma/(1+gh.mean.gamma)
fit.kappa1 <- kmeans(kappa1,centers=2)
plot(kappa1,col=fit.kappa1$cluster,pch=17)

library(ggplot2)
source("C:/Users/jyotishka/OneDrive/Documents/R/Count2/plot_utils.R")


plot.data = rbind(
  data.frame(class="X",type="observations",value =npts1, x=seq(1:100)),
  data.frame(class="Y",type="observations",value =npts2, x=seq(1:100)),
  data.frame(class="Shrinkage",type="P(at least one non-zero)",value = 1-(kappa1)*(1-kappa2)*(1-kappa3), x=seq(1:100)),
  data.frame(class="Shrinkage",type="Posterior mean (group 1)",value = gh.mean.theta1, x=seq(1:100)),
  data.frame(class="Shrinkage",type="differential abundance",value = kappa2*kappa3, x=seq(1:100)))

(var.plot <- ggplot(data=plot.data,aes(x=x,y=value,group=class, color = interaction(type,class)))+
    geom_line(aes(colour = interaction(type,class)), alpha=0.75,size=1,show_guide=FALSE)+
    facet_wrap(~type,ncol=2,scale="free")+
    theme_bw()+
    xlab(expression(theta)) + labs(title = "Two-sample count data with GH prior"))

var.plot <- var.plot + theme(axis.title.y = element_text(size = rel(1.2), angle = 90))+
  theme(axis.title.x = element_text(size = rel(1.2)))
var.plot<- var.plot+ theme(axis.text = element_text(size = rel(1.2))) #+
#theme()
#theme(legend.title=element_text(size=15, face="bold"),legend.text=element_text(size=15))
var.plot <- var.plot+theme(strip.text.x = element_text(size=15, face="bold"),strip.text.y = element_text(size=15, face="bold"))
print(var.plot)

ggsave("inclusion_separation.pdf",var.plot, height = 7, width = 7)

ggsave("inclusion_separation.eps",var.plot, height = 7, width = 7, device = cairo_ps)

