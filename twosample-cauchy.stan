//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=0> n; 
  real<lower=0> alpha;
  int x[n];
  int y[n];
  }
  parameters {
  vector<lower=0>[n] theta1_step;
  vector<lower=0>[n] lambda_step; 
  vector<lower=0>[n] gamma_step;
  real<lower=0,upper=1>tau;
  real<lower=0,upper=1> phi;
  }
transformed parameters {
  vector[n] lambda;   
  vector<lower=0>[n] gamma;  
  vector[n] theta1; 
  vector[n] theta2;
  lambda = lambda_step * tau;
  theta1 = theta1_step .* (lambda .* lambda) * tau^2;
  gamma = gamma_step * phi;
  for (i in 1:n)
  theta2[i] = theta1[i]*(gamma[i]^2);
  }
model {
  tau ~ uniform(0,1);
  phi ~ uniform(0,1);
  lambda_step ~ cauchy(0, 1);
  theta1_step ~ gamma(alpha, 1);
  gamma_step ~ cauchy(0,1);
  x ~ poisson(theta1);
  y ~ poisson(theta2);
  }  