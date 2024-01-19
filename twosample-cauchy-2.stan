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
  real<lower=0,upper=1> tau;
  real<lower=0,upper=1> phi;
}

transformed parameters {
  vector[n] lambda;   
  vector[n] gamma;  // Removed the lower=0 constraint
  vector[n] theta1; 
  vector[n] theta2;

  lambda = lambda_step * tau;
  theta1 = theta1_step .* (lambda .* lambda) * tau^2;
  gamma = gamma_step * phi;
  theta2 = theta1 .* (gamma .* gamma);
}

model {
  // Removed the tau ~ uniform(0,1) and phi ~ uniform(0,1) lines
  lambda_step ~ cauchy(0, 1);
  theta1_step ~ gamma(alpha, 1);
  gamma_step ~ cauchy(0, 1);
  x ~ poisson(theta1);
  y ~ poisson(theta2);
}