//Code for Lee-Carter model with Negative-Binomial error

data {
  int<lower = 1> J;                    // number of age categories
  int<lower = 1> T;                    // number of years
  int d[J*T];                          // vector of deaths
  vector[J* T] e;                      // vector of exposures
  int<lower = 1> Tfor;                  // number of forecast years
}
transformed data {
  vector[J * T] offset = log(e);
  int<lower = 1> L;
  L=J*Tfor;
}
parameters {
  real<lower=0> phi; // neg. binomial dispersion parameter
  vector[J] a;                          // alpha_x
  simplex[J] b;                        // beta_x, strictly positive and sums to 1
  
  real alpha;                          // drift in the random walk
  vector[T-1] ks;                   // vector of kappa
  
  real<lower = 0> sigma[1];            // standard deviation of the random walk
}
transformed parameters {        // This block defines a new vector where the first component is zero, this is required for identifiability of the Lee-Carter model. Otherwise, the chains will not converge.
  vector[T] k;
  k[1] = 0;
  k[2:T]=ks;
}

model {
  vector[J * T] mu;
  int pos = 1;
  for (t in 1:T) for (x in 1:J) {
    mu[pos] = offset[pos]+ a[x] + b[x] * k[t];      // Predictor dynamics
    pos += 1;
  }
  
  ks[1] ~ normal(alpha,sigma[1]);
  ks[2:(T-1)] ~ normal(alpha + ks[1:(T- 2)],sigma[1]);    // Normal dynamics for the random walk
  
  target +=neg_binomial_2_log_lpmf (d|mu,phi);      // Negative-Binomial log model
  target += exponential_lpdf(sigma[1] | 1);         // Exponential prior for sigma
}

generated quantities {
  vector[Tfor] k_p;
  vector[L] mufor;
  vector[J*T] log_lik;
  int pos = 1;
  int pos2= 1;
  k_p[1] = k[T]+alpha + sigma[1] * normal_rng(0,1);
  for (t in 2:Tfor) k_p[t] = k_p[t - 1] + alpha + sigma[1] * normal_rng(0,1);
  for (t in 1:Tfor) for (x in 1:J) {
    mufor[pos] = gamma_rng(phi,phi/exp(a[x] + b[x] * k_p[t]));
    pos += 1;
  }
 for (t in 1:T) for (x in 1:J) {
    log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | offset[pos2]+ a[x] + b[x] * k[t],phi);
     pos2 += 1;
}
  }
  
  // The generated quantities block generates forecasts but also the matrices of log-likelihood that we need to compute the optimal weights
