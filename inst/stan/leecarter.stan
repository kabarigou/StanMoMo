//Code for Lee-Carter model with Poisson or Negative Binomial error

data {
  int<lower = 1> J;                    // number of age categories
  int<lower = 1> T;                    // number of years
  int d[J*T];                          // vector of deaths
  vector[J* T] e;                      // vector of exposures
  int<lower = 1> Tfor;                  // number of forecast years
  int<lower = 0> Tval;                  // number of validation years
  int dval[J*Tval];                     // vector of deaths for validation
  vector[J* Tval] eval;                 // vector of exposures for validation
  int<lower=0,upper=1> family;          // family = 0 for Poisson, 1 for NB
}
transformed data {
  vector[J * T] offset = log(e);        // log exposures
  vector[J * Tval] offset2 = log(eval);     // log exposures for validation
  int<lower = 1> L;                     // size of prediction vector
  L=J*Tfor;
}
parameters {
  real<lower=0> aux[family > 0];       // neg. binomial inverse dispersion parameter
  vector[J] a;                          // alpha_x
  simplex[J] b;                        // beta_x, strictly positive and sums to 1

  real c;                           // drift term
  vector[T-1] ks;                   // vector of kappa

  real<lower = 0> sigma;            // standard deviation of the random walk
}
transformed parameters {        // This block defines a new vector where the first component is zero, this is required for identifiability of the Lee-Carter model. Otherwise, the chains will not converge.
  vector[T] k;
  real phi = negative_infinity();       // neg. binomial dispersion parameter
  k[1] = 0;
  k[2:T]=ks;
  if (family > 0) phi = inv(aux[1]);
}

model {
  vector[J * T] mu;           //force of mortality
  int pos = 1;
  for (t in 1:T) for (x in 1:J) {
    mu[pos] = offset[pos]+ a[x] + b[x] * k[t];      // Predictor dynamics
    pos += 1;
  }

  target += normal_lpdf(ks[1]|c,sigma);
  target += normal_lpdf(ks[2:(T-1)]|c+ks[1:(T- 2)],sigma);    // Random walk with drift prior
  
  if (family ==0){
    target += poisson_log_lpmf(d |mu);                // Poisson log model
  }
  else {
    target +=neg_binomial_2_log_lpmf (d|mu,phi);      // Negative-Binomial log model
  }
  target += normal_lpdf(a|0,10);              // Prior on alpha_x
  target += normal_lpdf(c|0,sqrt(10));              // Prior on drift
  target += dirichlet_lpdf(b|rep_vector(1, J));              // Prior on beta_x
  target += exponential_lpdf(sigma | 0.1);         // Exponential prior for sigma
  if (family > 0) target += normal_lpdf(aux|0,1)- normal_lcdf(0 | 0, 1);
}

generated quantities {
  vector[Tfor] k_p;
  vector[L] mufor; // predicted death rates
  vector[J*T] log_lik; // log likelihood on data
  vector[J*Tval] log_lik2; // log likelihood on validation data
  int pos = 1;
  int pos2= 1;
  int pos3= 1;
  k_p[1] = c+k[T]+sigma * normal_rng(0,1);
  for (t in 2:Tfor) k_p[t] = c+k_p[t - 1] + sigma * normal_rng(0,1);
  if (family==0){
    for (t in 1:Tfor) for (x in 1:J) {
    mufor[pos] = a[x] + b[x] * k_p[t];
    pos += 1;
  }
  mufor = exp(mufor);
  for (t in 1:T) for (x in 1:J) {
    log_lik[pos2] = poisson_log_lpmf (d[pos2] | offset[pos2]+ a[x] + b[x] * k[t]);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
    log_lik2[pos3] = poisson_log_lpmf(dval[pos3] | offset2[pos3]+ a[x] + b[x] * k_p[t]);
     pos3 += 1;
}
  }
  else if (family > 0){
    for (t in 1:Tfor) for (x in 1:J) {
      if ( fabs(a[x] + b[x] * k_p[t])>15   ){
         mufor[pos] = 0;
    pos += 1;
      } else {
         mufor[pos] = gamma_rng(phi,phi/exp(a[x] + b[x] * k_p[t]));
    pos += 1;
          }
  }
  for (t in 1:T) for (x in 1:J) {
     log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | offset[pos2]+ a[x] + b[x] * k[t],phi);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
     log_lik2[pos3] = neg_binomial_2_log_lpmf (dval[pos3] | offset2[pos3]+ a[x] + b[x] * k_p[t],phi);
     pos3 += 1;
}
  }
    }

  // The generated quantities block generates forecasts but also the matrices of log-likelihood that we need to compute the optimal weights for validation

