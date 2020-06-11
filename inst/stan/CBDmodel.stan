//Code for the CBD model with Poisson error

data {
  int<lower = 1> J;                    // number of age categories
  int<lower = 1> T;                    // number of years
  int d[J*T];                          // vector of deaths
  vector[J* T] e;                      // vector of exposures
  vector[J] age;                        // vector of ages
  int<lower = 1> Tfor;                  // number of forecast years
  int<lower = 0> Tval;                  // number of forecast years
  int dval[J*Tval];
  vector[J* Tval] eval;
  int<lower=0,upper=1> family;
}
transformed data {
  vector[J * T] offset = log(e);
  vector[J * Tval] offset2 = log(eval);
  int<lower = 1> L;
  L=J*Tfor;
}
parameters {
  real<lower=0.1> aux[family > 0]; // neg. binomial dispersion parameter
  real alpha;                          // drift in random walk for kappa_1
  real alpha2;                          // drift in random walk for kappa_2
  real<lower = 0> sigma[2];            // standard deviations for kappa_1 and kappa_2
  vector[T] k;                          // vector of kappa_1
  vector[T] k2;                          // vector of kappa_2
}

transformed parameters {
  real phi = negative_infinity();
  if (family > 0) phi = aux[1];
}

model {
  vector[J * T] mu;
  int pos = 1;
  for (t in 1:T) for (x in 1:J) {
    mu[pos] = offset[pos]+k[t]+(age[x]-mean(age))*k2[t];      //Predictor dynamics
    pos += 1;
  }

  k[1] ~ normal(-3,1);
  k[2:T] ~ normal(alpha + k[1:(T- 1)],sigma[1]);            //Random walk dynamics for kappa_1
  k2[1] ~ normal(0,1);
  k2[2:T] ~ normal(alpha2 + k2[1:(T- 1)],sigma[2]);           //Random walk dynamics for kappa_2
   if (family ==0){
    target += poisson_log_lpmf(d |mu);                // Poisson log model
  }
  else {
    target +=neg_binomial_2_log_lpmf (d|mu,phi);      // Negative-Binomial log model
  }
  target += exponential_lpdf(sigma[2] | 1);
}

generated quantities {
  vector[Tfor] k_p;
  vector[Tfor] k2_p;
  vector[L] mufor;
  vector[J*T] log_lik;
  vector[J*Tval] log_lik2;
  int pos = 1;
  int pos2= 1;
  int pos3= 1;


  k_p[1] = k[T]+alpha + sigma[1] * normal_rng(0,1);
  for (t in 2:Tfor) k_p[t] = k_p[t - 1] + alpha + sigma[1] * normal_rng(0,1);
  k2_p[1] = k2[T]+alpha2 + sigma[2] * normal_rng(0,1);
  for (t in 2:Tfor) k2_p[t] = k2_p[t - 1] + alpha2 + sigma[2] * normal_rng(0,1);

  if (family==0){
    for (t in 1:Tfor) for (x in 1:J) {
    mufor[pos] = k_p[t]+(age[x]-mean(age))*k2_p[t];
    pos += 1;
  }
  mufor=exp(mufor);
  for (t in 1:T) for (x in 1:J) {
    log_lik[pos2] = poisson_log_lpmf (d[pos2] | offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t]);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
    log_lik2[pos3] = poisson_log_lpmf(dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t]);
     pos3 += 1;
}
  }
  else if (family > 0){
    for (t in 1:Tfor) for (x in 1:J) {
     mufor[pos] = gamma_rng(phi,phi/exp(k_p[t]+(age[x]-mean(age))*k2_p[t]));
    pos += 1;
  }
  for (t in 1:T) for (x in 1:J) {
     log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t],phi);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
     log_lik2[pos3] = neg_binomial_2_log_lpmf (dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t],phi);
     pos3 += 1;
}
  }
  }


