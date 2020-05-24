//Code for APC model with Negative-Binomial error

data {
 int<lower = 1> J;                    // number of age categories
  int<lower = 1> T;                    // number of years
  int d[J*T];                          // vector of deaths
  vector[J* T] e;                      // vector of exposures
  int<lower = 1> Tfor;                  // number of forecast years
}
transformed data {
  vector[J * T] offset = log(e);
  int<lower = 1> C;                     // index for the cohort parameter
  int<lower = 1> L;                     // index for the projections
  C=J+T-1;
  L=J*Tfor;
}
parameters {
   real<lower=0> phi; // neg. binomial dispersion parameter
  vector[J] a;                         // alpha_x
  
  real alpha;                          // drift in random walk
  vector[T-1] ks;                       // vector of kappa
  
  real<lower=-1,upper=1> psi;                          // Parameters for the AR(2) process
  real<lower=-1,upper=1> psi2;
  vector[C-2] gs;                     //  vector of gamma
  
  real<lower = 0> sigma[2];            // standard deviation for the random walk and AR(2) process
}
transformed parameters {  // This block defines a new vector where the first component of kappa and gamma is zero, this is required for identifiability of the RH model. Otherwise, the chains will not converge.
  vector[T] k;
  vector[C] g;
  k[1] = 0;
  k[2:T]=ks;
  g[1]=0;
  g[2:(C-1)]=gs;
  g[C]=0;
    }

model {
  vector[J * T] mu;
  int pos = 1;
  for (t in 1:T) for (x in 1:J) {
    mu[pos] = offset[pos]+ a[x] + k[t]+g[t-x+J];          // Predictor dynamics
    pos += 1;
  }
  
  ks[1] ~ normal(alpha,sigma[1]);
  ks[2:(T-1)] ~ normal(alpha + ks[1:(T- 2)],sigma[1]);      // Dynamics of the random walk
  gs[1] ~ normal(0,1);
  gs[2] ~ normal((psi+psi2)*gs[1],sigma[2]);
  gs[3:(C-2)] ~ normal((psi+psi2)*gs[2:(C- 3)]-psi*psi2*gs[1:(C- 4)],sigma[2]);   //Dynamics of the AR(2) process
  
   target +=neg_binomial_2_log_lpmf (d|mu,phi); //Negative Binomial model
  target += normal_lpdf(a|-3,2);              // Prior on alpha_x
  target += exponential_lpdf(sigma | 1);        //Prior on sigma
}

generated quantities {
  vector[Tfor] k_p;
  vector[Tfor] g_p;
  vector[C+Tfor] gf;
  vector[L] mufor;
  vector[J*T] log_lik;
  int pos = 1;
  int pos2= 1;
  k_p[1] = k[T]+alpha + sigma[1] * normal_rng(0,1);
  for (t in 2:Tfor) k_p[t] = k_p[t - 1] + alpha + sigma[1] * normal_rng(0,1);
  g_p[1]=(psi+psi2)*g[C]-psi*psi2*g[C-1]+sigma[2]*normal_rng(0,1);
  g_p[2]=(psi+psi2)*g_p[1]-psi*psi2*g[C]+sigma[2]*normal_rng(0,1);
  for (t in 3:Tfor) g_p[t]=(psi+psi2)*g_p[t-1]-psi*psi2*g_p[t-2]+sigma[2]*normal_rng(0,1);
  gf=append_row(g,g_p);
  for (t in 1:Tfor) for (x in 1:J) {
    mufor[pos] = gamma_rng(phi,phi/exp(a[x] + k_p[t]+gf[T+t-x+J]));
    pos += 1;
  }
 for (t in 1:T) for (x in 1:J) {
    log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | offset[pos2]+ a[x] + k[t]+g[t-x+J],phi);
     pos2 += 1;
}
}
