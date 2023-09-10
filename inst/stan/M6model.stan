//Code for the M6 model with Poisson and Negative Binomial error

data {
  int<lower = 1> J;                    // number of age categories
  int<lower = 1> T;                    // number of years
  array[J*T] int d;                          // vector of deaths
  vector[J* T] e;                      // vector of exposures
  vector[J] age;                        // vector of ages
  int<lower = 1> Tfor;                  // number of forecast years
   int<lower = 0> Tval;                  // number of validation years
  array[J*Tval] int dval;                     // vector of deaths for validation
  vector[J* Tval] eval;                 // vector of exposures for validation
  int<lower=0,upper=1> family;          // family = 0 for Poisson, 1 for NB
}
transformed data {
  vector[J * T] input_offset = log(e);
  vector[J * Tval] offset2 = log(eval);
  int<lower = 1> C;                     // index for the cohort parameter
  int<lower = 1> L;
   C=J+T-1;
  L=J*Tfor;
}
parameters {
  array[family > 0] real<lower=0> aux; // neg. binomial dispersion parameter
  real c1;                           // drift term for kappa_1
  real c2;                           // drift term for kappa_2
  array[3] real<lower = 0> sigma;            // standard deviations for kappa_1 and kappa_2
  vector[T] k;                          // vector of kappa_1
  vector[T] k2;                          // vector of kappa_2
  real<lower=-1,upper=1> rho;

  real psi;                          // Parameters for the AR(2) process
  real psi2;
  vector[C-2] gs;                                     //vector of gamma
}

transformed parameters {
  vector[C] g;
  real phi = negative_infinity();
  if (family > 0) phi =  inv(aux[1]);
  g[1]=0;
  g[2:(C-1)]=gs;
  g[C]=0;
}

model {
  vector[J * T] mu;
  int pos = 1;
  for (t in 1:T) for (x in 1:J) {
    mu[pos] = input_offset[pos]+k[t]+(age[x]-mean(age))*k2[t]+g[t-x+J];      //Predictor dynamics
    pos += 1;
  }

  target += normal_lpdf(k[1]|c1,sqrt(10));
  target += normal_lpdf(k2[1]|c2,sqrt(10));

  target += normal_lpdf(k[2:T] | c1+k[1:(T- 1)], sigma[1]);
  target += normal_lpdf(k2[2:T] | c2+k2[1:(T- 1)] + rho * sigma[2]*(k[2:T]-c1-k[1:(T- 1)])/ sigma[1],
  sigma[2] * sqrt(1 - square(rho)));
  
 target+=normal_lpdf(gs[1]|0,sigma[3]) ;
target+=normal_lpdf(gs[2]|psi*gs[1],sigma[3]) ;
target+=normal_lpdf(gs[3:(C-2)]|psi*gs[2:(C- 3)]+psi2*gs[1:(C- 4)],sigma[3]) ;
  
   if (family ==0){
    target += poisson_log_lpmf(d |mu);                // Poisson log model
  }
  else {
    target +=neg_binomial_2_log_lpmf (d|mu,phi);      // Negative-Binomial log model
  }
  target += exponential_lpdf(sigma| 0.1);
  target += normal_lpdf(c1|0,sqrt(10));
  target += normal_lpdf(c2|0,sqrt(10));
  target += uniform_lpdf(rho|-1,1);
  target += normal_lpdf(psi|0,sqrt(10)); 
  target += normal_lpdf(psi2|0,sqrt(10)); 
  if (family > 0) target += normal_lpdf(aux|0,1)- normal_lcdf(0 | 0, 1); // prior on overdispersion parameter
}

generated quantities {
  vector[Tfor] k_p;
  vector[Tfor] k2_p;
  vector[Tfor] g_p;
  vector[C+Tfor] gf;
  vector[L] mufor;
  vector[J*T] log_lik;
  vector[J*Tval] log_lik2;
  int pos = 1;
  int pos2= 1;
  int pos3= 1;
  
  k_p[1] = c1+k[T]+sigma[1] * normal_rng(0,1);
  for (t in 2:Tfor) k_p[t] = c1+k_p[t - 1] + sigma[1] * normal_rng(0,1);

  k2_p[1] = c2+k2[T]+ rho*sigma[2]/sigma[1]*(k_p[1]-c1-k[T])+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
  
  for (t in 2:Tfor) k2_p[t] = c2+k2_p[t - 1]+ rho*sigma[2]/sigma[1]*(k_p[t]-c1-k_p[t - 1])+sigma[2] * sqrt(1 - square(rho))*normal_rng(0,1);
  
  g_p[1]=psi*g[C]+psi2*g[C-1]+sigma[3]*normal_rng(0,1);
  g_p[2]=psi*g_p[1]+psi2*g[C]+sigma[3]*normal_rng(0,1);
  for (t in 3:Tfor) g_p[t]=psi*g_p[t-1]+psi2*g_p[t-2]+sigma[3]*normal_rng(0,1);
  gf=append_row(g,g_p);

  if (family==0){
    for (t in 1:Tfor) for (x in 1:J) {
    mufor[pos] = k_p[t]+(age[x]-mean(age))*k2_p[t]+gf[T+t-x+J];
    pos += 1;
  }
  mufor=exp(mufor);
  for (t in 1:T) for (x in 1:J) {
    log_lik[pos2] = poisson_log_lpmf (d[pos2] | input_offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t]+g[t-x+J]);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
    log_lik2[pos3] = poisson_log_lpmf(dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t]+gf[T+t-x+J]);
     pos3 += 1;
}
  }
  else if (family > 0){
     for (t in 1:Tfor) for (x in 1:J) {
      if ( abs(k_p[t]+(age[x]-mean(age))*k2_p[t]+gf[T+t-x+J])>15   ){
         mufor[pos] = 0;
    pos += 1;
      } else {
         mufor[pos] = gamma_rng(phi,phi/exp(k_p[t]+(age[x]-mean(age))*k2_p[t]+gf[T+t-x+J]));
    pos += 1;
          }
  }
  for (t in 1:T) for (x in 1:J) {
     log_lik[pos2] = neg_binomial_2_log_lpmf (d[pos2] | input_offset[pos2]+ k[t]+(age[x]-mean(age))*k2[t]+g[t-x+J],phi);
     pos2 += 1;
}

 for (t in 1:Tval) for (x in 1:J) {
     log_lik2[pos3] = neg_binomial_2_log_lpmf (dval[pos3] | offset2[pos3]+ k_p[t]+(age[x]-mean(age))*k2_p[t]+gf[T+t-x+J],phi);
     pos3 += 1;
}
  }
  }


