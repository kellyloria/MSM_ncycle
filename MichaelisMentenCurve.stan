data {
  int<lower=1> N;                  // Total number of observations
  int<lower=1> D;                  // Number of dates
  int<lower=1, upper=D> date[N];   // Date index for each observation
  vector[N] TMR_NO3;               // Nitrate concentration
  vector[N] Uadd;                  // Uptake values
}

parameters {
  real<lower=0> mu_Vmax;           // Population-level mean for Vmax
  real<lower=0> mu_Km;             // Population-level mean for Km
  real<lower=0> sigma_Vmax;        // Variability of Vmax across dates
  real<lower=0> sigma_Km;          // Variability of Km across dates
  real<lower=0> sigma;             // Residual error
  
  vector<lower=0>[D] Vmax;         // Date-specific Vmax
  vector<lower=0>[D] Km;           // Date-specific Km
}

model {
  // Priors
  mu_Vmax ~ normal(0.001, 10);         // Prior for population-level Vmax
  mu_Km ~ normal(0.001, 10);           // Prior for population-level Km
  sigma_Vmax ~ cauchy(0, 2.5);       // Prior for Vmax variability
  sigma_Km ~ cauchy(0, 2.5);         // Prior for Km variability
  sigma ~ cauchy(0, 2);            // Prior for residual error

  Vmax ~ normal(mu_Vmax, sigma_Vmax); // Date-specific Vmax
  Km ~ normal(mu_Km, sigma_Km);       // Date-specific Km

  // Likelihood
  for (n in 1:N) {
    Uadd[n] ~ normal((Vmax[date[n]] * TMR_NO3[n]) / (Km[date[n]] + TMR_NO3[n]), sigma);
  }
}
