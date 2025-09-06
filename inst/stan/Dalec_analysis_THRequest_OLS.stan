//Varying intercepts OLS model - estimates different alpha and beta values per regime
//Single es parameter, single alpha and beta for Greenhouse data


data {
  int<lower=1> N_pops;       // Number of populations - 10 for greenhouse
  vector[N_pops] z_obs;    // Observed trait means - witout missing values - Field is first here, then greenhouse
  vector<lower=0>[N_pops] z_error_obs; // SD for measurement error in z - staring with complete case
  vector[N_pops] x_obs;    // Environmental variables
  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma; // Residual standard deviation
  vector [N_pops] z_true;     //Measurement error corrected
}

transformed parameters {
}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  beta ~ normal(beta_prior[1], beta_prior[2]);

  // Measurement error - also combines z from field and greenhouse into one vector 20 populations long
  z_obs ~ normal(z_true, z_error_obs);

  // Likelihood
  z_true ~ normal(alpha+beta*x_obs, sigma);
  //sigma ~ lognormal(log(0.5), 1);  //
  sigma ~ lognormal(log(1), 0.5);  //

  //z_true_pops ~ normal(mu_pops, sigma_e);
  //target += normal_lpdf(z_true_pops_raw | 0, 1);



}

generated quantities {
}



