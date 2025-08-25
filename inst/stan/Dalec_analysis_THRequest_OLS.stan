//Varying intercepts OLS model - estimates different alpha and beta values per regime
//Single es parameter, single alpha and beta for Greenhouse data


functions{
      vector merge_missing( int[] miss_indexes , vector x_obs , vector x_miss ) {
        int N = dims(x_obs)[1];
        int N_miss = dims(x_miss)[1];
        vector[N] merged;
        merged = x_obs;
        for ( i in 1:N_miss )
            merged[ miss_indexes[i] ] = x_miss[i];
        return merged;
    }
}
data {
  int<lower=1> N_pops;       // Number of populations - 10 for field

  vector[N_pops] z_obs;    // Observed trait means - witout missing values - Field is first here, then greenhouse
  int N_z_miss; //Number of missing trait means
  int z_missidx[N_z_miss]; //Index of which z's are missing from full dataset

  vector[N_pops] x_obs;    // Environmental variables
  int N_x_miss; //Number of missing trait means
  int x_missidx[N_x_miss]; //Index of which x's are missing from full dataset

  vector<lower=0>[N_pops] z_error_obs; // SD for measurement error in z - staring with complete case

  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
}

parameters {
  real alpha;
  real beta;
  real<lower=0> sigma_e; // Residual standard deviation
  vector [N_pops] z_true_pops_raw;     // Non-centered version
}

transformed parameters {
  vector[N_pops] mu_pops = alpha + x_obs * beta;
  vector [N_pops] z_true_pops = mu_pops + z_true_pops_raw * sigma_e;
}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  beta ~ normal(beta_prior[1], beta_prior[2]);

  // Measurement error - also combines z from field and greenhouse into one vector 20 populations long
  z_obs ~ normal(z_true_pops, z_error_obs);

  // Likelihood
  //z_true ~ multi_normal(mu, Sigma);
  sigma_e ~ lognormal(log(0.5), 0.5);  //

  z_true_pops ~ normal(mu_pops, sigma_e);
  target += normal_lpdf(z_true_pops_raw | 0, 1);

}

generated quantities {
  vector[N_pops] log_lik;
  for (n in 1:N_pops) {
    log_lik[n] = normal_lpdf(z_obs[n] | mu_pops[n], sqrt(sigma_e^2 + z_error_obs[n]^2));
  }
}



