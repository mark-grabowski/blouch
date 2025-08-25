//Varying intercepts model - estimates different alpha and beta values per regime
//Single es parameter, single alpha and beta for Greenhouse data
//Prior to implementing hierarchical priors


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

  vector[N_pops] z_error_obs; // SD for measurement error in z - staring with complete case

  matrix[N_pops, N_pops] M;

  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
  vector[2] es_prior;
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)

  real<lower=0> es;      // Siingle population wide es

}

parameters {
  real alpha;
  real beta;

  real<lower=0> sd_theta; // Standard deviationsfor theta deviations
  real<lower=0> sd_e;     // Standard deviation for residual deviations

  vector [N_pops] z_true_pops;     // True z values for 10 populations, MigSel version taking account measurement error due
}

transformed parameters {
  matrix[N_pops,N_pops] A;       // Adaptation matrix
  matrix[N_pops,N_pops] AM_term; // I + A^-1 M
  vector[N_pops] mu_pops;        // Mean of y for each population, MigSel Version

  matrix[N_pops,N_pops] Sigma;   // Variance-covariance matrix of y
  matrix[N_pops, N_pops] Omega_theta = diag_matrix(rep_vector(1.0, N_pops)); //Diagonal matrix
  matrix[N_pops, N_pops] Omega_e   = diag_matrix(rep_vector(1.0, N_pops)); //Diagonal matrix

  // Construct A
  for (i in 1:N_pops) {
    A[i, i] = -es * (1 + M[i, i]); //Single es
    for (j in 1:N_pops) {
      if (i != j) {
        A[i, j] = 0;
      }
    }
  }

  AM_term = diag_matrix(rep_vector(1, N_pops)) + inverse(A) * M;

  mu_pops = alpha + x_obs * beta; //Assume single alpha and beta across Field and Greenhouse data

  matrix[N_pops, N_pops] diag_sd_theta = diag_matrix(rep_vector(sd_theta,N_pops));
  matrix[N_pops, N_pops] diag_sd_e = diag_matrix(rep_vector(sd_e,N_pops));

  Sigma = inverse(AM_term) * (diag_sd_theta * Omega_theta * diag_sd_theta) * (inverse(AM_term)') + (diag_sd_e * Omega_e * diag_sd_e);
  matrix[N_pops, N_pops] L_Sigma = cholesky_decompose(Sigma);//For speed increase

  vector[N_pops] mu_eq_pops = inverse(AM_term) * mu_pops;
}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  beta ~ normal(beta_prior[1], beta_prior[2]);

  sd_theta ~ lognormal(log(0.5), 0.2);  // Mean around 0.1 on the original scale, smaller SD on log scale
  sd_e ~ lognormal(log(0.5), 0.2);  // Mean around 0.1 on the original scale, smaller SD on log scale

  // Measurement error
  z_obs ~ normal(z_true_pops, z_error_obs);

  // Likelihood
  //z_true ~ multi_normal(mu, Sigma);
  z_true_pops ~ multi_normal_cholesky(mu_eq_pops, L_Sigma);

}

generated quantities {
  // Keep these as they were:
  vector[N_pops] log_lik;

  for (n in 1:N_pops) {
    log_lik[n] = normal_lpdf(z_obs[n] | mu_eq_pops[n], sqrt(Sigma[n,n] + z_error_obs[n]^2));
  }

  real half_life = log(2)/es;
}



