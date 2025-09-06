functions{
}
data {
  int<lower=1> N_pops;       // Number of populations - 10 for field
  vector[N_pops] z_obs;    // Observed trait means - witout missing values - Field is first here, then greenhouse
  vector[N_pops] z_error_obs; // SD for measurement error in z - staring with complete case
  vector[N_pops] x_obs;    // Environmental variables
  matrix[N_pops, N_pops] M;
  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)

}

parameters {
  real alpha;
  real beta;
  //real<lower=0> es;      // Siingle population wide es
  real<lower=0> sd_theta; // Standard deviationsfor theta deviations
  real<lower=0> sd_e;     // Standard deviation for residual deviations
  vector [N_pops] z_true_pops;     // True z values for 10 populations, MigSel version taking account measurement error due
  vector[N_pops] x_true_pops;
  vector<lower=0>[N_pops] x_obs_sd;
}

transformed parameters {
  real<lower=0> es = 5;      // Siingle population wide es
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

  mu_pops = alpha + x_true_pops * beta; //Assume single alpha and beta across Field and Greenhouse data

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

  es ~ lognormal(log(0.1),0.75); //Ensures positive es

  sd_theta ~ lognormal(log(0.5), 0.2);  // Mean around 0.1 on the original scale, smaller SD on log scale
  sd_e ~ lognormal(log(0.5), 0.2);  // Mean around 0.1 on the original scale, smaller SD on log scale

  // Measurement error
  z_obs ~ normal(z_true_pops, z_error_obs);
  x_obs_sd ~ lognormal(log(1),0.5);
  x_obs ~ normal(x_true_pops, x_obs_sd);


  // Likelihood
  //z_true ~ multi_normal(mu, Sigma);
  z_true_pops ~ multi_normal_cholesky(mu_eq_pops, L_Sigma);

}

generated quantities {
  real half_life = log(2)/es;
}



