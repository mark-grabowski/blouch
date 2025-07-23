//Varying intercepts model - estimates different alpha and beta values per regime
//Population-specific thetas, e, es
//Prior to implementing hierarchical priors

data {
  int<lower=1> N;       // Number of populations
  int<lower=1> n_theta;
  vector[N] z_obs;    // Observed trait means
  vector[N] x_obs;    // Environmental variables
  vector<lower=0>[N] z_error; // SD for measurement error in y
  vector<lower=0>[N] x_error; // SD for measurement error in x
  matrix[N,N] M;       // Migration matrix
  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
  vector[2] es_prior;
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)
  array[N] int<lower=1, upper=n_theta> reg_assign;
}

parameters {
  vector[n_theta] alpha;         // Intercept
  vector[n_theta] beta;           // Slope
  //real<lower=0> es;      // Individual-level es
  vector<lower=0>[N] es;      // Population-specific es
  vector<lower=0.00001>[N] sd_theta; // Standard deviations for theta deviations
  vector<lower=0.00001>[N] sd_e;     // Standard deviations for residual deviations
  //real<lower=0.00001> sd_theta; // Standard deviations for theta deviations
  //real<lower=0.00001> sd_e;     // Standard deviations for residual deviations

  cholesky_factor_corr[N] L_Omega_theta; // Cholesky factor of correlation matrix for theta
  cholesky_factor_corr[N] L_Omega_e;     // Cholesky factor of correlation matrix for e
  vector[N] z_true;     // True y values
  vector[N] x_true;     // True x values

}

transformed parameters {
  matrix[N,N] A;       // Adaptation matrix
  matrix[N,N] AM_term; // I + A^-1 M
  vector[N] mu;        // Mean of y
  matrix[N,N] Sigma;   // Variance-covariance matrix of y
  matrix[N, N] Omega_theta = L_Omega_theta * L_Omega_theta'; // Recover correlation matrix for theta
  matrix[N, N] Omega_e = L_Omega_e * L_Omega_e';             // Recover correlation matrix for e
  //matrix[N, N] Omega_theta = diag_matrix(rep_vector(1.0, N)); //Diagonal matrix
  //matrix[N, N] Omega_e   = diag_matrix(rep_vector(1.0, N)); //Diagonal matrix

  // Construct A
  for (i in 1:N) {
    //A[i, i] = -es * (1 + M[i, i]); //Single es
    A[i, i] = -es[i] * (1 + M[i, i]); //Population-specific es
    for (j in 1:N) {
      if (i != j) {
        A[i, j] = 0;
      }
    }
  }
  AM_term = diag_matrix(rep_vector(1, N)) + inverse(A) * M;

  for(i in 1:N){
    mu[i] = alpha[reg_assign[i]] + x_true[i] * beta[reg_assign[i]];
  }

  matrix[N, N] diag_sd_theta = diag_matrix(sd_theta);
  matrix[N, N] diag_sd_e = diag_matrix(sd_e);
  //matrix[N, N] diag_sd_theta = diag_matrix(rep_vector(sd_theta, N));
  //matrix[N, N] diag_sd_e = diag_matrix(rep_vector(sd_e, N));

  Sigma = inverse(AM_term) * (diag_sd_theta * Omega_theta * diag_sd_theta) * (inverse(AM_term)') + (diag_sd_e * Omega_e * diag_sd_e);

  vector[N] mu_eq = inverse(AM_term) * mu;

}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  beta ~ normal(beta_prior[1], beta_prior[2]);
  //alpha ~ normal(0,1);
  //beta ~ normal(0,1);
  //es ~ exponential(4); //  mu=1/scale
  //es ~ exponential(es_prior); //  mu=1/scale
  es ~ normal(es_prior[1], es_prior[2]) T[0,]; //Prior to implementing hierarchical prior on es; T[0,] truncates the normal distribution at zero, ensuring es remains positive

  // Priors for standard deviations (e.g., half-normal or exponential)
  //exp(0+0.2^2) - lognoral mean
  //sd_theta ~ exponential(0.5);
  //sd_e ~ exponential(0.5);
  //sd_theta ~ lognormal(log(0.25), 1);  // Mean around 0.25 on the original scale, large SD on log scale
  //sd_e ~ lognormal(log(0.25), 1);  // Mean around 0.25 on the original scale, large SD on log scale
  sd_theta ~ lognormal(log(0.1), 0.5);  // Mean around 0.1 on the original scale, smaller SD on log scale
  sd_e ~ lognormal(log(0.1), 0.5);  // Mean around 0.1 on the original scale, smaller SD on log scale


  // Prior on the Cholesky factor of the correlation matrix
  L_Omega_theta ~ lkj_corr_cholesky(nu_cor);
  L_Omega_e ~ lkj_corr_cholesky(nu_cor);

  // Measurement error
  x_obs ~ normal(x_true, x_error);
  z_obs ~ normal(z_true, z_error);

  // Likelihood
  //z_true ~ multi_normal(mu, Sigma);
  matrix[N, N] L_Sigma = cholesky_decompose(Sigma);//For speed increase
  z_true ~ multi_normal_cholesky(mu_eq, L_Sigma);

}

generated quantities {
  vector[N] log_lik;
  vector[N] var_z;
  matrix[N, N] Cor_theta = Omega_theta;
  matrix[N, N] Cor_e = Omega_e;

  for (n in 1:N) {
    var_z[n] = Sigma[n, n];
    log_lik[n] = normal_lpdf(z_obs[n] | z_true[n], z_error[n]); //The log_lik calculation then simply evaluates how likely the actual observation z_obs[n] was, given that specific draw's value for z_true[n] and the known measurement noise z_error[n].

  }


}
