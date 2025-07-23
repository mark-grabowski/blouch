//Varying alpha model - estimates population alpha along with varying alpha and betas
//Single es value for all populations

data {
  int<lower=1> N;          // Number of populations
  vector[N] z_obs;        // Observed trait means
  vector[N] x_obs;        // Environmental variables
  vector<lower=0>[N] z_error; // SD for measurement error in y
  vector<lower=0>[N] x_error; // SD for measurement error in x
  matrix[N,N] M;          // Migration matrix
  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;   // Prior for beta (mean, sd)
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)
  int<lower=1> max_n_theta; // Maximum possible number of alpha (fixed in data)
  vector<lower=0>[max_n_theta] reg_prior; // Concentration parameters for Dirichlet prior on prob_reg


}

parameters {
  simplex[max_n_theta] prob_reg[N];        // Probability of each population belonging to each of the max_n_theta alpha
  vector[max_n_theta] raw_alpha;       // Raw, un-ordered alpha locations (up to max_n_theta)
  vector[max_n_theta] raw_beta;         // Raw, un-ordered slopes (up to max_n_theta)

  real mu_log_es; //Mean of log(es)
  real<lower=0>sigma_log_es; //Overall SD of log(es) - among-populations
  vector[N] log_es_raw; //Non-centered raw effects
  //real<lower=0> es;      // Individual-level es
  //vector<lower=0>[N] es;      // Population-specific es

  real mu_log_sd_theta; //Mean of log(sd_theta)
  real<lower=0>sigma_log_sd_theta; //Overall SD of log(sd_theta) - among-populations
  vector[N] log_sd_theta_raw; //Non-centered raw effects
  //vector<lower=0.00001>[N] sd_theta; // Standard deviations for theta deviations

  real mu_log_sd_e; //Mean of log(sd_e)
  real<lower=0>sigma_log_sd_e; //Overall SD of log(sd_e) - among-populations
  vector[N] log_sd_e_raw; //Non-centered raw effects

  cholesky_factor_corr[N] L_Omega_theta;
  cholesky_factor_corr[N] L_Omega_e;
  vector[N] z_true;
  vector[N] x_true;
}

transformed parameters {
  vector[max_n_theta] alpha = sort_asc(raw_alpha); // Order all potential alpha
  vector[max_n_theta] beta = raw_beta;           // All potential slopes (order corresponds to sorted alpha)
  matrix[N,N] A;
  matrix[N,N] AM_term;
  vector[N] mu;
  matrix[N,N] Sigma;
  matrix[N, N] Omega_theta = L_Omega_theta * L_Omega_theta';
  matrix[N, N] Omega_e = L_Omega_e * L_Omega_e';

  vector[N] log_es = mu_log_es + log_es_raw * sigma_log_es; //Non-centered scaling
  vector[N] log_sd_theta = mu_log_sd_theta + log_sd_theta_raw * sigma_log_sd_theta; //Non-centered scaling
  vector[N] log_sd_e = mu_log_sd_e + log_sd_e_raw * sigma_log_sd_e; //Non-centered scaling

  vector<lower=0>[N] es;      // Population-specific es
  vector<lower=0>[N] sd_theta; // Population-specific standard deviations for theta deviations
  vector<lower=0>[N] sd_e;     // Population-specific standard deviations for residual deviations

  for (i in 1:N) {
    es[i] = exp(log_es[i]); //Transform back to normal scale for calculation of Sigma below
    sd_theta[i] = exp(log_sd_theta[i]); //Transform back to normal scale for calculation of Sigma below
    sd_e[i] = exp(log_sd_e[i]); //Transform back to normal scale for calculation of Sigma below

  }

  matrix[N, N] diag_sd_theta = diag_matrix(sd_theta);
  matrix[N, N] diag_sd_e = diag_matrix(sd_e);

  for (i in 1:N) {
    A[i, i] = -es[i] * (1 + M[i, i]);
    for (j in 1:N) {
      if (i != j) {
        A[i, j] = 0;
      }
    }
  }

  AM_term = diag_matrix(rep_vector(1, N)) + inverse(A) * M;

  for (i in 1:N) {
    mu[i] = sum(prob_reg[i] .* (alpha + x_true[i] * beta)); // Mixture of all potential alpha
  }

  Sigma = inverse(AM_term) * (diag_sd_theta * Omega_theta * diag_sd_theta) * (inverse(AM_term)') + (diag_sd_e * Omega_e * diag_sd_e);
  vector[N] mu_eq = inverse(AM_term) * mu;

}

model {
  // Priors
  mu_log_es ~ normal(log(0.5), 0.3);   // Prior for the mean of log_sd_theta values. Centered on your previous fixed value, but could be wider.
  sigma_log_es ~ normal(0,0.25) T[0,];      // Prior for the SD among log_es values (e.g., Half-Normal(0, 0.5) or Half-Cauchy could also be used)
  log_es_raw ~ std_normal();          // Each raw effect ~ N(0,1) - essentially Z-score

  mu_log_sd_theta ~ normal(log(0.1), 0.5);   // Prior for the mean of log_sd_theta values. Centered on your previous fixed value, but could be wider.
  sigma_log_sd_theta ~ exponential(2);      // Prior for the SD among log_sd_theta values (e.g., Half-Normal(0, 0.5) or Half-Cauchy could also be used)
  log_sd_theta_raw ~ std_normal();          // Each raw effect ~ N(0,1) - essentially Z-score

  mu_log_sd_e ~ normal(log(0.1), 0.5);   // Prior for the mean of log_sd_theta values. Centered on your previous fixed value, but could be wider.
  sigma_log_sd_e ~ exponential(2);      // Prior for the SD among log_sd_theta values (e.g., Half-Normal(0, 0.5) or Half-Cauchy could also be used)
  log_sd_e_raw ~ std_normal();          // Each raw effect ~ N(0,1) - essentially Z-score

  raw_alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  raw_beta ~ normal(beta_prior[1], beta_prior[2]);

  L_Omega_theta ~ lkj_corr_cholesky(nu_cor);
  L_Omega_e ~ lkj_corr_cholesky(nu_cor);

  // Prior on the assignment probabilities
  for (i in 1:N) {
    prob_reg[i] ~ dirichlet(reg_prior);
  }

  // Measurement error
  x_obs ~ normal(x_true, x_error);
  z_obs ~ normal(z_true, z_error);

  matrix[N, N] L_Sigma = cholesky_decompose(Sigma);//For speed increase
  z_true ~ multi_normal_cholesky(mu_eq, L_Sigma);

}

generated quantities {
  array[N] int<lower=1, upper=max_n_theta> estimated_opt_assign;
  vector[N] log_lik;
  vector[N] var_z;
  matrix[N, N] Cor_theta = Omega_theta;
  matrix[N, N] Cor_e = Omega_e;

  for (i in 1:N) {
    estimated_opt_assign[i] = categorical_rng(prob_reg[i]);
  }
  for (n in 1:N) {
    var_z[n] = Sigma[n, n];
    log_lik[n] = normal_lpdf(z_obs[n] | z_true[n], z_error[n]); //The log_lik calculation then simply evaluates how likely the actual observation z_obs[n] was, given that specific draw's value for z_true[n] and the known measurement noise z_error[n].

  }

}
