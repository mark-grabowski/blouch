data {
  int<lower=1> N;          // Number of populations
  vector[N] y_obs;        // Observed trait means
  vector[N] x_obs;        // Environmental variables
  vector<lower=0>[N] y_error; // SD for measurement error in y
  vector<lower=0>[N] x_error; // SD for measurement error in x
  matrix[N,N] M;          // Migration matrix
  vector[2] optima_prior; // Prior for optima (mean, sd)
  vector[2] beta_prior;   // Prior for beta (mean, sd)
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)
  int<lower=1> max_n_opt; // Maximum possible number of optima (fixed in data)
  vector<lower=0>[max_n_opt] alpha_prior; // Concentration parameters for Dirichlet prior on prob_opt
}

parameters {
  simplex[max_n_opt] prob_opt[N];        // Probability of each population belonging to each of the max_n_opt optima
  vector[max_n_opt] raw_optima;       // Raw, un-ordered optima locations (up to max_n_opt)
  vector[max_n_opt] raw_beta;         // Raw, un-ordered slopes (up to max_n_opt)
  vector<lower=0>[N] es_raw;
  real<lower=0> mu_es;
  real<lower=0> sigma_es;
  vector<lower=0.00001>[N] sd_theta;
  vector<lower=0.00001>[N] sd_e;
  cholesky_factor_corr[N] L_Omega_theta;
  cholesky_factor_corr[N] L_Omega_e;
  vector[N] y_true;
  vector[N] x_true;
}

transformed parameters {
  vector[max_n_opt] optima = sort_asc(raw_optima); // Order all potential optima
  vector[max_n_opt] beta = raw_beta;           // All potential slopes (order corresponds to sorted optima)
  matrix[N,N] A;
  matrix[N,N] AM_term;
  vector[N] mu;
  matrix[N,N] Sigma;
  matrix[N, N] Omega_theta = L_Omega_theta * L_Omega_theta';
  matrix[N, N] Omega_e = L_Omega_e * L_Omega_e';
  vector<lower=0>[N] es = mu_es + sigma_es * es_raw;

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
    mu[i] = sum(prob_opt[i] .* (optima + x_true[i] * beta)); // Mixture of all potential optima
  }

  matrix[N, N] diag_sd_theta = diag_matrix(sd_theta);
  matrix[N, N] diag_sd_e = diag_matrix(sd_e);

  Sigma = diag_sd_theta * Omega_theta * diag_sd_theta + AM_term * (diag_sd_e * Omega_e * diag_sd_e) * AM_term';
}

model {
  // Priors
  raw_optima ~ normal(optima_prior[1], optima_prior[2]);
  raw_beta ~ normal(beta_prior[1], beta_prior[2]);
  es_raw ~ exponential(1);
  mu_es ~ exponential(4);
  sigma_es ~ lognormal(0, 1);
  sd_theta ~ lognormal(log(0.25), 1);
  sd_e ~ lognormal(log(0.25), 1);
  L_Omega_theta ~ lkj_corr_cholesky(nu_cor);
  L_Omega_e ~ lkj_corr_cholesky(nu_cor);

  // Prior on the assignment probabilities
  for (i in 1:N) {
    prob_opt[i] ~ dirichlet(alpha_prior);
  }

  // Likelihood for assignment
  for (i in 1:N) {
    y_obs[i] ~ normal(mu[i], sqrt(Sigma[i, i])); // Likelihood based on mixture
  }

  x_obs ~ normal(x_true, x_error);
  y_obs ~ normal(y_true, y_error); // Measurement error (potentially redundant)
}

generated quantities {
  array[N] int<lower=1, upper=max_n_opt> estimated_opt_assign;
  for (i in 1:N) {
    estimated_opt_assign[i] = categorical_rng(prob_opt[i]);
  }
  vector[N] log_lik;
  vector[N] var_y;
  matrix[N, N] Cor_theta = Omega_theta;
  matrix[N, N] Cor_e = Omega_e;

  for (n in 1:N) {
    var_y[n] = Sigma[n, n];
    log_lik[n] = normal_lpdf(y_obs[n] | mu[n], sqrt(var_y[n]));
  }
}
