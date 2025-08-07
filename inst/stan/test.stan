functions {
  matrix compute_Omega(matrix H, real t) {
    return matrix_exp(-H * t);
  }

  matrix compute_V(matrix H, matrix Sigma, real t) {
    int k = dims(H)[1];
    matrix[k, k] integrand;
    matrix[k, k] V;
    matrix[k, k] H_inv = inverse(H); // assuming H is diagonal or easily invertible
    matrix[k, k] eHt = matrix_exp(-H * t);
    V = eHt * Sigma * eHt';
    return V;
  }

  real pruning_likelihood(
    int N,
    int n_nodes,
    int k,
    int[] parent,
    int[] lchild,
    int[] rchild,
    vector[] traits,
    real[] branch_lengths,
    matrix H,
    matrix Sigma,
    vector theta,
    matrix P,
    matrix[] meas_error
  ) {
    vector[n_nodes] loglik;
    vector[k] mu[n_nodes];
    matrix[k, k] V[n_nodes];

    for (i in 1:N) {
      int node = i;
      real t = branch_lengths[node];
      matrix[k, k] Omega = compute_Omega(H, t);
      matrix[k, k] V_i = compute_V(H, Sigma, t) + meas_error[i] + 1e-6 * diag_matrix(rep_vector(1.0, k));

      mu[node] = Omega * theta;
      V[node] = V_i;

      matrix[k, k] L = cholesky_decompose(V_i);
      vector[k] resid = traits[i] - mu[node];
      vector[k] tmp = mdivide_left_tri_low(L, resid);
      loglik[node] = -0.5 * dot_self(tmp) - sum(log(diagonal(L))) - 0.5 * k * log(2 * pi());
    }

    for (i in (N + 1):n_nodes) {
      int lc = lchild[i];
      int rc = rchild[i];
      real t_l = branch_lengths[lc];
      real t_r = branch_lengths[rc];

      matrix[k, k] Omega_l = compute_Omega(H, t_l);
      matrix[k, k] Omega_r = compute_Omega(H, t_r);

      matrix[k, k] V_l = Omega_l * V[lc] * Omega_l' + compute_V(H, Sigma, t_l);
      matrix[k, k] V_r = Omega_r * V[rc] * Omega_r' + compute_V(H, Sigma, t_r);

      V_l += 1e-6 * diag_matrix(rep_vector(1.0, k));
      V_r += 1e-6 * diag_matrix(rep_vector(1.0, k));

      vector[k] mu_l = Omega_l * mu[lc];
      vector[k] mu_r = Omega_r * mu[rc];

      matrix[k, k] V_i = V_l + V_r;
      matrix[k, k] L = cholesky_decompose(V_i);
      vector[k] resid = mu_l - mu_r;
      vector[k] tmp = mdivide_left_tri_low(L, resid);
      loglik[i] = loglik[lc] + loglik[rc] - 0.5 * dot_self(tmp) - sum(log(diagonal(L))) - 0.5 * k * log(2 * pi());

      V[i] = V_i;
      mu[i] = mdivide_left_spd(V_i, V_l * mu_r + V_r * mu_l);  // precision-weighted mean
    }

    return loglik[n_nodes];
  }
}

data {
  int<lower=1> N;
  int<lower=1> n_nodes;
  int<lower=1> k;
  int<lower=1> parent[n_nodes];
  int<lower=0> lchild[n_nodes];
  int<lower=0> rchild[n_nodes];
  real branch_lengths[n_nodes];
  vector[k] traits[N];
  matrix[k, k] meas_error[N];
}

parameters {
  // Trait optima
  vector[k] theta_raw;
  vector<lower=0>[k] sigma_theta;

  // Sigma: diffusion matrix
  vector<lower=0>[k] sigma_Sigma;
  cholesky_factor_corr[k] L_Omega_Sigma;

  // P: stationary covariance
  vector<lower=0>[k] sigma_P;
  cholesky_factor_corr[k] L_Omega_P;

  // H matrix (fixed or estimated)
  // Uncomment below to estimate H
  // matrix<lower=0>[k, k] H_raw;
}

transformed parameters {
  vector[k] theta = theta_raw .* sigma_theta;

  matrix[k, k] Sigma = diag_pre_multiply(sigma_Sigma, L_Omega_Sigma) *
                       diag_pre_multiply(sigma_Sigma, L_Omega_Sigma)';
  matrix[k, k] P = diag_pre_multiply(sigma_P, L_Omega_P) *
                   diag_pre_multiply(sigma_P, L_Omega_P)';

  // Use fixed H here — can be changed if estimating H
  matrix[k, k] H = diag_matrix(rep_vector(2.0, k));  // example: diagonal with 2s
}

model {
  // Priors
  sigma_theta ~ normal(0, 2);
  theta_raw ~ std_normal();

  sigma_Sigma ~ normal(0, 1);
  L_Omega_Sigma ~ lkj_corr_cholesky(2);

  sigma_P ~ normal(0, 2);
  L_Omega_P ~ lkj_corr_cholesky(2);

  // Likelihood
  target += pruning_likelihood(N, n_nodes, k, parent, lchild, rchild,
                               traits, branch_lengths, H, Sigma, theta, P, meas_error);
}

generated quantities {
  real log_lik;
  log_lik = pruning_likelihood(N, n_nodes, k, parent, lchild, rchild,
                               traits, branch_lengths, H, Sigma, theta, P, meas_error);
}
