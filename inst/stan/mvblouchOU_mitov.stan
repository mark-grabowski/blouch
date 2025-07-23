//mvBlouch mvOU - Mitov approach - fixed predictors model - 2025
//Multivariate OU model accounts for measurement error on X and/or Y variables
//After Bartoszek et al. 2012; Clavel et al. 2015; Mitov et al. 2019

functions {
  matrix calc_mu_y_mat(int N, int n_traits, vector path_nodes, array[] int path_nodes_idx, array[] int regimes, array[] int regimes_idx,
    vector branch_lengths, array[] int branch_lengths_idx, matrix F_mat, matrix OU_mat, vector y_0_ancestral){
      int pos = 1;
      int pos_reg = 1;
      int pos_branch_lengths = 1;
      matrix[n_traits,n_traits] I = diag_matrix(rep_vector(1,n_traits));
      matrix[N,n_traits] mu_y_mat;

      for(i in 1:N){ #For each path back to root, calculate E[Z]
        int n_path_nodes_idx = path_nodes_idx[i];
        vector[n_path_nodes_idx] path = segment(path_nodes, pos, path_nodes_idx[i]);
        pos = pos + path_nodes_idx[i];

        int n_regimes_idx = regimes_idx[i];
        array[regimes_idx[i]] int regime_path = segment(regimes, pos_reg, regimes_idx[i]);
        pos_reg = pos_reg + regimes_idx[i];

        int n_branch_lengths_idx = branch_lengths_idx[i];
        vector[branch_lengths_idx[i]] branch_length_path = segment(branch_lengths, pos_branch_lengths, branch_lengths_idx[i]);
        pos_branch_lengths = pos_branch_lengths + branch_lengths_idx[i];

        int n_seg = size(branch_length_path);

        vector[n_traits] current_expected_state = y_0_ancestral;

        for(j in 1:n_seg){ //For each segment of each branch, work form tip to root

          vector[n_traits] term1 = exp(-F_mat*branch_length_path[j]) * current_expected_state;
          matrix[n_traits,n_traits] term2 = I-exp(-F_mat*branch_length_path[j]);
          vector[n_traits] term3 = OU_mat[regimes[j],]';
          current_expected_state = term1 + term2 * term3;
          }
        mu_y_mat[i,1:n_traits] = current_expected_state';
        }
        return mu_y_mat;
      }
matrix calc_covM(int N, int n_traits, matrix t_MRCA_tips, matrix t_root_MRCA, matrix F_mat, matrix Sigma, matrix lambdas_sum_matrix, matrix P, matrix inv_P){
  matrix[N * n_traits, N * n_traits] cov_M;
  matrix[n_traits,n_traits] M_eigen;
  real tiny_positive_jitter_for_integral = 1e-6; // Ensure this is defined and used

  for(i in 1:N){
    for(j in 1:N){
      matrix[n_traits,n_traits] term1 = exp(-F_mat*t_MRCA_tips[i,j]);
      matrix[n_traits,n_traits] term2 = exp(-F_mat'*t_MRCA_tips[j,i]);

      if(t_root_MRCA[i,j] != 0){
        M_eigen = ((1 / lambdas_sum_matrix * (1 - exp(-lambdas_sum_matrix*t_root_MRCA[i,j]))));
      } else {
        // This is the crucial 'else' branch where M_eigen becomes zero
        M_eigen = rep_matrix(0, n_traits, n_traits);
      }

      matrix[n_traits,n_traits] M_transformed_Sigma = inv_P  * Sigma * inv_P';

      matrix[n_traits,n_traits] M_integral = P * (M_eigen .* M_transformed_Sigma) * P';

      // --- THIS IS THE NEW CRITICAL LINE WE ADDED ---
      // If M_eigen was zero, M_integral will be zero. Add a tiny base covariance.
      // This ensures the covariance contributed by common history is never exactly zero
      // and always has a minimum positive definite component.
      if (t_root_MRCA[i,j] == 0) {
          M_integral += diag_matrix(rep_vector(tiny_positive_jitter_for_integral, n_traits));
      }
      // ---------------------------------------------


      matrix[n_traits,n_traits] cov_traits_ij = term1 * M_integral * term2;

      int start_row = (i - 1) * n_traits + 1;
      int end_row = i * n_traits;
      int start_col = (j - 1) * n_traits + 1;
      int end_col = j * n_traits;

      cov_M[start_row:end_row,start_col:end_col] = cov_traits_ij;
      }
    }
    return cov_M;
  }
}
data {
  int<lower=1> N;
  int<lower=1>  n_traits;
  int<lower=1> n_branches;
  int<lower=1> n_path_nodes;
  int<lower=1> n_regs;

  matrix[N,n_traits] y_obs;
  matrix[N,n_traits] y_error;

  vector[n_path_nodes] path_nodes;
  array[N] int path_nodes_idx;
  array[n_branches] int regimes;
  array[N] int regimes_idx;
  vector[n_branches] branch_lengths;
  array[N] int branch_lengths_idx;
  matrix[n_regs,n_traits] OU_mat;

  matrix[N,N] t_MRCA_tips;
  matrix[N,N] t_root_MRCA;

  //real log_lambda_min_val; // New data variable for the lower bound of log_lambdas


}
parameters {
  matrix[N,n_traits] y_true_raw;     // True x values
  vector[n_traits] mu_y_true_param;     // True x values
  vector[n_traits] log_sigma_y_true;     // True x values

  vector[n_traits] y_0_ancestral_raw; //Ancestral value
  vector[n_traits] mu_y_0_ancestral;
  vector[n_traits] log_sigma_y_0_ancestral;
  //vector<lower=0>[n_traits] lambdas; //Eigenvalues of F_mat
  //vector<lower=0, upper=2>[n_traits] lambdas; // Constrain lambdas, e.g., max 5.0
                                            // Adjust upper=5 as a starting point.
                                            // If still failing, try smaller, like upper=2.
                                            // The exact value depends on your tree height.
  vector<lower=log(1e-6)>[n_traits] log_lambdas; // Apply the lower bound
  corr_matrix[n_traits] Omega_P;     // Correlation matrix for P (part of LKJ prior)
  vector<lower=1e-7>[n_traits] sigma_P; // Standard deviations for P (part of LKJ prior)
  corr_matrix[n_traits] Omega_Sigma;     // Correlation matrix for P (part of LKJ prior)
  vector<lower=1e-7>[n_traits] sigma_Sigma; // Standard deviations for P (part of LKJ prior)
}
transformed parameters{
  matrix[n_traits,n_traits] Sigma = diag_matrix(sigma_Sigma) * Omega_Sigma * diag_matrix(sigma_Sigma);
  matrix[n_traits,n_traits] P = diag_matrix(sigma_P) * Omega_P * diag_matrix(sigma_P);
  matrix[n_traits,n_traits] inv_P =  inverse(P);
  vector<lower=0>[n_traits] lambdas = exp(log_lambdas);
  matrix[n_traits,n_traits] F_mat = P * diag_matrix(lambdas) * inv_P;
  matrix[n_traits, n_traits] lambdas_sum_matrix;
  vector[n_traits] y_0_ancestral;
  vector<lower=0>[n_traits] sigma_y_true = exp(log_sigma_y_true);     // True x values
  vector<lower=0>[n_traits] sigma_y_0_ancestral = exp(log_sigma_y_0_ancestral);     // True x values

  matrix[N,n_traits] y_true; // The actual y_true values
  for(j in 1:n_traits){ // Apply transformation for each trait (column)
    y_true[,j] = mu_y_true_param[j] + y_true_raw[,j] * sigma_y_true[j];
  }

  for (i in 1:n_traits) {
    y_0_ancestral[i] = mu_y_0_ancestral[i] + y_0_ancestral_raw[i] * sigma_y_0_ancestral[i];
  }

  for (k in 1:n_traits) {
    for (l in 1:n_traits) {
      lambdas_sum_matrix[k, l] = lambdas[k] + lambdas[l];
    }
  }

  matrix[N,n_traits] mat_mu_y_true = calc_mu_y_mat(N, n_traits, path_nodes, path_nodes_idx, regimes, regimes_idx,
      branch_lengths, branch_lengths_idx, F_mat, OU_mat, y_0_ancestral);

  //print("--- DEBUGGING V CALCULATION AT INITIALIZATION ---");
//  print("Current lambdas: ", lambdas);
//  print("Current sigma_P: ", sigma_P);
//  print("Current sigma_Sigma: ", sigma_Sigma);
//  print("Current Omega_P: ", Omega_P);
//  print("Current Omega_Sigma: ", Omega_Sigma);
//  print("Derived P: ", P);
//  print("Derived inv_P: ", inv_P);
//  print("Derived F_mat: ", F_mat);
//  print("Derived Sigma: ", Sigma);
//  print("Derived lambdas_sum_matrix: ", lambdas_sum_matrix);

  matrix[N*n_traits,N*n_traits] V = calc_covM(N, n_traits, t_MRCA_tips, t_root_MRCA, F_mat, Sigma, lambdas_sum_matrix, P, inv_P);

  // Add a tiny positive value to the diagonal for numerical stability
  // Start with a small value like 1e-6, if it fails, try 1e-5, etc.
  // The scale of your 'cov_M' values (around 0.003-0.001) suggests 1e-6 or 1e-7 might be appropriate.
  for (k in 1:(N*n_traits)) {
    V[k, k] += 10; // Keep the large jitter (0.1) as we found it necessary for positive definiteness
    for (l in (k + 1):(N*n_traits)) { // Loop only over elements above the diagonal
      V[l, k] = V[k, l]; // Make V[l,k] (lower triangle) exactly equal to V[k,l] (upper triangle)
    }
  }


  // --- NEW DEBUG PRINTS FOR V BEFORE CHOLESKY ---
 // print("--- DEBUG FINAL V BEFORE CHOLESKY ---");
  //print("Dims of V: ", dims(V));
  //print("V values (first 6x6 block): ", block(V, 1, 1, min(N*n_traits, 6), min(N*n_traits, 6)));
  //print("V values (last 6x6 block): ", block(V, N*n_traits - 5, N*n_traits - 5, 6, 6));

  // --- IMPORTANT: Check eigenvalues of V *after* jitter and symmetry enforcement ---
  //print("--- DEBUG EIGENVALUES OF V (AFTER JITTER & SYMMETRY) ---");
  //print(eigenvalues_sym(V)); // This line prints all eigenvalues of the symmetric matrix V

  vector[N*n_traits] y_true_vec = to_vector(y_true); //
  vector[N*n_traits] expected_y_tip_values = to_vector(mat_mu_y_true); //Squash matrix
  matrix[N*n_traits, N*n_traits] L_V = cholesky_decompose(V);

}

model {
  mu_y_0_ancestral ~ normal(mean(y_obs), sd(y_obs));   // Prior for the mean of log_sd_theta values. Centered on your previous fixed value, but could be wider.
  //sigma_y_0_ancestral ~ exponential(2);
  y_0_ancestral_raw ~ std_normal();          // Each raw effect ~ N(0,1) - essentially Z-score
  log_sigma_y_0_ancestral ~ normal(log(0.5),1.0);

  to_vector(y_true_raw) ~ std_normal(); // Compact way to specify N(0,1) for all elements
  mu_y_true_param ~ normal(mean(y_obs), sd(y_obs));

  //sigma_y_true ~ exponential(2); // Prior for the standard deviations (mean of 1.0, positive)
  log_sigma_y_true ~ normal(log(0.5),1.0);
  log_lambdas ~ normal(log(0.5),1.0);
  //lambdas ~ exponential(2);
  Omega_P ~ lkj_corr(10); //Prior for correlation matrix for P
  //sigma_P ~ cauchy(0, 2.5); //Prior for standard deviations for P
  sigma_P ~ normal(1.0, 0.5); //Prior for standard deviations for P
  Omega_Sigma ~ lkj_corr(10); //Prior for correlation matrix for Sigma
  //sigma_Sigma ~ cauchy(0, 2.5); //Prior for standard deviations for Sigma
  sigma_Sigma ~ normal(1.0, 0.5); //Prior for standard deviations for Sigma

  for(n in 1:n_traits){ // Measurement error
    y_obs[,n] ~ normal(y_true[,n], y_error[,n]);}


  y_true_vec ~ multi_normal_cholesky(expected_y_tip_values, L_V); //Likelihood function - how likely is it that we observe y_true_vec given the parameters?
  //mu_y_true = expected trait values at the tips of the tree, given the ancestral state (often an optimum) and the evolutionary parameters of the OU process.

}
generated quantities {
  vector[N*n_traits] log_lik;
  int n = 1;
  for(j in 1:n_traits){
    for (i in 1:N) {
      log_lik[n] = normal_lpdf(y_obs[i,j] | y_true[i,j], y_error[i,j]); //For each specific observed data point x_obs[i,j],
      //it's calculating how well that observed value x_obs[i,j] matches the x_true[i,j] value that was sampled
      //from the posterior in that particular MCMC iteration, taking into account the known measurement error x_error[i,j].
      n = n + 1;
    }


  }
}

