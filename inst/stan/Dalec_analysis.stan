//Varying intercepts model - estimates different alpha and beta values per regime
//Single es parameter, population specific alphas and betas across regimmes (Greenhouse vs. field here)
//Prior to implementing hierarchical priors
//Deals with missing data

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
  int<lower=1> N_total;       // Total number of observations - 20 F, 20 GH
  int<lower=1> N_pops;       // Number of populations
  int<lower=1> N_regimes;

  vector[N_total] z_obs;    // Observed trait means - witout missing values - Field is first here, then greenhouse
  int N_z_miss; //Number of missing trait means
  int z_missidx[N_z_miss]; //Index of which z's are missing from full dataset

  vector[N_pops] x_obs;    // Environmental variables
  int N_x_miss; //Number of missing trait means
  int x_missidx[N_x_miss]; //Index of which x's are missing from full dataset

  vector[N_total] z_error_obs; // SD for measurement error in z - staring with complete case
  vector<lower=0>[N_pops] x_error; // SD for measurement error in x - starting with complete case

  //array[2] matrix[N_pops, N_pops] M; Migration matrices
  matrix[N_pops, N_pops] M;

  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
  vector[2] es_prior;
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)

  array[N_total] int<lower=1, upper=N_regimes> reg_idx;
  array[N_total] int<lower=1, upper=N_pops> pop_idx;

}

parameters {
  vector<lower=0>[N_z_miss] z_impute; //Imputed z values
  vector<lower=0>[N_z_miss] z_error_impute; //Iputed z error (SD) values
  vector[N_x_miss] x_impute; //Iputed x values
  real alpha;
  real beta;
  real alpha_OLS;
  real beta_OLS;

  //vector[N_regimes] alpha;         // Intercept
  //vector[N_regimes] beta;           // Slope
  real<lower=0> es;      // Individual-level es
  //vector<lower=0>[N] es;      // Population-specific es
  vector<lower=0.00001>[N_pops] sd_theta; // Standard deviations for theta deviations
  vector<lower=0.00001>[N_pops] sd_e;     // Standard deviations for residual deviations

  //cholesky_factor_corr[N_pops] L_Omega_theta; // Cholesky factor of correlation matrix for theta
  //cholesky_factor_corr[N_pops] L_Omega_e;     // Cholesky factor of correlation matrix for e
  vector<lower=0>[N_pops] z_true_pops;     // True z values for 10 populations, taking account measurement error due to Greenhouse vs. Field
  vector[N_pops] x_true_pops;     // True x values
  //vector<lower=0>[N_pops] z_obs_sd; // NOW A VECTOR, one SD per population
  vector<lower=0>[N_pops] x_obs_sd; // Assuming x also has pop-specific SDs


}

transformed parameters {
  vector[N_total] z_merge_obs;
  vector[N_total] z_merge_error;

  vector[N_pops] x_merge_obs;

  matrix[N_pops,N_pops] A;       // Adaptation matrix
  matrix[N_pops,N_pops] AM_term; // I + A^-1 M
  vector[N_pops] mu_pops;        // Mean of y for each population
  matrix[N_pops,N_pops] Sigma;   // Variance-covariance matrix of y
  //matrix[N_pops, N_pops] Omega_theta = L_Omega_theta * L_Omega_theta'; // Recover correlation matrix for theta
  //matrix[N_pops, N_pops] Omega_e = L_Omega_e * L_Omega_e';             // Recover correlation matrix for e
  matrix[N_pops, N_pops] Omega_theta = diag_matrix(rep_vector(1.0, N_pops)); //Diagonal matrix
  matrix[N_pops, N_pops] Omega_e   = diag_matrix(rep_vector(1.0, N_pops)); //Diagonal matrix

  z_merge_obs = merge_missing(z_missidx, to_vector(z_obs), z_impute); //Make full vector by merging imputed and real values
  z_merge_error = merge_missing(z_missidx, to_vector(z_error_obs), z_error_impute);
  x_merge_obs = merge_missing(x_missidx, to_vector(x_obs), x_impute);


  // Construct A
  for (i in 1:N_pops) {
    A[i, i] = -es * (1 + M[i, i]); //Single es
    //A[i, i] = -0.1 * (1 + M[i, i]); //es set at 0.1
    //A[i, i] = -5 * (1 + M[i, i]); //es set at 5
    //A[i, i] = -es[i] * (1 + M[i, i]); //Population-specific es
    for (j in 1:N_pops) {
      if (i != j) {
        A[i, j] = 0;
      }
    }
  }

  AM_term = diag_matrix(rep_vector(1, N_pops)) + inverse(A) * M;

  for(i in 1:N_pops){
    mu_pops[i] = alpha + x_true_pops[i] * beta; //Assume single alpha and beta across Field and Greenhouse data
    }


  matrix[N_pops, N_pops] diag_sd_theta = diag_matrix(sd_theta);
  matrix[N_pops, N_pops] diag_sd_e = diag_matrix(sd_e);
  //matrix[N, N] diag_sd_theta = diag_matrix(rep_vector(sd_theta, N));
  //matrix[N, N] diag_sd_e = diag_matrix(rep_vector(sd_e, N));

  Sigma = inverse(AM_term) * (diag_sd_theta * Omega_theta * diag_sd_theta) * (inverse(AM_term)') + (diag_sd_e * Omega_e * diag_sd_e);

  vector[N_pops] mu_eq_pops = inverse(AM_term) * mu_pops;

}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  beta ~ normal(beta_prior[1], beta_prior[2]);
  alpha_OLS ~ normal(alpha_prior[1], alpha_prior[2]);
  beta_OLS ~ normal(beta_prior[1], beta_prior[2]);

  //es ~ normal(es_prior[1], es_prior[2]) T[0,]; //Prior to implementing hierarchical prior on es; T[0,] truncates the normal distribution at zero, ensuring es remains positive
  //es ~ normal(0,1) T[0,];
  es ~ lognormal(log(0.1),0.5); //Ensures positive es
  for(i in 1:N_pops){
    sd_theta[i] ~ lognormal(log(0.5), 0.2);  // Mean around 0.1 on the original scale, smaller SD on log scale
    sd_e[i] ~ lognormal(log(0.5), 0.2);  // Mean around 0.1 on the original scale, smaller SD on log scale
  }


  //z_obs_sd ~ lognormal(log(1), 0.5); // Median 0.7, SD of log(SD) is 0.5
  x_obs_sd ~ lognormal(log(1), 0.5); // Median 0.7, SD of log(SD) is 0.5

  //sd_theta ~ normal(0,1) T[0,]; // Mean around 0.1 on the original scale, smaller SD on log scale
  //sd_e ~ normal(0,1) T[0,];  // Mean around 0.1 on the original scale, smaller SD on log scale

  // Prior on the Cholesky factor of the correlation matrix
  //L_Omega_theta ~ lkj_corr_cholesky(nu_cor);
  //L_Omega_e ~ lkj_corr_cholesky(nu_cor);

  //z_impute ~ lognormal(log(25),0.5); //Mean of GA
  //x_impute ~ lognormal(log(12),0.5); //Mean of slp

  z_impute ~ normal(25,0.5); //Mean of GA
  z_error_impute ~ lognormal(log(1), 0.5); // Median 0.7, SD of log(SD) is 0.5
  x_impute ~ normal(0,1); //Mean of sl

  // Measurement error - also combines z from field and greenhouse into one vector 20 populations long
  for(n in 1:N_total){
    z_merge_obs[n] ~ normal(z_true_pops[pop_idx[n]], z_merge_error[n]);
    //z_merge_obs[n] ~ normal(z_true_pops[pop_idx[n]], z_obs_sd[pop_idx[n]]);
  }

  for(n in 1:N_pops){
    x_merge_obs[n] ~ normal(x_true_pops[n], x_obs_sd[n]);
  }

  // Likelihood
  //z_true ~ multi_normal(mu, Sigma);
  matrix[N_pops, N_pops] L_Sigma = cholesky_decompose(Sigma);//For speed increase
  z_true_pops ~ multi_normal_cholesky(mu_eq_pops, L_Sigma);

}

generated quantities {
  // Keep these as they were:
  vector[N_total] log_lik;
  vector[N_pops] var_z_pops;
  matrix[N_pops, N_pops] Cor_theta = diag_matrix(rep_vector(1.0, N_pops));
  matrix[N_pops, N_pops] Cor_e = diag_matrix(rep_vector(1.0, N_pops));

  // Declare variables for determinants only
  real det_A;
  real det_AM_term;
  real det_K_term;
  real det_V_term;
  real det_Sigma_term1;
  real det_Sigma;

  // Recalculate components (same as before, these are crucial for the determinants)
  matrix[N_pops,N_pops] A_gq;
  matrix[N_pops,N_pops] AM_term_gq;
  // Omega_theta_gq and Omega_e_gq are already diagonal as per your transformed parameters
  matrix[N_pops, N_pops] Omega_theta_gq = diag_matrix(rep_vector(1.0, N_pops));
  matrix[N_pops, N_pops] Omega_e_gq = diag_matrix(rep_vector(1.0, N_pops));
  matrix[N_pops, N_pops] diag_sd_theta_gq = diag_matrix(sd_theta);
  matrix[N_pops, N_pops] diag_sd_e_gq = diag_matrix(sd_e);

  for (i in 1:N_pops) {
    A_gq[i, i] = -es * (1 + M[i, i]);
    for (j in 1:N_pops) {
      if (i != j) {
        A_gq[i, j] = 0;
      }
    }
  }
  AM_term_gq = diag_matrix(rep_vector(1, N_pops)) + inverse(A_gq) * M;

  matrix[N_pops, N_pops] K_term_gq = diag_sd_theta_gq * Omega_theta_gq * diag_sd_theta_gq;
  matrix[N_pops, N_pops] V_term_gq = diag_sd_e_gq * Omega_e_gq * diag_sd_e_gq;

  // IMPORTANT: The following lines will still cause an error if A_gq or AM_term_gq
  // are singular when inverse() is called. If that's the case, the sampler will halt
  // and you will see the error related to inverse().
  matrix[N_pops, N_pops] Sigma_term1_gq = inverse(AM_term_gq) * K_term_gq * inverse(AM_term_gq)';
  matrix[N_pops, N_pops] Sigma_gq = Sigma_term1_gq + V_term_gq;

  // Now assign the determinant diagnostics
  det_A = determinant(A_gq);
  det_AM_term = determinant(AM_term_gq);
  det_K_term = determinant(K_term_gq);
  det_V_term = determinant(V_term_gq);
  det_Sigma_term1 = determinant(Sigma_term1_gq);
  det_Sigma = determinant(Sigma_gq);



  // Other generated quantities calculations
  for(n in 1:N_pops){
    var_z_pops[n] = Sigma_gq[n, n];
  }

  for (n in 1:N_total) {
    //log_lik[n] = normal_lpdf(z_merge_obs[n] | z_true_pops[pop_idx[n]], z_obs_sd[pop_idx[n]]);
    log_lik[n] = normal_lpdf(z_merge_obs[n] | z_true_pops[pop_idx[n]], z_merge_error[pop_idx[n]]);
  }
}



