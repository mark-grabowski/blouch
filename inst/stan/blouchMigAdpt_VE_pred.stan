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

  vector[N_total] x_obs;    // Environmental variables
  int N_x_miss; //Number of missing trait means
  int x_missidx[N_x_miss]; //Index of which x's are missing from full dataset

  vector<lower=0>[N_total] z_error; // SD for measurement error in z - staring with complete case
  vector<lower=0>[N_total] x_error; // SD for measurement error in x - starting with complete case

  //array[2] matrix[N_pops, N_pops] M; #Migration matrices
  matrix[N_pops, N_pops] M;

  vector[2] alpha_prior; // Prior for alpha (mean, sd)
  vector[2] beta_prior;    // Prior for beta (mean, sd)
  vector[2] es_prior;
  real<lower=1> nu_cor;   // Shape parameter for LKJ prior on Omega (must be >= 1)

  array[N_total] int<lower=1, upper=N_regimes> reg_idx;
  array[N_total] int<lower=1, upper=N_pops> pop_idx;

}

parameters {
  vector[N_z_miss] z_impute; #Imputed z values
  vector[N_x_miss] x_impute; #Iputed x values
  real alpha;
  real beta;
  //vector[N_regimes] alpha;         // Intercept
  //vector[N_regimes] beta;           // Slope
  real<lower=0> es;      // Individual-level es
  //vector<lower=0>[N] es;      // Population-specific es
  vector<lower=0.00001>[N_pops] sd_theta; // Standard deviations for theta deviations
  vector<lower=0.00001>[N_pops] sd_e;     // Standard deviations for residual deviations

  cholesky_factor_corr[N_pops] L_Omega_theta; // Cholesky factor of correlation matrix for theta
  cholesky_factor_corr[N_pops] L_Omega_e;     // Cholesky factor of correlation matrix for e
  vector[N_pops] z_true_pops;     // True z values for 20 populations, taking account measurement error due to Greenhouse vs. Field
  vector[N_pops] x_true;     // True x values

}

transformed parameters {
  vector[N_total] z_merge;
  vector[N_total] x_merge;

  matrix[N_pops,N_pops] A;       // Adaptation matrix
  matrix[N_pops,N_pops] AM_term; // I + A^-1 M
  vector[N_pops] mu_pops;        // Mean of y
  matrix[N_pops,N_pops] Sigma;   // Variance-covariance matrix of y
  matrix[N_pops, N_pops] Omega_theta = L_Omega_theta * L_Omega_theta'; // Recover correlation matrix for theta
  matrix[N_pops, N_pops] Omega_e = L_Omega_e * L_Omega_e';             // Recover correlation matrix for e
  //matrix[N, N] Omega_theta = diag_matrix(rep_vector(1.0, N)); //Diagonal matrix
  //matrix[N, N] Omega_e   = diag_matrix(rep_vector(1.0, N)); //Diagonal matrix

  z_merge = merge_missing(z_missidx, to_vector(z_obs), z_impute); #Make full vector by merging imputed and real values
  x_merge = merge_missing(x_missidx, to_vector(x_obs), x_impute);


  // Construct A
  for (i in 1:N_pops) {
    A[i, i] = -es * (1 + M[i, i]); //Single es
    //A[i, i] = -es[i] * (1 + M[i, i]); //Population-specific es
    for (j in 1:N_pops) {
      if (i != j) {
        A[i, j] = 0;
      }
    }
  }

  AM_term = diag_matrix(rep_vector(1, N_pops)) + inverse(A) * M;

  for(i in 1:N_pops){
    mu_pops[i] = alpha + x_true[i] * beta;
    }


  matrix[N_pops, N_pops] diag_sd_theta = diag_matrix(sd_theta);
  matrix[N_pops, N_pops] diag_sd_e = diag_matrix(sd_e);
  //matrix[N, N] diag_sd_theta = diag_matrix(rep_vector(sd_theta, N));
  //matrix[N, N] diag_sd_e = diag_matrix(rep_vector(sd_e, N));

  Sigma = inverse(AM_term) * (diag_sd_theta * Omega_theta * diag_sd_theta) * (inverse(AM_term)') + (diag_sd_e * Omega_e * diag_sd_e);
  vector[N_pops] mu_eq_pops = inverse(AM_term) * mu;

}

model {
  // Priors
  alpha ~ normal(alpha_prior[1], alpha_prior[2]);
  beta ~ normal(beta_prior[1], beta_prior[2]);
  es ~ normal(es_prior[1], es_prior[2]) T[0,]; //Prior to implementing hierarchical prior on es; T[0,] truncates the normal distribution at zero, ensuring es remains positive

  sd_theta ~ lognormal(log(0.1), 0.5);  // Mean around 0.1 on the original scale, smaller SD on log scale
  sd_e ~ lognormal(log(0.1), 0.5);  // Mean around 0.1 on the original scale, smaller SD on log scale

  // Prior on the Cholesky factor of the correlation matrix
  L_Omega_theta ~ lkj_corr_cholesky(nu_cor);
  L_Omega_e ~ lkj_corr_cholesky(nu_cor);

  // Measurement error - also combines z from field and greenhouse into one vector 20 populations long
  for(n in 1:N_total){
    z_merge_obs[n] ~ normal(z_true_pops[pops_idx[n]], z_error[n]);
    #Taking all oberved values and combining to distrubution with error the result of Field versus greenhouse
  }


  x_merge ~ normal(x_true, x_error); #Estimate x_true based on merged data and known SE

  // Likelihood
  //z_true ~ multi_normal(mu, Sigma);
  matrix[N, N] L_Sigma = cholesky_decompose(Sigma);//For speed increase
  z_true_pops ~ multi_normal_cholesky(mu_eq_pops, L_Sigma);

}

generated quantities {
  vector[N_total] log_lik;
  vector[N_total] var_z;
  matrix[N_total, N_total] Cor_theta = Omega_theta;
  matrix[N_total, N_total] Cor_e = Omega_e;

  for (n in 1:N_total) {
    var_z[n] = Sigma[n, n];
    log_lik[n] = normal_lpdf(z_obs[n] | z_true[n], z_error[n]); //The log_lik calculation then simply evaluates how likely the actual observation z_obs[n] was, given that specific draw's value for z_true[n] and the known measurement noise z_error[n].

  }


}
