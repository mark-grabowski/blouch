//mvBlouch mvOUOU model - 2025
//Multivariate Brownian motion model accounts for measurement error on X and/or Y variables
//After Revell and Harmon 2008, Clavel et al. 2015

functions {
  matrix kprod(matrix A, matrix B) {
    int m_A = rows(A);
    int n_A = cols(A);
    int m_B = rows(B);
    int n_B = cols(B);

    matrix[m_A * m_B, n_A * n_B] Z; // Correct dimensions for the result matrix

    for (i in 1:m_A) {
      for (j in 1:n_A) { // Loop over columns of A
        for (p in 1:m_B) {
          for (q in 1:n_B) {
            int row_in_Z = (i - 1) * m_B + p;
            int col_in_Z = (j - 1) * n_B + q;
            Z[row_in_Z, col_in_Z] = A[i, j] * B[p, q];
          }
        }
      }
    }
    return Z;
  }
}
data {
  int N;
  int n_traits;
  matrix[N,n_traits] x_obs;
  matrix[N,n_traits] x_error;
  matrix[N,N] C; // ta is time from the base of the phylogeny to the most re-cent common ancestor of the two species.
  vector[n_traits] a_hat_prior;
  real sigma_R_prior_sd;
}
parameters {
  matrix[N,n_traits] x_true;     // True x values
  vector[n_traits] a_hat; //Estimated phylogenetic mean
  corr_matrix[n_traits] Omega;     // Correlation matrix for R (part of LKJ prior)
  vector<lower=0>[n_traits] sigma_R; // Standard deviations for R (part of LKJ prior)

}
transformed parameters{
  matrix[n_traits,n_traits] R = diag_matrix(sigma_R) * Omega * diag_matrix(sigma_R);
}
model {
  //Priors
  for(i in 1:n_traits){
    a_hat[i] ~ normal(a_hat_prior[i], 1); //Prior for a_hat based on column means of X traits
  }

  Omega ~ lkj_corr(2); //Prior for correlation matrix for R
  sigma_R ~ normal(0, sigma_R_prior_sd); //Prior for standard deviations for R

  // Measurement error
  for(n in 1:n_traits){
    x_obs[,n] ~ normal(x_true[,n], x_error[,n]);

  }

  vector[N*n_traits] x_true_vec = to_vector(x_true);

  matrix[N, N] L_C = cholesky_decompose(C);
  matrix[n_traits, n_traits] L_R = cholesky_decompose(R);
  matrix[N*n_traits,N*n_traits] L_V = kprod(L_R , L_C);


  matrix[N, n_traits] mat_mu_x_true; // Matrix for mu_x values
  for(i in 1:n_traits){
    for(j in 1:N){
      mat_mu_x_true[j,i] = a_hat[i]; #Fill mu x vactor with a_hat first N are for trait 1, second N for trait 2

    }
  }

  vector[N*n_traits] mu_x_true = to_vector(mat_mu_x_true); //Squash matrix
  x_true_vec ~ multi_normal_cholesky(mu_x_true, L_V);

}
generated quantities {
  vector[N*n_traits] log_lik;
  int n = 1;
  for(j in 1:n_traits){
    for (i in 1:N) {
      log_lik[n] = normal_lpdf(x_obs[i,j] | x_true[i,j], x_error[i,j]); //For each specific observed data point x_obs[i,j],
      //it's calculating how well that observed value x_obs[i,j] matches the x_true[i,j] value that was sampled
      //from the posterior in that particular MCMC iteration, taking into account the known measurement error x_error[i,j].
      n = n + 1;
    }


  }
}

