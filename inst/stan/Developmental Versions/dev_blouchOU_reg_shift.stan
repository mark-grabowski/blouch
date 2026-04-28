//Blouch OU model - 2025 - Blouch v2.0
//Regime model - for single regime painting and SIMMAPS

functions { //Calculate covariance matrix
//Elements of covariance matrix for GLS depends on rate of adaptation, time separating the two species ($t_{ij}$), and time between their MRCA and the root of the tree ($t_{ra}$)
  matrix calc_V(real a,real sigma2_y,matrix ta, matrix tij) {
    int N = dims(ta)[1];
    matrix[N, N] Vt;
    Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
    return(Vt);
  }


data {
  int N; //Number of tips
  int M; //Number of edges
  matrix[N,M] path_matrix; //1s if branch J is ancestor to tip i, 0 otherwise
  array[M] j //Index for branches
  matrix[N,N] t_MRCA_root; //Time from tip to MRCA of two tips
  matrix[N,N] t_sep_tips; //Time separating tips

  vector[N] Y_obs; //Y observed
  vector[N] Y_error; //Y observed

  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
}

parameters {
  real<lower=0> hl;
  vector[n_reg] optima; //Regime Coefficients
  real <lower=0> vy;
  vector[N] Y;
}
transformed parameters{
  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  matrix[N,N] V = calc_V(a, sigma2_y,ta, tij);
  matrix[N,N] L_v = cholesky_decompose(V);
  matrix[N,n_reg] dmX = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes);
  vector[N] mu = dmX*optima;
}
model {
  target += lognormal_lpdf(hl|hl_prior[1],hl_prior[2]);
  target += exponential_lpdf(vy|vy_prior);
  target += normal_lpdf(optima|optima_prior[1],optima_prior[2]);
  target += multi_normal_cholesky_lpdf(Y | mu , L_v);
  target += normal_lpdf(Y_obs | Y, Y_error);
}
generated quantities {

}
