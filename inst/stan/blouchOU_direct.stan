//Blouch OU model - 2025 - Blouch v2.0
//Direct effect model
//Using Hansen (1997)
//With Statistical Rethinking ME

functions {
  matrix calc_V( real a,real sigma2_y,matrix ta, matrix tij) {
    int N = dims(ta)[1];
    matrix[N, N] Vt;
    Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
    return Vt;
  }
}
data {
  int N;
  int Z_direct;
  vector[N] Y_obs;
  matrix[N,Z_direct] X_obs;
  vector[N] Y_error;
  matrix[N,Z_direct] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
  vector[2] beta_prior;
}
parameters {
  real <lower = 0> hl;
  vector[Z_direct] beta;
  real optima;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z_direct] X;
}
transformed parameters{
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  matrix[N,N] V = calc_V(a, sigma2_y,ta, tij);
  matrix[N,N] L_v = cholesky_decompose(V);
  vector[N] mu = optima+X*beta;//Given measurement error in X variable, uncomment this statement
}
model {
  target += lognormal_lpdf(hl|hl_prior[1],hl_prior[2]);
  target += exponential_lpdf(vy|vy_prior);
  target += normal_lpdf(optima|optima_prior[1],optima_prior[2]);
  target += normal_lpdf(beta|beta_prior[1],beta_prior[2]);
  for(i in 1:Z_direct){//Given measurement error in X variable, uncomment this nested statement
    target += normal_lpdf(X[,i]|0,1);
    target += normal_lpdf(X_obs[,i]|X[,i],X_error[,i]);
  }
  target += multi_normal_cholesky_lpdf(Y | mu , L_v);
  target += normal_lpdf(Y_obs | Y, Y_error);

}
generated quantities {
  real g_i;
  real sigma_ii;
  real sigma_i;
  real u_i;
  vector[N] log_lik;
  //Based on https://cran.r-project.org/web/packages/loo/vignettes/loo2-non-factorized.htmlloo-cv-for-multivariate-normal-models
  //LOO-CV for multivariate normal models
  matrix[N,N] V_total = V + diag_matrix(square(Y_error));
  matrix[N,N] inv_V = inverse(V_total);
  for(i in 1:N){
      g_i = (inv_V*(Y_obs-mu))[i];
      sigma_ii = inv_V[i,i];
      u_i = Y_obs[i]-g_i/sigma_ii;
      sigma_i = 1/sigma_ii;

      log_lik[i] = -0.5*log(2*pi()*sigma_i)-0.5*(square(Y_obs[i]-u_i)/sigma_i);
      }
}
