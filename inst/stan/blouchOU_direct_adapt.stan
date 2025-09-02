//Blouch OU model - 2025 - Blouch v2.0
//Combination of direct effect and adaptive predictors
//Model accounts for measurement error on X and/or Y variables
//See below for which lines to comment/uncomment based on whether measurement error is present

functions {
  matrix calc_mixed_dmX(real a, vector T_term, matrix X, int Z_direct, int Z_adaptive){
    int N = dims(X)[1];
    int Z = dims(X)[2];
    vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
    matrix[N,Z_adaptive] rhos = rep_matrix(rho,Z_adaptive);
    matrix[N,Z] dmX = append_col(X[,1:Z_direct],X[,Z_direct+1:Z_adaptive+Z_direct] .* rhos);
    return(dmX);
  }
  matrix calc_direct_V( real a,real sigma2_y,matrix ta, matrix tij) {
        int N = dims(ta)[1];
        matrix[N, N] Vt;
        Vt = sigma2_y /( 2 * a) * ((1 - exp(-2 * a * ta)) .* exp(-a * tij)); //From Hansen (1997)
        return Vt;
  }
  matrix calc_adaptive_V(real a,real sigma2_y,matrix ta, matrix tij, matrix tja, vector T_term, vector beta, matrix sigma2_x) {
    int N = dims(ta)[1];
    int Z = dims(beta)[1];
    vector[Z] ones = rep_vector(1,Z);
    matrix[N,N] ti = rep_matrix(T_term,N);
    matrix[N,N] term0;
    matrix[N,N] term1;
    matrix[N,N] term2;
    matrix[N,N] Vt;
    real var_opt;
    if(Z==1){var_opt = beta[1] * beta[1] * sigma2_x[1,1];
    }else{var_opt = beta[1:Z]' * sigma2_x * ones;}
    term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
    term1 = (1 - exp(-a * ti)) ./ (a * ti);
    term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
    Vt = term0 + var_opt * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
    return Vt;
  }
}
data {
  int N;
  int Z_direct;
  int Z_adaptive;
  int Z_X_error;
  vector[N] Y_obs;
  matrix[N,Z_direct+Z_adaptive] X_obs;
  vector[N] Y_error;
  matrix[N,Z_X_error] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[Z_adaptive,Z_adaptive] sigma2_x;
  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
  vector[2] beta_prior;

}
parameters {
  real <lower = 0> hl;
  vector[Z_direct+Z_adaptive] beta; //Assuming a positive relationship among traits
  real optima;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z_direct+Z_adaptive] X;
}
transformed parameters{
  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  matrix[N,Z_direct+Z_adaptive] dmX = calc_mixed_dmX(a,T_term,X,Z_direct,Z_adaptive); //Given measurement error in X variable, uncomment this statement
  matrix[N,N] V = calc_adaptive_V(a,sigma2_y,ta,tij,tja,T_term,beta[(Z_direct+1):(Z_adaptive+Z_direct)],sigma2_x);
  matrix[N,N] L_v = cholesky_decompose(V);
  vector[N] mu = optima+dmX*beta;
}
model {
  target += lognormal_lpdf(hl|hl_prior[1],hl_prior[2]);
  target += exponential_lpdf(vy|vy_prior);
  target += normal_lpdf(optima|optima_prior[1],optima_prior[2]);
  target += normal_lpdf(beta|beta_prior[1],beta_prior[2]);
  for(i in 1:(Z_direct+Z_adaptive)){//Given measurement error in X variable, uncomment this nested statement
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
  vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
  vector[Z_adaptive] beta_e;
  for(i in 1:Z_adaptive){
    beta_e[i] = beta[Z_direct+i]* rho[i];
    }
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

