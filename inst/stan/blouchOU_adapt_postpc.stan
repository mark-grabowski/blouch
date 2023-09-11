//Blouch OU model reprogrammed - 2023
//Adaptive model accounts for measurement error on X and/or Y variables
//See below for which lines to comment/uncomment based on whether measurement error is present

functions {
  matrix calc_dmX(real a, vector T_term, matrix X){
    int N = dims(X)[1];
    int Z = dims(X)[2];
    vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
    matrix[N,Z] rhos = rep_matrix(rho,Z);
    matrix[N,Z] dmX = X .* rhos;
    return(dmX);
  }
  matrix calc_V(real a,real sigma2_y,matrix ta, matrix tij, matrix tja, vector T_term, vector beta, matrix sigma2_x) {
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
//    }else{var_opt = (beta[1:Z] .* beta[1:Z])' * sigma2_x * ones;}
    }else{var_opt = (beta[1:Z])' * sigma2_x * ones;}
    term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
    term1 = (1 - exp(-a * ti)) ./ (a * ti);
    term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
    Vt = term0 + var_opt * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
    return Vt;
  }
}
data {
  int N;
  int Z;
  vector[N] Y_obs;
  matrix[N,Z] X_obs;
  vector[N] Y_error;
  matrix[N,Z] X_error;
  matrix[N,N] ta;
  matrix[N,N] tij;
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[Z,Z] sigma2_x;

}
parameters {
  real <lower = 0> hl;
  vector<lower=0>[Z] beta;
  real alpha;
  //real <lower=0> sigma2_y;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z] X;
}
transformed parameters{
}
model {
  vector[N] mu;
  real a;
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z] dmX;
  real sigma2_y = vy*(2*(log(2)/hl));
  hl ~ lognormal(log(0.25),0.75);
  vy ~ exponential(5);
  alpha ~ normal(2,0.2); //intercept from OLS
  beta ~ normal(0,0.25);
  a = log(2)/hl;
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta,sigma2_x);
  L_v = cholesky_decompose(V);
  for(i in 1:Z){ //Given measurement error in X variable, uncomment this nested statement
    X[,i] ~ normal(0,1);
    X_obs[,i] ~ normal(X[,i], X_error[,i]);
  }
  dmX = calc_dmX(a,T_term,X);//Given measurement error in X variable, uncomment this statement
  //dmX = calc_dmX(a,T_term,X_obs);//Given no measurement error in X variable, uncomment this statement
  mu = alpha+dmX*beta;
  Y ~ multi_normal_cholesky(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  Y_obs ~ normal(Y,Y_error); //Given measurement error in Y variable, uncomment this statement
  //Y_obs ~ multi_normal_cholesky(mu , L_v); //Given no measurement error in Y variable, uncomment this statement
}
generated quantities {
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z] dmX;
  vector[N] mu;
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta,sigma2_x);
  L_v = cholesky_decompose(V);
  dmX = calc_dmX(a,T_term,X);//Given measurement error in X variable, uncomment this statement
  mu = alpha+dmX*beta;
  Y_sim = multi_normal_cholesky_rng(mu , L_v);//Given measurement error in Y variable, uncomment this statement
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]); //Given measurement error in Y variable, uncomment this statement
  }
}

