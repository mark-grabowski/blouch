//Blouch OU model reprogrammed - 2023
//Combination of regime model with direct effect and adaptive predictors
//Regime model - for single regime painting and SIMMAPS
//Multilevel model - varying effects - varying intercepts and slopes across regimes
//See below for which lines to comment/uncomment based on whether measurement error is present

functions {
  int num_matches(vector x, real y) { //Thanks to Stan Admin Jonah -https://discourse.mc-stan.org/t/how-to-find-the-location-of-a-value-in-a-vector/19768/2
  int n = 0;
  for (i in 1:rows(x))
    if (x[i] == y)
      n += 1;
  return(n);
  }
  int[] which_equal(vector x, real y) {
    int match_positions[num_matches(x, y)];
    int pos = 1;
    //for (i in 1:size(x)) {
    for (i in 1:(dims(x)[1])) {
      if (x[i] == y) {
        match_positions[pos] = i;
        pos += 1;
        }
      }
    return(match_positions);
  }
  vector weight_segments(real a, vector t_beginning, vector t_end, real time, int nodes){//Individual lineage, calculate weights per segment
    vector[nodes] weights = append_row(exp(-a * t_beginning) - exp(-a * t_end),exp(-a * time));
    return(weights);
  }

  row_vector weights_regimes(int n_reg, real a, vector t_beginning, vector t_end, real time, vector reg_match, int nodes){//
    //Individual lineage, calculate weights for regimes on each segement
    vector[nodes] weight_seg = weight_segments(a, t_beginning[1:(nodes-1)], t_end[1:(nodes-1)], time, nodes);
    //print(weight_seg);
    vector[n_reg] reg_weights = rep_vector(0,n_reg);
    for(i in 1:n_reg){//reg_match should have values 1,2,3 denoting different regimes
      int ids[num_matches(reg_match, i)] = which_equal(reg_match, i); //Returns indixes of matching regimes in weight_segments vector
      //print(ids);
      //print(weight_seg[ids]);
      reg_weights[i] = sum(weight_seg[ids]);
      //print(reg_weights[i]);
      }
    return(reg_weights');
  }
  matrix calc_optima_matrix(int N, int n_reg, real a, matrix t_beginning, matrix t_end, matrix times, matrix reg_match, int[] nodes){
    matrix[N,n_reg] optima_matrix = rep_matrix(0,N,n_reg);
    for(i in 1:N){ //For each tip/lineage, figure out weighting of regimes
      optima_matrix[i,] = weights_regimes(n_reg, a, t_beginning[i,]', t_end[i,]', times[i,1], reg_match[i,]', nodes[i]);
      //print(i);
      //print(optima_matrix[i,]);
      }
    return(optima_matrix);
  }
  matrix calc_mixed_dmX(real a, vector T_term, matrix X, int Z_direct, int Z_adaptive){
    int N = dims(X)[1];
    int Z = dims(X)[2];
    vector[N] rho = (1 - (1 - exp(-a * T_term))./(a * T_term));
    matrix[N,Z_adaptive] rhos = rep_matrix(rho,Z_adaptive);
    matrix[N,Z] dmX = append_col(X[,1:Z_direct],X[,Z_direct+1:Z_adaptive+Z_direct] .* rhos);
    return(dmX);
  }

  matrix calc_V(real a,real sigma2_y,matrix ta, matrix tij, matrix tja, vector T_term, matrix beta, matrix sigma2_x, int Z_adaptive, int n_reg) {
    int N = dims(ta)[1];
    matrix[Z_adaptive,1] ones = rep_matrix(1,Z_adaptive,1);
    matrix[N,N] ti = rep_matrix(T_term,N);
    matrix[N,N] term0;
    matrix[N,N] term1;
    matrix[N,N] term2;
    matrix[N,N] Vt;
    real var_opt;
    var_opt=sum(square(beta)*sigma2_x);
    term0 = ((var_opt + sigma2_y) / (2 * a)) * (1 - exp( -2 * a * ta)) .* exp(-a * tij);
    term1 = (1 - exp(-a * ti)) ./ (a * ti);
    term2 = exp(-a * tja) .* (1 - exp(-a * ti)) ./ (a * ti);
    Vt = term0 + var_opt * (ta .* term1 .* (term1') - ((1 - exp(-a * ta)) ./ a) .* (term2 + (term2'))); //From Hansen et al. (2008)
    return Vt;
  }
}
data {
  int N; //Numnber of tips
  int n_reg; //Number of regimes
  int Z_direct;
  int Z_adaptive;
  int Z_X_error;
  int max_node_num; //Max number of nodes in lineage
  vector[N] Y_obs; //Y observed
  matrix[N,Z_direct+Z_adaptive] X_obs;
  vector[N] Y_error;
  matrix[N,Z_X_error] X_error;
  matrix[Z_adaptive,Z_adaptive] sigma2_x;
  matrix[N,N] ta; //Time from tip to ancestor
  matrix[N,N] tij;
  matrix[N,N] tja;
  vector[N] T_term;
  matrix[N, max_node_num] t_beginning; //Matrix of times for beginning of segments to node
  matrix[N, max_node_num] t_end; //Matrix of times for end of segments to
  matrix[N, max_node_num] times; //Matrix of root to node times
  matrix[N, max_node_num] reg_match; //Matrix of 1,2,3 denoting each regime for each node in a lineage. 0 if no node
  int nodes[N]; //Vector of number of nodes per lineage
  int reg_tips[N]; //Regimes at the tips
  vector[2] hl_prior;
  real vy_prior;
  vector[2] optima_prior;
  vector[2] beta_prior;
  vector[2] sigma_prior;

}

parameters {
  real<lower=0> hl;
  real <lower=0> vy;
  vector[N] Y;
  matrix[N,Z_direct+Z_adaptive] X;
  real optima_bar;
  real beta_bar;
  vector<lower=0>[1+Z_direct+Z_adaptive] sigma;
  cholesky_factor_corr[(1+Z_direct+Z_adaptive)] L_Rho;
  matrix[(1+Z_direct+Z_adaptive),n_reg] z;
}
transformed parameters{
  vector[n_reg] optima;
  matrix[n_reg,Z_direct+Z_adaptive] beta;
  matrix[n_reg,(1+Z_direct+Z_adaptive)] v;
  v = (diag_pre_multiply(sigma, L_Rho) * z)';
  beta = beta_bar + v[, 2:(Z_direct+Z_adaptive+1)];
  optima = optima_bar + v[, 1];
}
model {
  matrix[N,N] V;
  vector[N] mu;
  matrix[N,N] L_v;
  matrix[N,Z_direct+Z_adaptive] pred_X;
  real a = log(2)/hl;
  real sigma2_y = vy*(2*(log(2)/hl));
  matrix[N,n_reg] optima_matrix;
  //hl ~ lognormal(log(0.25),0.25);
  target += lognormal_lpdf(hl|hl_prior[1],hl_prior[2]);
  //vy ~ exponential(20);
  target += exponential_lpdf(vy|vy_prior);
  //optima_bar ~ normal(mean(Y),1);
  //beta_bar ~ normal(0,0.25);
  target += normal_lpdf(optima_bar|optima_prior[1],optima_prior[2]);
  target += normal_lpdf(beta_bar|beta_prior[1],beta_prior[2]);
  //L_Rho ~ lkj_corr_cholesky(2);
  target += lkj_corr_cholesky_lpdf(L_Rho|2);
  //sigma ~ exponential(5);
  //sigma ~ normal(0,1);
  target += normal_lpdf(sigma|sigma_prior[1],sigma_prior[2]);
  for(i in 1:(Z_direct+Z_adaptive)){//Given measurement error in X variable, uncomment this nested statement
    //X[,i] ~ normal(0,1);
    target += normal_lpdf(X[,i]|0,1);
    //X_obs[,i] ~ normal(X[,i], X_error[,i]);
    target += normal_lpdf(X_obs[,i]|X[,i],X_error[,i]);
  }
  for(i in 1:(1+Z_direct+Z_adaptive)){
    for(j in 1:n_reg){
       target += normal_lpdf(z[i,j]|0,1);
    }
  }
  optima_matrix = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes);
  pred_X = calc_mixed_dmX(a,T_term,X,Z_direct,Z_adaptive);//Given measurement error in X variable, uncomment this statement
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta[,(Z_direct+1):(Z_direct+Z_adaptive)],sigma2_x,Z_adaptive,n_reg);
  L_v = cholesky_decompose(V);
  for(i in 1:N){
    mu[i] = optima_matrix[i,]*optima+pred_X[i,]*beta[reg_tips[i],]';
    }
  target += multi_normal_cholesky_lpdf(Y | mu , L_v);
  target += normal_lpdf(Y_obs | Y, Y_error);

}
generated quantities {
  vector[N] Y_sim;
  vector[N] Y_sim_obs;
  matrix[N,N] V;
  matrix[N,N] L_v;
  matrix[N,Z_direct+Z_adaptive] pred_X;
  matrix[N,n_reg] optima_matrix;
  vector[N] mu;

  real sigma2_y = vy*(2*(log(2)/hl));
  real a = log(2)/hl;
  real rho = (1 - (1 - exp(-a * T_term))./(a * T_term))[1];
  matrix[n_reg,Z_adaptive] beta_e = beta[,(Z_direct+1):(Z_direct+Z_adaptive)]*rho;
  V = calc_V(a,sigma2_y,ta,tij,tja,T_term,beta[,(Z_direct+1):(Z_direct+Z_adaptive)],sigma2_x,Z_adaptive,n_reg);
  L_v = cholesky_decompose(V);
  optima_matrix = calc_optima_matrix(N, n_reg, a, t_beginning, t_end, times, reg_match, nodes);
  pred_X = calc_mixed_dmX(a,T_term,X,Z_direct,Z_adaptive);//Given measurement error in X variable, uncomment this statement

  for(i in 1:N){
    mu[i] = optima_matrix[i,]*optima+pred_X[i,]*beta[reg_tips[i],]';
    }

  Y_sim = multi_normal_cholesky_rng(mu , L_v);
  for(i in 1:N){
    Y_sim_obs[i] = normal_rng(Y_sim[i],Y_error[i]);
  }
}
