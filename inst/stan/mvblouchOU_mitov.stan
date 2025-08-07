//mvBlouch mvOU - Mitov approach - fixed predictors model - 2025
//Multivariate OU model accounts for measurement error on X and/or Y variables
//After Bartoszek et al. 2012; Clavel et al. 2015; Mitov et al. 2019

functions {
//Pruning algorithm
  real calculate_log_likelihood_lp(
    int n_traits,
    int n_post_order_path_nodes,
    matrix y_true,
    vector y_root,
    array[] int post_order_path_nodes,
    vector branch_lengths,
    array[] int branch_regime_idx, //Path of nodes from tip to root
    array[] int parent_of_node, //Path of nodes from tip to root
    array[] int node_types,
    array[] matrix H_mats,
    array[] matrix P_mats,
    array[] matrix inv_P_mats,
    array[] matrix Sigma_mats,
    array[] vector theta_mats,
    array[] matrix lambdas_sum_mats){

    array[n_post_order_path_nodes] matrix[n_traits,n_traits] L_kernals = rep_array(rep_matrix(0.0, n_traits, n_traits), n_post_order_path_nodes);
    array[n_post_order_path_nodes] vector[n_traits] m_kernals = rep_array(rep_vector(0.0,n_traits),n_post_order_path_nodes);
    array[n_post_order_path_nodes] real r_kernals = rep_array(0.0, n_post_order_path_nodes);

    for(i in 1:n_post_order_path_nodes){ //Loop through nodes in postorder
      int current_node = post_order_path_nodes[i];
      int parent = parent_of_node[current_node];
      real branch_length = branch_lengths[current_node];
      int current_reg = branch_regime_idx[current_node];
      int current_node_type = node_types[current_node];
      matrix[n_traits,n_traits] Phi;
      vector[n_traits] omega;
      matrix[n_traits,n_traits] V_i;
      matrix[n_traits,n_traits] inv_V_i;
      matrix[n_traits,n_traits] C_i;
      vector[n_traits] d_i;
      real f_i;

      if(current_node_type != 0){ //Skip for root node
          Phi = compute_Phi(
          branch_length,
          H_mats[current_reg]);

        omega = compute_omega(
          n_traits,
          branch_length,
          H_mats[current_reg],
          theta_mats[current_reg]);

        V_i = compute_V_i(
          n_traits,
          branch_length,
          P_mats[current_reg],
          inv_P_mats[current_reg],
          Sigma_mats[current_reg],
          lambdas_sum_mats[current_reg]);


   // --- Debugging prints ---
        //print("Iter: ", i);

        //real det_V_i = determinant(V_i);
        //print("Node: ", current_node, " Branch length: ", branch_length, " Node type: ",
        //  current_node_type, " Currrent regime: ", current_reg);
        //print("Determinant V_i: ", det_V_i);

        //if (det_V_i <= 0) {
        //  print("WARNING: V_i is not positive definite at node ", current_node);
        //}

        //vector[n_traits] eig_vals = eigenvalues_sym(V_i);
        //print("Eigenvalues V_i at node ", current_node, ": ", eig_vals);

        inv_V_i = inverse(V_i);
      }


    //omega + Phi *x_j
    if (current_node_type == 1){// Is current node is a tip, use known values
      vector[n_traits] x_i = y_true[current_node,]'; //Observed tip value
      vector[n_traits] resid = x_i - omega; //=Phi*x_j Reguime dependent offset

      //real logdet_term = log_determinant(2 * pi() * V_i);
      //print("Log determinant term at tip node ", current_node, ": ", logdet_term);

      C_i = -0.5 * Phi' * inv_V_i * Phi; //Quadratic coefficient matrix
      d_i = -Phi' * inv_V_i * resid; //Linear coefficient vector
      f_i = -0.5 * dot_product(resid, inv_V_i* resid) - 0.5 * log_determinant(2 * pi() * V_i);
      //= Mahalanobis distance term
      //A log-determinant normalization constant - different from paper? w = x_i - omega?
      L_kernals[parent] += C_i;
      m_kernals[parent] += d_i;
      r_kernals[parent] += f_i;
    }

    if (current_node_type == 2){// Is current node is a internal, use previously estimated values from all children
      matrix[n_traits,n_traits] L = L_kernals[current_node];
      vector[n_traits] m = m_kernals[current_node];
      real r = r_kernals[current_node];

      matrix[n_traits,n_traits] A_i = -0.5 * inv_V_i;
      vector[n_traits] b_i = inv_V_i * omega;
      matrix[n_traits,n_traits] Q_i = A_i + L;

      //real det_Q_i = determinant(Q_i);
      //print("Determinant Q_i at internal node ", current_node, ": ", det_Q_i);

      //if (det_Q_i <= 0){
      //  print("WARNING: Q_i is not positive definite at internal node ", current_node);
      //}

      matrix[n_traits,n_traits] inv_Q_i = inverse(Q_i);
      vector[n_traits] z = b_i+ m;

      C_i = -0.5 * Phi' * inv_Q_i* Phi; //Quadratic coefficient matrix
      d_i = to_vector(-Phi' * inv_Q_i * z); //Linear coefficient vector
      f_i = r - 0.25 * dot_product(z, inv_Q_i * z)
        + 0.5 * log(2 * pi())
        - 0.5* log_determinant(-2 * Q_i);

      L_kernals[parent] += C_i;
      m_kernals[parent] += d_i;
      r_kernals[parent] += f_i;
      }
    }

    //Return the final log-likelihood
    //if (node_type[current_node]==0){// Is current node is the root ...
    int root_node_idx = post_order_path_nodes[n_post_order_path_nodes];
    matrix[n_traits,n_traits] L_root = L_kernals[root_node_idx];
    vector[n_traits] m_root = m_kernals[root_node_idx];
    real r_root = r_kernals[root_node_idx];

    //L_root = L_kernals[current_node];
    //m_root = m_kernals[current_node];
    //r_root = r_kernals[current_node];
    //}
    return -0.5 * quad_form(L_root, y_root) + dot_product(m_root, y_root) + r_root; //Sum log-likelihood information
  }
  matrix compute_Phi(real branch_length, matrix H_mat){
    return matrix_exp(-H_mat * branch_length);
  }

  vector compute_omega(int n_traits, real branch_length, matrix H_mat, vector theta_mat){
    matrix[n_traits,n_traits] I = diag_matrix(rep_vector(1.0, n_traits));
    return ((I - matrix_exp(-H_mat * branch_length)) * theta_mat);
  }

  matrix compute_V_i(int n_traits, real branch_length, matrix P_branch, matrix inv_P_branch,  matrix Sigma_branch, matrix lambdas_sum_matrix_branch){
    matrix[n_traits,n_traits] term1;
    for(i in 1:n_traits){ //Hadamard fraction elementwise
      for(j in 1:n_traits){
        real lambdas_sum = lambdas_sum_matrix_branch[i,j];
        if (fabs(lambdas_sum) > 1e-10) {
          term1[i,j] = (1-exp(-lambdas_sum*branch_length))/lambdas_sum;
        }
        else{
          term1[i,j] = branch_length; //When lambda_sum = 0
        }
      }
    }
    matrix[n_traits,n_traits] term2 = inv_P_branch  * (Sigma_branch * inv_P_branch');
    matrix[n_traits,n_traits] V_i = P_branch * (term1 .* term2) * P_branch';
    return V_i;
    }
}

data {
  int<lower=1> N;
  int<lower=1>  n_traits;
  int<lower=1> n_regs;
  int<lower=1> n_post_order_path_nodes;
  int<lower=1> n_parent_of_nodes;

  matrix[N,n_traits] y_obs;
  matrix[N,n_traits] y_error;

  array[n_post_order_path_nodes] int post_order_path_nodes; //Path of nodes from tip to root
  vector[n_post_order_path_nodes] branch_lengths;
  array[n_post_order_path_nodes] int branch_regime_idx; //Path of nodes from tip to root
  array[n_parent_of_nodes] int parent_of_node; //Path of nodes from tip to root
  array[n_post_order_path_nodes] int node_types;
}
parameters {
    matrix[N,n_traits] y_true;
    vector[n_traits] y_root;
    array[n_regs] cholesky_factor_corr[n_traits] L_Omega_P;     // Correlation matrix for P (part of LKJ prior)
    array[n_regs] vector<lower=0>[n_traits] sigma_P; // Standard deviations for P (part of LKJ prior)
    array[n_regs] cholesky_factor_corr[n_traits] L_Omega_Sigma;     // Correlation matrix for P (part of LKJ prior)
    array[n_regs] vector<lower=0>[n_traits] sigma_Sigma; // Standard deviations for P (part of LKJ prior)
    array[n_regs] vector<lower=0>[n_traits] log_lambdas; // Apply the lower bound
    array[n_regs] vector<lower=0>[n_traits] theta_mats; // Apply the lower bound


}
transformed parameters{
  array[n_regs] matrix[n_traits,n_traits] Sigma_mats;
  array[n_regs] matrix[n_traits,n_traits] P_mats;
  array[n_regs] matrix[n_traits,n_traits] inv_P_mats;
  array[n_regs] vector<lower=0>[n_traits] lambdas_mats;
  array[n_regs] matrix[n_traits,n_traits] H_mats;
  array[n_regs] matrix[n_traits, n_traits] lambdas_sum_mats;

  for(i in 1:n_regs){
    lambdas_mats[i] = exp(log_lambdas[i]);
    matrix[n_traits, n_traits] Omega_Sigma_i = multiply_lower_tri_self_transpose(L_Omega_Sigma[i]);
    Sigma_mats[i] = quad_form_diag(Omega_Sigma_i, sigma_Sigma[i]);

    matrix[n_traits, n_traits] L = diag_pre_multiply(sigma_P[i], L_Omega_P[i]);
    matrix[n_traits, n_traits] L_inv = mdivide_left_tri_low(identity_matrix(n_traits),L);

    P_mats[i] = L * L';             //
    inv_P_mats[i] = inverse(P_mats[i]); //
    H_mats[i] = L * diag_matrix(lambdas_mats[i]) * L_inv;
    for (k in 1:n_traits){
      for (l in 1:n_traits) {
      lambdas_sum_mats[i][k, l] = lambdas_mats[i][k] + lambdas_mats[i][l];
      }
    }
  }
}
model {
  //Priors
  y_true[,1] ~ normal(3.8,1);
  y_true[,2] ~ normal(6.2,1);
  y_root[1] ~ normal(3.8,1);
  y_root[2] ~ normal(6.2,1);

  for(i in 1:n_regs){
    L_Omega_P[i] ~ lkj_corr_cholesky(2); //Prior for correlation matrix for P
    //sigma_P[i] ~ cauchy(0, 2); //Prior for standard deviations for P
    sigma_P[i] ~ lognormal(0, 0.2); //Prior for standard deviations for P
    L_Omega_Sigma[i] ~ lkj_corr_cholesky(2); //Prior for correlation matrix for Sigma
    //sigma_Sigma[i] ~ cauchy(0, 2); //Prior for standard deviations for Sigma
    sigma_Sigma[i] ~ lognormal(0, 0.2); //Prior for standard deviations for Sigma
    theta_mats[i] ~ normal(0,1);
    log_lambdas[i] ~ normal(log(1),1.0);
    theta_mats[i] ~ normal(0,1);

  }


  // Measurement error
  for(n in 1:n_traits){
    y_obs[,n] ~ normal(y_true[,n], y_error[,n]);
    }

  target += calculate_log_likelihood_lp(
    n_traits,
    n_post_order_path_nodes,
    y_true,
    y_root,
    post_order_path_nodes,
    branch_lengths,
    branch_regime_idx, //Path of nodes from tip to root
    parent_of_node, //Path of nodes from tip to root
    node_types,
    H_mats,
    P_mats,
    inv_P_mats,
    Sigma_mats,
    theta_mats,
    lambdas_sum_mats);

}
generated quantities {
  //real g_i;
  //real sigma_ii;
  //real sigma_i;
  //real u_i;
  //vector[N] log_lik;
  //Based on https://cran.r-project.org/web/packages/loo/vignettes/loo2-non-factorized.htmlloo-cv-for-multivariate-normal-models
  //LOO-CV for multivariate normal models
  //matrix[N,N] V_total = V + diag_matrix(square(y_error));
  //matrix[N,N] inv_V = inverse(V_total);
  //for(i in 1:N){
      //g_i = (inv_V*(Y_obs-mu))[i];
      //sigma_ii = inv_V[i,i];
      //u_i = Y_obs[i]-g_i/sigma_ii;
      //sigma_i = 1/sigma_ii;
      //log_lik[i] = -0.5*log(2*pi()*sigma_i)-0.5*(square(Y_obs[i]-u_i)/sigma_i);
      //}
}
