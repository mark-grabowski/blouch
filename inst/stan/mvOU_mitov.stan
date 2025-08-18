//08/16/2025 Programming mvOUOU model using Mitov et al. (2020) approach

functions {
   matrix compute_Phi(real branch_length, matrix H_mat){
    return matrix_exp(-H_mat * branch_length);
  }

  vector compute_omega(int n_traits, real branch_length, matrix H_mat, vector theta_mat){
    matrix[n_traits,n_traits] I = identity_matrix(n_traits);
    return ((I - matrix_exp(-H_mat * branch_length)) * theta_mat);
  }

  matrix compute_V_i(int n_traits, real branch_length, matrix P_branch, matrix inv_P_branch,  matrix Sigma_branch, matrix lambdas_sum_matrix_branch){
   // matrix[n_traits,n_traits] term1 = (1 - exp(-lambdas_sum_matrix_branch * branch_length)) / lambdas_sum_matrix_branch;
    matrix[n_traits,n_traits] term1;
    for(i in 1:n_traits){
      for(j in 1:n_traits){
        real lambdas_sum = lambdas_sum_matrix_branch[i,j];
        term1[i,j] = (1 - exp(-lambdas_sum * branch_length)) / lambdas_sum;
        }
    }
    matrix[n_traits,n_traits] term2 = inv_P_branch  * (Sigma_branch * inv_P_branch');
    matrix[n_traits,n_traits] V_i = P_branch * (term1 .* term2) * P_branch';

    return V_i;
    }

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
    matrix H_mats,
    matrix P_mats,
    matrix inv_P_mats,
    matrix Sigma_mats,
    vector theta_mats,
    matrix lambdas_sum_mats){
    array[n_post_order_path_nodes] matrix[n_traits,n_traits] L_kernals = rep_array(rep_matrix(0.0, n_traits, n_traits), n_post_order_path_nodes);
    array[n_post_order_path_nodes] vector[n_traits] m_kernals = rep_array(rep_vector(0.0,n_traits),n_post_order_path_nodes);
    array[n_post_order_path_nodes] real r_kernals = rep_array(0.0, n_post_order_path_nodes);

    for(i in 1:n_post_order_path_nodes){ //Loop through nodes in postorder
      int current_node = post_order_path_nodes[i]; //Get current node in postorder path
      int parent = parent_of_node[current_node]; //Get current node's parent
      real branch_length = branch_lengths[current_node]; //Get branch length leading to current node - from parent
      int current_node_type = node_types[current_node]; //Get current node type - 0 for root, 1 for tips, 2 for internal nodes
      matrix[n_traits,n_traits] Phi;
      vector[n_traits] omega;
      matrix[n_traits,n_traits] V_i;
      matrix[n_traits,n_traits] inv_V_i;

      if(current_node_type != 0){ //Skip for root node
          Phi = compute_Phi(
          branch_length,
          H_mats);

        omega = compute_omega(
          n_traits,
          branch_length,
          H_mats,
          theta_mats);

        V_i = compute_V_i(
          n_traits,
          branch_length,
          P_mats,
          inv_P_mats,
          Sigma_mats,
          lambdas_sum_mats);
      }
      inv_V_i = inverse(V_i);
      matrix[n_traits,n_traits] A_i = -0.5 * inv_V_i;
      vector[n_traits] b_i = inv_V_i* omega;
      matrix[n_traits,n_traits] C_i = -0.5 * Phi' *inv_V_i * Phi;
      vector[n_traits] d_i = -Phi' * inv_V_i * omega;
      matrix[n_traits,n_traits] E_i = Phi' * inv_V_i;
      real det_V = determinant(V_i);
      real f_i = -0.5 * omega' * inv_V_i * omega - (n_traits/2) * log(2*pi()) - 0.5 * log(det_V);

      if(current_node_type == 1){ //If type = tips
        vector[n_traits] x_i = y_true[current_node,]'; //Observed tip value
        L_kernals[parent] += C_i;
        m_kernals[parent] += d_i + E_i * x_i;
        r_kernals[parent] += x_i' * A_i * x_i + x_i' * b_i + f_i;
      }
      if(current_node_type == 2){ //If type = internals
        matrix[n_traits,n_traits] L_i = L_kernals[current_node];
        vector[n_traits] m_i = m_kernals[current_node];
        real r_i = r_kernals[current_node];
        matrix[n_traits,n_traits] AL = A_i + L_i;
        matrix[n_traits,n_traits] inv_AL = inverse(AL);
        L_kernals[parent] += C_i - 0.25 * E_i * inv_AL * E_i';
        m_kernals[parent] += d_i - 0.5 * E_i * inv_AL * (b_i + m_i);
        r_kernals[parent] += f_i + r_i + (n_traits/2) * log(2*pi())
          - 0.5 * log(determinant(-2*AL))
          - 0.25 * (b_i + m_i)' * inv_AL
          * (b_i + m_i);
      }
    }
    //Exiting loop at root, last in post_order_path
    int root_node = post_order_path_nodes[n_post_order_path_nodes];
    //return exp(y_root' * L_kernals[root_node] * y_root + y_root' * m_kernals[root_node] + r_kernals[root_node]);
    return (y_root' * L_kernals[root_node] * y_root + y_root' * m_kernals[root_node] + r_kernals[root_node]);
  }
}
data {
  int<lower=1> N;
  int<lower=1> n_nodes;
  int<lower=1>  n_traits;
  int<lower=1> n_regs;
  int<lower=1> n_post_order_path_nodes;
  int<lower=1> n_parent_of_nodes;

  matrix[N,n_traits] y_obs;
  matrix[N,n_traits] y_error;

  array[n_post_order_path_nodes] int post_order_path_nodes;
  vector[n_nodes] branch_lengths;
  array[n_nodes] int branch_regime_idx;
  array[n_nodes] int parent_of_node;
  array[n_nodes] int node_types;
}
parameters {
  matrix[N,n_traits] y_true;
  vector[n_traits] y_root;
  vector[n_traits] theta_mats;
  vector<lower=0.001>[n_traits] lambdas_mats; //Eigenvalues - must be > 0
  corr_matrix[n_traits] Omega_Sigma; //Correlation component of Sigma matrix
  vector<lower=0.001>[n_traits] sigma_Sigma; //Standard deviations for Sigma (part of LKJ prior)
  matrix[n_traits,n_traits] P_mats;

}
transformed parameters{
  matrix[n_traits, n_traits] lambdas_sum_mats;
  matrix[n_traits,n_traits] Sigma_mats = diag_matrix(sigma_Sigma) * Omega_Sigma * diag_matrix(sigma_Sigma);   //Original
  matrix[n_traits,n_traits] inv_P_mats = inverse(P_mats);
  matrix[n_traits,n_traits] H_mats = P_mats * diag_matrix(lambdas_mats) * inv_P_mats; //P * lambdas * P^-1 - Eigenvectors * eigenvales * inverse(Eigenvectors)


  // 6. Precompute λ_i + λ_j for computing V_i
  for (i in 1:n_traits){
    for (j in 1:n_traits) {
      lambdas_sum_mats[i, j] = lambdas_mats[i] + lambdas_mats[j];
    }
  }
}
model {
  //Priors
  for (i in 1:N) {
    for (k in 1:n_traits) {
      y_obs[i, k] ~ normal(y_true[i, k], y_error[i, k]);
    }
  }
  for (k in 1:n_traits){
    y_true[,k] ~ normal(0, 1.0);
  }
  y_root ~ normal(0, 1.0);
  theta_mats ~ normal(0, 1.0);
  lambdas_mats ~ lognormal(1.0,0.5); //Normal distribution centered on 1, SD = 0.5;
  sigma_Sigma ~ lognormal(1, 0.5); // Centered around true sigma_Sigma values
  Omega_Sigma ~ lkj_corr(4); //4 = Peaked at 0 - extreme correlations less possible


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
}
