//mvBlouch mvOU - Mitov approach - fixed predictors model - 2025
//Multivariate OU model accounts for measurement error on X and/or Y variables
//After Bartoszek et al. 2012; Clavel et al. 2015; Mitov et al. 2019
//Single regime version,
functions {
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
        if (fabs(lambdas_sum) < 1e-8) {
          term1[i,j] = branch_length;}
        else {
          term1[i,j] = (1 - exp(-lambdas_sum * branch_length)) / lambdas_sum;
        }
      }
    }
    matrix[n_traits,n_traits] term2 = inv_P_branch  * (Sigma_branch * inv_P_branch');
    matrix[n_traits,n_traits] V_i = P_branch * (term1 .* term2) * P_branch';

    return V_i;
    }

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
      int current_node = post_order_path_nodes[i];
      int parent = parent_of_node[current_node];
      real branch_length = branch_lengths[current_node];
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


   // --- Debugging prints ---
        //print("Iter: ", i);

        //real det_V_i = determinant(V_i);
        //print("Node: ", current_node, " Branch length: ", branch_length, " Node type: ",
        //  current_node_type);
        //print("Node: ", current_node, " Branch length: ", branch_length, " Node type: ",
          //current_node_type, " Currrent regime: ", current_reg);

        //print("Determinant V_i: ", det_V_i);

        //if (det_V_i <= 0) {
        //  print("WARNING: V_i is not positive definite at node ", current_node);
        //}

        //vector[n_traits] eig_vals = eigenvalues_sym(V_i);
        //print("Eigenvalues V_i at node ", current_node, ": ", eig_vals);

      // Add jitter to ensure V_i is SPD

      // Diagnostic print statements
      //print("Checking V_i at node ", current_node, " (before jitter)...");
      //print("V_i[1,2]: ", V_i[1,2], ", V_i[2,1]: ", V_i[2,1]);
      //print("Eigenvalues of V_i: ", eigenvalues_sym(V_i));
      //if (fabs(V_i[1,2] - V_i[2,1]) > 1e-9) {
      //    print("Warning: V_i is not symmetric at node ", current_node);
      //}

      V_i = 0.5 * (V_i + V_i');  // force symmetry
      //V_i += 1e-3 * diag_matrix(rep_vector(1.0, n_traits)); // jitter
      V_i += 1e-6 * diag_matrix(rep_vector(1.0, n_traits));
      inv_V_i = mdivide_left_spd(V_i, identity_matrix(n_traits));
      //print(eigenvalues_sym(V_i));

      //inv_V_i = inverse(V_i);
      }

    //omega + Phi *x_j
    if (current_node_type == 1){// Is current node is a tip, use known values
      //print("Phi = ", Phi);
      //print("omega = ", omega);

      //print("Sigma_mats = ", Sigma_mats);
      //print("V_i = ", V_i);
      //print("eigenvalues of V_i = ", eigenvalues_sym(V_i));

      vector[n_traits] x_i = y_true[current_node,]'; //Observed tip value
      vector[n_traits] resid = x_i - omega; //

      //real logdet_term = log_determinant(2 * pi() * V_i);
      //print("Log determinant term at tip node ", current_node, ": ", logdet_term);

      C_i = -0.5 * Phi' * inv_V_i * Phi; //Quadratic coefficient matrix
      d_i = -Phi' * inv_V_i * resid; //Linear coefficient vector
      f_i = -0.5 * dot_product(resid, inv_V_i* resid) - //0.5 * log_determinant(2 * pi() * V_i);
        0.5 * n_traits * log(2 * pi()) - 0.5 * log_determinant(V_i);
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

      // Diagnostic print statements
      //print("Checking Q_i at node ", current_node, " (before jitter)...");
      //print("Q_i[1,2]: ", Q_i[1,2], ", Q_i[2,1]: ", Q_i[2,1]);
      //print("Eigenvalues of Q_i: ", eigenvalues_sym(Q_i));
      //if (fabs(Q_i[1,2] - Q_i[2,1]) > 1e-9) {
      //    print("Warning: Q_i is not symmetric at node ", current_node);
      //}

      matrix[n_traits,n_traits] inv_neg_Q_i = mdivide_left_spd(-Q_i, identity_matrix(n_traits));
      matrix[n_traits,n_traits] inv_Q_i = -inv_neg_Q_i; // Because inv(Q_i) = -inv(-Q_i)

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
    //print("i=", i,
    //  " node=", current_node,
    //  " parent=", parent,
    //  " branch_length=", branch_length,
    //  " Phi=", Phi,
    //  " Omega=", omega);

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

  array[n_post_order_path_nodes] int post_order_path_nodes; //Path of nodes from tip to root
  vector[n_nodes] branch_lengths;
  array[n_nodes] int branch_regime_idx; //Path of nodes from tip to root
  array[n_nodes] int parent_of_node; //Path of nodes from tip to root
  array[n_nodes] int node_types;
}
parameters {
  matrix[N,n_traits] y_true;
  vector[n_traits] y_root;
  //cholesky_factor_corr[n_traits] L_Omega_P;     // Correlation matrix for P (part of LKJ prior)
  //vector<lower=0>[n_traits] sigma_P; // Standard deviations for P (part of LKJ prior)
  //cholesky_factor_corr[n_traits] L_Omega_Sigma;     // Correlation matrix for Sigma (part of LKJ prior)
  //vector [n_traits] log_sigma_Sigma; // log-scale standard deviations for Sigma (part of LKJ prior) - can have negative log values
  //real<lower=-0.95, upper=0.95> rho;  // prevent singularities
  vector[n_traits] theta_mats;
  //vector[n_traits] theta_raw_z;
  //real<lower=0>theta_sd;
  //vector[n_traits] theta_mu;
}
transformed parameters{
  vector[n_traits] lambdas_mats;
  matrix[n_traits, n_traits] lambdas_sum_mats;
  //vector[n_traits] theta_mats;
  // 1. Convert log-lambdas to lambdas
  //vector[n_traits] log_lambdas = log_lambdas_mu + log_lambdas_sd * log_lambdas_raw_z;

  //lambdas_mats = exp(log_lambdas);
  lambdas_mats[1] = 2; //Fixed lambdas at known values
  lambdas_mats[2] = 2;

  // 2. Construct P matrix from LKJ prior components
  //vector[n_traits] sigma_P = exp(sigma_P_mu_log + sigma_P_sd_log * sigma_P_raw_z); //Non-centered

  //matrix[n_traits, n_traits] L_P = diag_pre_multiply(sigma_P, L_Omega_P); //diag_matrix(sigma_P) * L_Omega_P
  //P_mats = L * L';
  //P_mats = multiply_lower_tri_self_transpose(L_P);  // P = L * L'
  matrix[n_traits,n_traits] P_mats = diag_matrix(rep_vector(1.0, n_traits)); //Fix P at identity matrix

  // 3. Invert P
  matrix[n_traits,n_traits] inv_P_mats = mdivide_left_spd(P_mats, diag_matrix(rep_vector(1.0,n_traits)));

  // 4. Construct H = P * Lambda * P⁻¹
  matrix[n_traits,n_traits] H_mats = P_mats * diag_matrix(lambdas_mats) * inv_P_mats; // = H = [2,0,0,2]

  // 5. Construct Σ = diag(sigma_Sigma) * Omega_Sigma * diag(sigma_Sigma)
  //Fixed sigma_Sigma version, estimate Omega_Sigma
  //sigma_Sigma
  //vector<lower=0>[n_traits] sigma_Sigma; //
  //sigma_Sigma[1] = 0.5; //Standard deviations
  //sigma_Sigma[2] = 1.0; //Standard deviations

  //Omega_Sigma
  //corr_matrix[n_traits] Omega_Sigma;
  //matrix[2,2] Omega_Sigma;
  //Omega_Sigma[1,1] = 1;
  //Omega_Sigma[1,2] = 0.2;
  //Omega_Sigma[2,1] = 0.2;
  //Omega_Sigma[2,2] = 1;

  //matrix[n_traits,n_traits] Sigma_mats = diag_matrix(sigma_Sigma) * Omega_Sigma * diag_matrix(sigma_Sigma);   //Original
  //vector[n_traits] sigma_Sigma = exp(log_sigma_Sigma); //Log-scale for SD
  //vector[n_traits] sigma_Sigma = sigma_Sigma_mu + sigma_Sigma_sd * sigma_Sigma_raw_z; //Non-centered, non-log scale
  //vector[n_traits] sigma_Sigma = 0.5 + 0.5 * sigma_Sigma_raw_z; //Non-centered, non-log scale

  //vector[n_traits] sigma_Sigma = exp(sigma_Sigma_mu_log + sigma_Sigma_sd_log * sigma_Sigma_raw_z); //Non-centered
  //matrix[n_traits, n_traits] Omega_Sigma_i = multiply_lower_tri_self_transpose(L_Omega_Sigma); //Omega_Sigma_i = correlation matrix


  //matrix[n_traits,n_traits] Sigma_mats = diag_pre_multiply(sigma_Sigma, L_Omega_Sigma) * Sigma_Z;
  //matrix[n_traits,n_traits] Sigma_mats = quad_form_diag(Omega_Sigma_i, sigma_Sigma);
  //Using rho for diagonals
  //matrix[2,2] Omega_Sigma;
  //matrix[2,2] Sigma_mats;

  //Omega_Sigma[1,1] = 1;
  //Omega_Sigma[2,2] = 1;
  //Omega_Sigma[1,2] = rho;
  //Omega_Sigma[2,1] = rho;

  matrix[n_traits, n_traits] Sigma_mats;
  Sigma_mats[1,1] = 0.5;
  Sigma_mats[2,1] = 0.1414214;
  Sigma_mats[1,2] = 0.1414214;
  Sigma_mats[2,2] = 1.0;
  //Sigma_mats = diag_pre_multiply(sigma_Sigma, L_Omega_Sigma)
  //        * diag_pre_multiply(sigma_Sigma, L_Omega_Sigma)';

  //Sigma_mats = quad_form_diag(Omega_Sigma, sigma_Sigma);

  //L_Omega_Sigma
  //matrix[n_traits, n_traits] L_Omega_Sigma = diag_matrix(rep_vector(1.0, n_traits));//Identity matrix

  //matrix[n_traits, n_traits] Sigma_mats = diag_pre_multiply(sigma_Sigma, L_Omega_Sigma) *
  //                                      diag_pre_multiply(sigma_Sigma, L_Omega_Sigma)';

  //Sigma_mats = 0.5 * (Sigma_mats + Sigma_mats');

  //matrix[n_traits, n_traits] Sigma_mats; //Fix Sigma_mats at known values
  //Sigma_mats[1,1] = 0.1;
  //Sigma_mats[2,2] = 0.25;
  //Sigma_mats[1,2] = 0;
  //Sigma_mats[2,1] = 0;


  // 6. Precompute λ_i + λ_j for computing V_i
  for (k in 1:n_traits){
    for (l in 1:n_traits) {
    lambdas_sum_mats[k, l] = lambdas_mats[k] + lambdas_mats[l];
    }
  }

  // 7. Non-centered theta_mats
  //theta_mats = theta_mu +  theta_raw_z * theta_sd;
  //theta_mats[1] = 3;
  //theta_mats[2] = 8;
  //theta_mats[1] = 0.8978490; //Scaled versions based on scaled y
  //theta_mats[2] = 0.9064088;
}

model {
  //Priors
  for (i in 1:N) {
    for (k in 1:n_traits) {
      y_obs[i, k] ~ normal(y_true[i, k], y_error[i, k]);
    }
  }
  for(k in 1:n_traits){
    y_true[,k] ~ normal(0, 0.5);
  }

  y_root ~ normal(0, 1);

  //sigma_P ~ normal(0.5, 0.25);             // tighten: default was too wide

  //Log-scale Sigma SD
  //log_sigma_Sigma ~ normal(0,1);

  //Non-centered P
  //L_Omega_P ~ lkj_corr_cholesky(2); //Prior for correlation matrix for P
  //sigma_P_raw_z ~ std_normal();
  //sigma_P_mu_log ~ normal(log(1),0.2);
  //sigma_P_sd_log ~ normal(0,0.1);

  //Non-centered Sigmas
  //rho ~ normal(0, 1);
  //L_Omega_Sigma ~ lkj_corr_cholesky(2); //Prior for correlation matrix for Sigma
  //sigma_Sigma ~ normal(0, 1.0) T[0, ];
  //rho ~ uniform(-0.95, 0.95);                // or normal(0, 0.5) for stronger regularization
  //rho ~ normal(0, 0.5);                // or normal(0, 0.5) for stronger regularization
  //sigma_Sigma ~ lognormal(log(0.8), 0.2); // Centered around true sigma_Sigma values
  //sigma_Sigma ~ lognormal(0, 0.25); // Centered around true sigma_Sigma values
  //Omega_Sigma ~ lkj_corr(2); //4 = Peaked at 0 - extreme correlations less possible
  //to_vector(Sigma_Z) ~ std_normal();

  //sigma_Sigma_mu ~ normal(0,0.2);
  //sigma_Sigma_sd ~ normal(0,0.2) T[0 , ];

  //Non-centered log-lambdas
  //log_lambdas_raw_z ~ normal(0, 0.3);//std_normal();
  //log_lambdas_mu ~ normal(log(1),0.2);
  //log_lambdas_sd ~ normal(0,0.1);

  //Thetas
  theta_mats ~ normal(0, 1.0); // Strong regularization around 0
  //Non-centered thetas
  //theta_sd ~ normal(0,0.1);
  //theta_mu ~ normal(4,2);
  //to_vector(theta_raw_z) ~ std_normal();

  //log_lambdas ~ normal(log(2), 0.25); //Centered at 0.69 - centered

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
  //corr_matrix[n_traits] Omega_Sigma = L_Omega_Sigma * L_Omega_Sigma';
  vector[n_traits] Sigma_eigenvalues = eigenvalues_sym(Sigma_mats);


}
