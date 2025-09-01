//08/16/2025 Programming mvOUOU model using Mitov et al. (2020) approach - working from
//08/23/2025 - Appears to be working version
//Single regime version

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

  matrix build_kron_sum(int k, matrix H) {
    int K2 = k * k;
    matrix[K2, K2] L;
    // initialize to zero
    for (i in 1:K2) for (j in 1:K2) L[i,j] = 0.0;

    // blockwise fill: block (a,b) of size k x k corresponds to:
    // L[( (a-1)*k+1 : a*k ), ( (b-1)*k+1 : b*k )]
    for (a in 1:k) {
      for (b in 1:k) {
        // position in big matrix
        int row0 = (a - 1) * k;
        int col0 = (b - 1) * k;
        // contribution from I ⊗ H: if a == b add H
        if (a == b) {
          for (r in 1:k) for (c in 1:k)
            L[row0 + r, col0 + c] += H[r, c];
        }
        // contribution from H ⊗ I: add H[a,b] * I_k
        real h_ab = H[a, b];
        if (h_ab != 0) {
          for (r in 1:k) L[row0 + r, col0 + r] += h_ab;
        }
      }
    }
    return L;
  }

  // compute V_i = integral_0^t exp(-H s) Sigma exp(-H' s) ds
  // using Lyapunov solution and the identity V_i = V_inf - exp(-H t) V_inf exp(-H' t)
  matrix compute_V_from_HSigma(int n_traits, real branch_length, matrix H, matrix Sigma) {
    int k = n_traits;
    // handle degenerate tiny branch lengths
    if (branch_length < 1e-6) {
      // extremely small -> return tiny SPD matrix
      matrix[k, k] Vi = 1e-12 * diag_matrix(rep_vector(1.0, k));
      return Vi;
    }

    // Build Kron-sum L = I ⊗ H + H ⊗ I
    matrix[k*k, k*k] Lbig = build_kron_sum(k, H);

    // vec(Sigma)
    vector[k*k] vecS;
    for (i in 1:k) for (j in 1:k) vecS[(i - 1) * k + j] = Sigma[i, j];

    // Add tiny jitter to diagonal of Lbig to avoid near-singularity
    for (d in 1:(k*k)) Lbig[d, d] += 1e-10;

    // Solve vec(V_inf) = inv(Lbig) * vecS
    // For small k, inverse() is OK; if k larger consider more stable solvers
    matrix[k*k, k*k] Linv = inverse(Lbig);
    vector[k*k] vecVinf = Linv * vecS;

    // reshape vecVinf -> V_inf
    matrix[k, k] Vinf;
    for (i in 1:k) for (j in 1:k) Vinf[i, j] = vecVinf[(i - 1) * k + j];

    // compute exp(-H t) and final V_i
    matrix[k, k] expHt = matrix_exp(-H * branch_length);
    matrix[k, k] Vi = Vinf - expHt * Vinf * expHt';

    // symmetrize + tiny jitter for numerical SPD
    Vi = 0.5 * (Vi + Vi');
    for (d in 1:k) Vi[d,d] += 1e-12;

    return Vi;
  }


  matrix build_G(int n_traits, int i, int j, real angle){
    matrix[n_traits,n_traits] G = identity_matrix(n_traits);
    real cos_angle = cos(angle);
    real sin_angle = sin(angle);
    G[i,i] = cos_angle;
    G[i,j] = -sin_angle;
    G[j,i] = sin_angle;
    G[j,j] = cos_angle;
    return G;
  }

  matrix build_Q(int n_traits, vector Givens_angles){ //Iterate through Givens_angles to update Q
    matrix[n_traits,n_traits] Q = identity_matrix(n_traits); //Start with identity matrix for Q
    int pos = 1;
    for(j in 2:n_traits){
      for(i in 1:(j-1)){
        Q = Q * build_G(n_traits, i, j, Givens_angles[pos]);
        pos = pos + 1;
      }
    }
    return Q;
  }

  matrix buld_lower_tri_T(int n_traits, vector T_lower_tri){
    matrix[n_traits,n_traits] T_lt = identity_matrix(n_traits); //Start with identity matrix for Q
    int pos = 1;
    for(j in 1:(n_traits-1)){ //Controls columns
      for(i in (j+1):n_traits){ //Controls rows
        T_lt[i,j] = T_lower_tri[pos];
        pos = pos + 1;
        }
      }
    return T_lt;
  }


  real calculate_log_likelihood_lp(
    int n_nodes,
    int n_post_order_path_nodes,
    int n_traits,
    matrix y_true,
    vector y_root,
    array[] int post_order_path_nodes,
    vector branch_lengths,
    array[] int branch_regime_idx, //Path of nodes from tip to root
    array[] int parent_of_node, //Path of nodes from tip to root
    array[] int node_types,
    matrix H_mats,
    //matrix P_mats,
    //matrix inv_P_mats,
    matrix Sigma_mats,
    vector theta_mats,
    matrix lambdas_sum_mats){
    array[n_nodes] matrix[n_traits,n_traits] L_kernals = rep_array(rep_matrix(0.0, n_traits, n_traits), n_post_order_path_nodes);
    array[n_nodes] vector[n_traits] m_kernals = rep_array(rep_vector(0.0,n_traits),n_post_order_path_nodes);
    array[n_nodes] real r_kernals = rep_array(0.0, n_post_order_path_nodes);

    for(i in 1:n_post_order_path_nodes){ //Loop through nodes in postorder
      int current_node = post_order_path_nodes[i]; //Get current node in postorder path
      int parent = parent_of_node[current_node]; //Get current node's parent
      real branch_length = branch_lengths[current_node]; //Get branch length leading to current node - from parent
      int current_node_type = node_types[current_node]; //Get current node type - 0 for root, 1 for tips, 2 for internal nodes
      matrix[n_traits,n_traits] Phi;
      vector[n_traits] omega;
      matrix[n_traits,n_traits] V_i;
      matrix[n_traits,n_traits] inv_V_i;

      matrix[n_traits,n_traits] A_i = rep_matrix(0.0, n_traits, n_traits);
      vector[n_traits] b_i = rep_vector(0.0, n_traits);
      matrix[n_traits,n_traits] C_i = rep_matrix(0.0, n_traits, n_traits);
      vector[n_traits] d_i = rep_vector(0.0, n_traits);
      matrix[n_traits,n_traits] E_i = rep_matrix(0.0, n_traits, n_traits);
      real f_i = 0;


      if(current_node_type != 0){ //Skip for root node
          Phi = compute_Phi(
          branch_length,
          H_mats);

        omega = compute_omega(
          n_traits,
          branch_length,
          H_mats,
          theta_mats);

        V_i = compute_V_from_HSigma(
          n_traits,
          branch_length,
          H_mats,
          Sigma_mats);


        //V_i = compute_V_i(
        //  n_traits,
        //  branch_length,
        //  P_mats,
        //  inv_P_mats,
        //  Sigma_mats,
        //  lambdas_sum_mats);

        inv_V_i = inverse(V_i);
        A_i = -0.5 * inv_V_i;
        b_i = inv_V_i* omega;
        C_i = -0.5 * Phi' *inv_V_i * Phi;
        d_i = -Phi' * inv_V_i * omega;
        E_i = Phi' * inv_V_i;
        f_i = -0.5 * omega' * inv_V_i * omega - (n_traits/2) * log(2*pi()) - 0.5 * log_determinant(V_i);
      }


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
          - 0.5 * log_determinant(-2*AL)
          - 0.25 * (b_i + m_i)' * inv_AL
          * (b_i + m_i);
      }
    }
    //Exiting loop at root, last in post_order_path
    int root_node = post_order_path_nodes[n_post_order_path_nodes];
    //return log of exp(y_root' * L_kernals[root_node] * y_root + y_root' * m_kernals[root_node] + r_kernals[root_node]);
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
  vector<lower=-pi()/2,upper=pi()/2>[n_traits*(n_traits-1)/2] Givens_angles; //Angles for G matrix, Givens rotation
  vector[n_traits * (n_traits - 1) / 2] T_lower_tri;
}

transformed parameters{
  matrix[n_traits, n_traits] lambdas_sum_mats;
  matrix[n_traits,n_traits] Sigma_mats = diag_matrix(sigma_Sigma) * Omega_Sigma * diag_matrix(sigma_Sigma);   //Original
  //matrix[n_traits,n_traits] H_mats = P_mats * diag_matrix(lambdas_mats) * inv_P_mats; //P * lambdas * P^-1 - Eigenvectors * eigenvales * inverse(Eigenvectors)

  for (i in 1:n_traits){//Calculate lambda sums matrix - lambda_i + lambda_j
    for (j in 1:n_traits) {
      lambdas_sum_mats[i, j] = lambdas_mats[i] + lambdas_mats[j];
    }
  }

  //Constructing Q
  matrix[n_traits,n_traits] Q = build_Q(n_traits, Givens_angles);

  //Constructing T
  matrix[n_traits,n_traits] T = buld_lower_tri_T(n_traits, T_lower_tri); //Build T and assign lower tri values

  for (i in 1:n_traits){
    T[i,i] = lambdas_mats[i]; //Assign diagonal of T
  }

  //Constructing H = Q * T * t(Q)
  matrix[n_traits,n_traits] H_mats = Q * T * Q';
  //complex_matrix[n_traits,n_traits] P = eigenvectors(H_mats);
  //matrix[n_traits,n_traits] P_mats = get_real(P);
  //matrix[n_traits,n_traits] inv_P_mats = inverse(P_mats);
}


model {
  //Priors
  for (i in 1:N) {
    for (k in 1:n_traits) {
      y_obs[i, k] ~ normal(y_true[i, k], y_error[i, k]);
    }
  }
  //for (k in 1:n_traits){ //Allowing y_true to be determined by OU model and observation data
  //  y_true[,k] ~ normal(0, 1.0);
  //}

  y_root ~ normal(0, 1.0);
  theta_mats ~ normal(0, 1.0);
  lambdas_mats ~ lognormal(1.0,0.5); //Normal distribution centered on 1, SD = 0.5; assume positive eigenvalues as in Mitov et al. 2019
  sigma_Sigma ~ lognormal(log(0.7), 0.25); // Centered around true sigma_Sigma values
  Omega_Sigma ~ lkj_corr(2); //4 = Peaked at 0 - extreme correlations less possible
  Givens_angles ~ uniform(-pi()/2,pi()/2);
  T_lower_tri ~ normal(0,0.5);

  target += calculate_log_likelihood_lp(
    n_nodes,
    n_post_order_path_nodes,
    n_traits,
    y_true,
    y_root,
    post_order_path_nodes,
    branch_lengths,
    branch_regime_idx,
    parent_of_node, //Path of nodes from tip to root
    node_types,
    H_mats,
    //P_mats,
    //inv_P_mats,
    Sigma_mats,
    theta_mats,
    lambdas_sum_mats);

}
generated quantities {
  }
