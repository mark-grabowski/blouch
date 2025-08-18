//mvBlouch mvOU - Mitov approach - fixed predictors model - 2025
//Multivariate OU model accounts for measurement error on X and/or Y variables
//After Bartoszek et al. 2012; Clavel et al. 2015; Mitov et al. 2019
//Single regime version,

functions {
  // --- numerics helpers ------------------------------------------------------

  matrix spd_regularize(matrix A, real jitter) {
    int k = rows(A);
    return 0.5 * (A + A') + jitter * diag_matrix(rep_vector(1.0, k));
  }

  // Invert SPD using its Cholesky factor: inv(A) = L^{-T} L^{-1}
  matrix inv_from_chol(matrix L) {
    int k = rows(L);
    matrix[k,k] I = diag_matrix(rep_vector(1.0, k));
    matrix[k,k] Linv = mdivide_left_tri_low(L, I);
    return Linv' * Linv;
  }

  // log|A| from its Cholesky diag
  real logdet_from_chol(matrix L) {
    vector[rows(L)] d = diagonal(L);
    return 2 * sum(log(d));
  }

  // --- OU pieces -------------------------------------------------------------

  matrix compute_Phi(real t, matrix H) {
    return matrix_exp(-H * t);
  }

  vector compute_omega(int k, real t, matrix H, vector theta) {
    matrix[k,k] I = diag_matrix(rep_vector(1.0, k));
    return (I - matrix_exp(-H * t)) * theta;
  }

  // V_i via Mitov (Hadamard fraction with λ_i + λ_j) and P/Σ decomposition
  matrix compute_V_i(
      int k,
      real t,
      matrix P, matrix invP,
      matrix Sigma,                 // process covariance (k x k)
      matrix lambdas_sum) {         // precomputed (λ_i + λ_j)
    matrix[k,k] term1;
    for (i in 1:k) {
      for (j in 1:k) {
        real s = lambdas_sum[i, j];
        term1[i, j] = (fabs(s) < 1e-12) ? t : (1 - exp(-s * t)) / s;
      }
    }
    // V = P * ((term1 .* (invP * Sigma * invP'))) * P'
    matrix[k,k] term2 = invP * (Sigma * invP');
    matrix[k,k] V = P * (term1 .* term2) * P';
    return 0.5 * (V + V'); // symmetrize
  }

  // --- Pruning (Mitov Theorem 1), canonical kernel form ----------------------
  // We store at each parent a kernel for its x_parent:
  //   log φ(x) = -0.5 x' L x + m' x + r
  // but we *store* C := -0.5 L (for stability), so L = -2C.

  real calculate_log_likelihood_lp(
    int k,
    int n_post,
    matrix y_true,              // rows indexed by node id; tips are 1..N
    vector y_root,              // parameter
    array[] int post_order_nodes,
    vector branch_lengths,      // indexed by child node id
    array[] int parent_of,
    array[] int node_types,     // 0=root, 1=tip, 2=internal
    matrix H,                   // single-regime H
    matrix P,
    matrix invP,
    matrix Sigma,
    vector theta,
    matrix lambdas_sum
  ) {
    // storage sized by *max node id* so we can accumulate by node id
    int max_id = max(post_order_nodes);
    array[max_id] matrix[k,k] C_store;
    array[max_id] vector[k]   m_store;
    array[max_id] real        r_store;

    for (nid in 1:max_id) {
      C_store[nid] = rep_matrix(0.0, k, k);
      m_store[nid] = rep_vector(0.0, k);
      r_store[nid] = 0.0;
    }

    // constants
    real jitter = 1e-9;
    matrix[k,k] I = diag_matrix(rep_vector(1.0, k));

    // postorder: children before parents
    for (idx in 1:n_post) {
      int node  = post_order_nodes[idx];
      int p     = parent_of[node];
      int ntype = node_types[node];

      if (ntype != 0) { // skip the root (no branch above it)
        real t = branch_lengths[node];

        // edge matrices
        matrix[k,k] Phi;
        vector[k]   omega;
        matrix[k,k] V;

        if (t < 1e-12) {
          Phi   = I;
          omega = rep_vector(0.0, k);
          V     = 1e-12 * I;
        } else {
          Phi   = compute_Phi(t, H);
          omega = compute_omega(k, t, H, theta);
          V     = compute_V_i(k, t, P, invP, Sigma, lambdas_sum);
        }

        // regularize & invert V
        matrix[k,k] V_reg = spd_regularize(V, jitter);
        matrix[k,k] L_V   = cholesky_decompose(V_reg);
        matrix[k,k] Vinv  = inv_from_chol(L_V);
        real logdet_V     = logdet_from_chol(L_V);

        if (ntype == 1) {
          // TIP: x_tip ~ N(omega + Phi * x_parent, V)
          vector[k] x_tip = y_true[node, ]';
          vector[k] resid = x_tip - omega;                 // <-- this is Mitov’s (x_i - ω)
          matrix[k,k] L_tip = Phi' * Vinv * Phi;          // L_tip = Φ' V^{-1} Φ

          // convert to stored canonical form at the **parent**
          matrix[k,k] C_tip = -0.5 * L_tip;               // C = -0.5 L
          vector[k]   m_tip = (Phi' * Vinv) * resid;      // m = Φ' V^{-1} (x_tip - ω)
          real        r_tip = -0.5 * dot_product(resid, Vinv * resid)
                              - 0.5 * k * log(2 * pi())
                              - logdet_V / 2;

          C_store[p] += C_tip;
          m_store[p] += m_tip;
          r_store[p] += r_tip;
        }

        if (ntype == 2) {
          // INTERNAL: combine child-kernel-at-node with this edge; lift to parent
          matrix[k,k] L_child = -2.0 * C_store[node];     // L_child from stored C
          vector[k]   m_child = m_store[node];
          real        r_child = r_store[node];

          // edge quadratic form pieces
          matrix[k,k] A = -0.5 * Vinv;                    // A = -1/2 V^{-1}
          vector[k]   b = Vinv * omega;                   // b = V^{-1} ω
          matrix[k,k] Q = 0.5 * (A + L_child + (A + L_child)');  // symmetrize

          // We need (-2Q) SPD for the theorem’s log| -2Q |
          matrix[k,k] M = -2.0 * Q;
          matrix[k,k] M_reg = spd_regularize(M, jitter);
          matrix[k,k] L_M   = cholesky_decompose(M_reg);
          matrix[k,k] Minv  = inv_from_chol(L_M);         // Minv = (-2Q)^{-1}
          real logdet_M     = logdet_from_chol(L_M);
          matrix[k,k] Qinv  = -0.5 * Minv;                // Q^{-1} = -1/2 * (-2Q)^{-1}

          vector[k] z = b + m_child;

          matrix[k,k] C_up = -0.5 * Phi' * Qinv * Phi;
          vector[k]   m_up = (Phi' * Qinv) * z;
          real        r_up = r_child
                             - 0.25 * dot_product(z, Qinv * z)
                             + 0.5 * k * log(2 * pi())
                             - 0.5 * logdet_M;

          C_store[p] += C_up;
          m_store[p] += m_up;
          r_store[p] += r_up;
        }
      }
    }

    // evaluate root kernel at y_root
    int root = post_order_nodes[n_post];
    matrix[k,k] C_root = C_store[root];
    vector[k]   m_root = m_store[root];
    real        r_root = r_store[root];

    // -0.5 * y' L y + m' y + r  with L = -2C
    return -0.5 * quad_form_sym(-2.0 * C_root, y_root)
           + dot_product(m_root, y_root)
           + r_root;
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
  vector[n_traits] theta_mats;
}
transformed parameters{
  vector[n_traits] lambdas_mats;
  matrix[n_traits, n_traits] lambdas_sum_mats;

  lambdas_mats[1] = 2; //Fixed lambdas at known values
  lambdas_mats[2] = 2;

  matrix[n_traits,n_traits] P_mats = diag_matrix(rep_vector(1.0, n_traits)); //Fix P at identity matrix
  matrix[n_traits,n_traits] inv_P_mats = mdivide_left_spd(P_mats, diag_matrix(rep_vector(1.0,n_traits)));
  matrix[n_traits,n_traits] H_mats = P_mats * diag_matrix(lambdas_mats) * inv_P_mats; // = H = [2,0,0,2]

matrix[n_traits, n_traits] Sigma_mats;
Sigma_mats[1,1] = 0.5;
Sigma_mats[1,2] = 0.1414214;
Sigma_mats[2,1] = 0.1414214;
Sigma_mats[2,2] = 1.0;
// (Optionally) enforce symmetry:
Sigma_mats = 0.5 * (Sigma_mats + Sigma_mats');

  // 6. Precompute λ_i + λ_j for computing V_i
  for (k in 1:n_traits){
    for (l in 1:n_traits) {
      lambdas_sum_mats[k, l] = lambdas_mats[k] + lambdas_mats[l];
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
    y_true[,k] ~ normal(0, 0.5);
  }

  y_root ~ normal(0, 1);

  theta_mats ~ normal(0, 1.0); // regularize thetas

target += calculate_log_likelihood_lp(
  n_traits,
  n_post_order_path_nodes,
  y_true,
  y_root,
  post_order_path_nodes,
  branch_lengths,
  parent_of_node,
  node_types,
  H_mats,
  P_mats,
  inv_P_mats,
  Sigma_mats,
  theta_mats,
  lambdas_sum_mats
);

}
generated quantities {
}
