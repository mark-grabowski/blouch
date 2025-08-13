# Test 1: Compare Mitov (Theorem 1) formulas to your Stan algebra (tip + one parent)
library(Matrix)    # for symmetric matrices etc

# helper: logdet that returns numeric
logdet_num <- function(M){
  # Use eigen to avoid sign issues
  ev <- eigen(M, symmetric = TRUE)$values
  if(any(ev <= 0)) return(NA_real_)
  return(sum(log(ev)))
}

# 1) Set up a small example (2 traits)
k <- 2

# choose H diagonal (lambda), use scaling such that Phi = exp(-H * t)
H <- diag(c(2, 2))   # as in your simulation
t_branch <- 0.3

Phi <- expm::expm(-H * t_branch)   # compute matrix exponential
# choose theta (optima)
theta <- c(3, 8)

# choose Sigma (process scale) in the Mitov formulation:
sigma_Sigma <- c(0.5, 1.0)        # SDs
rho <- 0.2                        # correlation
Omega <- matrix(c(1, rho, rho, 1), 2, 2)
Sigma <- diag(sigma_Sigma) %*% Omega %*% diag(sigma_Sigma) # full Sigma
# Note: Mitov sometimes uses Sigma_sq etc; here Sigma = process noise covariance

# Compute omega (mean shift) per your function: omega = (I - exp(-H t)) * theta
I <- diag(k)
omega <- (I - Phi) %*% theta

# Compute V_i the branch variance using Mitov's closed form for simple diagonal H
# Here we compute V_stationary via solving Lyapunov: vec(Vs) = solve(L) %*% vec(Sigma*Sigma')
# Simpler: for the 2x2 diagonal H, we can just use the formula implemented in your Stan compute_V_i.
# We'll compute it the same as in Stan: term1 and term2 approach.
lambdas <- diag(H)
lambdas_sum <- outer(lambdas, lambdas, "+")
term1 <- matrix(0, k, k)
for(i in 1:k) for(j in 1:k){
  s <- lambdas_sum[i,j]
  if (abs(s) < 1e-12) term1[i,j] <- t_branch
  else term1[i,j] <- (1 - exp(-s * t_branch)) / s
}
# Here P = I (you used P_mats = I in many tests); inv_P = I
P <- diag(k)
invP <- diag(k)
term2 <- invP %*% (Sigma %*% t(invP))
V <- P %*% (term1 * term2) %*% t(P)
# Force symmetry
V <- 0.5*(V + t(V))

# Check PD
cat("Eigenvalues V:", eigen(V, symmetric=TRUE)$values, "\n")
if(any(eigen(V, symmetric=TRUE)$values <= 0)) stop("V not PD")

invV <- solve(V)

# Now take an observed tip x_i (simulate or choose)
# We'll choose x_i consistent with the model: sample from N(omega + Phi * parent_state, V)
# choose parent_state (x_j), say root at zero for simplicity
x_parent <- c(0, 0)
x_mean <- as.numeric(omega + Phi %*% x_parent)
set.seed(123)
x_i <- as.numeric(MASS::mvrnorm(1, mu = x_mean, Sigma = V))

cat("Phi:\n"); print(Phi)
cat("omega:\n"); print(omega)
cat("V:\n"); print(V)
cat("x_i (sim):\n"); print(x_i)
cat("x_parent (chosen):\n"); print(x_parent)
cat("x_mean:\n"); print(x_mean)

# ----------------------------
# Theorem 1 explicit terms for node i (tip)
# Ai = -1/2 V^{-1}
Ai <- -0.5 * invV

# bi = V^{-1} * omega
bi <- invV %*% omega

# Ci = -1/2 Phi' V^{-1} Phi
Ci <- -0.5 * t(Phi) %*% invV %*% Phi

# di = -Phi' V^{-1} omega
di <- - t(Phi) %*% invV %*% omega

# Ei = Phi' V^{-1}
Ei <- t(Phi) %*% invV

# fi = -1/2 omega' V^{-1} omega - k/2 log(2π) - 1/2 log|V|
fi <- -0.5 * t(omega) %*% invV %*% omega - (k/2)*log(2*pi) - 0.5 * logdet_num(V)

cat("\n--- Theorem 1 quantities (Ai, bi, Ci, di, fi) ---\n")
print(Ai); print(bi); print(Ci); print(di); print(fi)

# ----------------------------
# Your Stan tip computations (as in your code)
resid <- x_i - as.numeric(omega)   # x_i - omega
C_stan_tip <- -0.5 * t(Phi) %*% invV %*% Phi
d_stan_tip <- - t(Phi) %*% invV %*% resid
f_stan_tip <- -0.5 * as.numeric(t(resid) %*% invV %*% resid) - 0.5 * k * log(2*pi) - 0.5 * logdet_num(V)

cat("\n--- Stan-style tip quantities (C, d, f) ---\n")
print(C_stan_tip); print(d_stan_tip); print(f_stan_tip)

# Now show algebraic relation between theorem di/fi and stan d/f:
# Expand theorem 'di' and 'fi' to show they match the way Stan splits terms:
# stan d = -Phi' V^{-1} (x_i - omega) = -Phi' V^{-1} x_i + Phi' V^{-1} omega
# So stan d = (-Phi' V^{-1} x_i) + (+ something). Mitov's di = -Phi' V^{-1} omega, so difference is -Phi' V^{-1} x_i
cat("\nCheck difference between Stan d_stan_tip and ( -Phi' V^{-1} x_i + Mitov_di ):\n")
lhs <- as.numeric(d_stan_tip)
rhs <- as.numeric(- t(Phi) %*% invV %*% x_i + di)
print(cbind(lhs=lhs, rhs=rhs, diff=lhs-rhs))

# fi comparison: Stan uses full quadratic in resid; theorem uses -1/2 omega' invV omega + terms.
cat("\nCheck fi equality (Stan f vs theorem f adjusted for observed x_i term):\n")
# Expand the quadratic: -0.5 (x - omega)' V^{-1} (x - omega) = -0.5 x' V^{-1} x + x' V^{-1} omega - 0.5 omega' V^{-1} omega
# So Stan's f = -0.5 x'V^{-1}x + x'V^{-1}omega -0.5 omega'V^{-1}omega - k/2 log(2pi) - 0.5 log |V|
# Theorem fi = -0.5 omega'V^{-1}omega - k/2 log(2pi) - 0.5 log |V|
# So Stan f = theorem fi + (-0.5 x'V^{-1}x + x'V^{-1}omega)
lhs_f <- as.numeric(f_stan_tip)
rhs_f <- as.numeric(fi + (-0.5 * t(x_i) %*% invV %*% x_i + t(x_i) %*% invV %*% omega))
print(c(lhs_f, rhs_f, lhs_f - rhs_f))

# ----------------------------
# Now simulate combining a single child contribution into parent:
# In Stan: you accumulate child's L (call it L_child = Ci_child), m_child = d_child,
# and r_child = f_child, then for the internal node you do:
L_child <- C_stan_tip        # matrix
m_child <- d_stan_tip        # vector
r_child <- f_stan_tip        # scalar

# Your A_i = -0.5 invV, b_i = invV * omega (already above)
Q_i <- Ai + L_child          # note Ai is -0.5 invV; L_child is C from tip
# For numerical safety check -Q_i SPD
eig_negQ <- eigen(-Q_i, symmetric = TRUE)$values

cat("\nEigenvalues of -Q_i (should be positive):\n"); print(eig_negQ)

inv_negQ <- tryCatch(solve(-Q_i), error=function(e) { NULL })
if(is.null(inv_negQ)) stop("(-Q_i) not invertible in this example; check signs/jitter")

inv_Q  <- -inv_negQ
z_vec <- as.numeric(bi + m_child)   # b_i + m

# Stan internal update:
C_upd <- -0.5 * t(Phi) %*% inv_Q %*% Phi
d_upd <- - t(Phi) %*% inv_Q %*% z_vec
f_upd <- as.numeric(r_child - 0.25 * t(z_vec) %*% inv_Q %*% z_vec + 0.5 * log(2*pi) - 0.5 * logdet_num(-2 * Q_i))

cat("\n--- Internal-node (Stan completed-square) updates ---\n")
print(C_upd); print(d_upd); print(f_upd)

# For sanity, no theorem explicit expression for internal combined step printed here,
# but you can verify symmetry and smallness
cat("\nSymmetry checks (C_upd symmetric?):\n"); print(max(abs(C_upd - t(C_upd))))

# Done
