# ==== R SCRIPT: Mitov blocks vs Stan-style pruning, node-by-node ====

library(MASS)
library(expm)

logdet_num <- function(M) {
  ev <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  sum(log(ev))
}

# ---- Model setup (2 traits) ----
k <- 2
I2 <- diag(2)

H <- diag(c(2, 2))        # OU rate matrix (diagonal for clarity)
theta <- c(3, 8)          # optima
sigma_Sigma <- c(0.5, 1)  # process SDs
rho <- 0.2
Omega <- matrix(c(1, rho, rho, 1), 2, 2)
Sigma <- diag(sigma_Sigma) %*% Omega %*% diag(sigma_Sigma)  # process covariance

# helper: Phi(t), omega(t), V(t)
Phi_of_t <- function(t) expm(-H * t)
omega_of_t <- function(t) (I2 - Phi_of_t(t)) %*% theta

# V(t): solve Lyapunov stationary, then V(t) = Vs - exp(-H t) Vs exp(-H' t)
V_of_t <- function(t) {
  # stationary covariance: solve (H Vs + Vs H') = Sigma
  # vectorized solution: vec(Vs) = solve(kron(I,H)+kron(t(H),I)) vec(Sigma)
  L <- kronecker(I2, H) + kronecker(t(H), I2)
  Vs_vec <- solve(L, as.vector(Sigma))
  Vs <- matrix(Vs_vec, k, k)
  Vt <- Vs - Phi_of_t(t) %*% Vs %*% t(Phi_of_t(t))
  0.5 * (Vt + t(Vt))
}

# Mitov blocks for a branch of length t
mitov_blocks <- function(t) {
  Phi <- Phi_of_t(t)
  omega <- omega_of_t(t)
  V <- V_of_t(t)
  V <- V + 1e-10 * I2
  invV <- solve(V)

  Ai <- -0.5 * invV
  bi <- invV %*% omega
  Ci <- -0.5 * t(Phi) %*% invV %*% Phi
  di <- - t(Phi) %*% invV %*% omega
  Ei <- t(Phi) %*% invV
  fi <- as.numeric(-0.5 * t(omega) %*% invV %*% omega -
                     (k/2) * log(2*pi) - 0.5 * logdet_num(V))
  list(Phi=Phi, omega=omega, V=V, invV=invV, Ai=Ai, bi=bi, Ci=Ci, di=di, Ei=Ei, fi=fi)
}

# ---- 1) Single tip under a parent ----
set.seed(123)
t_tip <- 0.3

blk <- mitov_blocks(t_tip)
Phi <- blk$Phi; omega <- blk$omega; V <- blk$V; invV <- blk$invV
Ai <- blk$Ai; bi <- blk$bi; Ci <- blk$Ci; di <- blk$di; fi <- blk$fi

x_parent <- c(0, 0)  # choose a parent state
x_mean <- as.numeric(omega + Phi %*% x_parent)
x_tip <- as.numeric(mvrnorm(1, mu = x_mean, Sigma = V))

cat("\n== TIP STEP: Mitov blocks vs Stan-style ==\n")
# Stan-style tip terms
resid <- x_tip - as.numeric(omega)
C_stan <- -0.5 * t(Phi) %*% invV %*% Phi
d_stan <- - t(Phi) %*% invV %*% resid
f_stan <- as.numeric(-0.5 * t(resid) %*% invV %*% resid -
                       (k/2) * log(2*pi) - 0.5 * logdet_num(V))

# Show equality:
# C: identical to Ci
cat("||C_stan - Ci||_max =", max(abs(C_stan - Ci)), "\n")

# d: d_stan  ==  (- Phi' invV x_tip + di)
d_mitov_full <- as.numeric(- t(Phi) %*% invV %*% x_tip + di)
cat("||d_stan - d_mitov_full||_max =", max(abs(d_stan - d_mitov_full)), "\n")

# f: f_stan == fi + [ -0.5 x'invV x + x'invV omega ]
rhs_f <- as.numeric(fi + (-0.5 * t(x_tip) %*% invV %*% x_tip + t(x_tip) %*% invV %*% omega))
cat("abs(f_stan - rhs_f) =", abs(f_stan - rhs_f), "\n")

# ---- 2) Two tips under one internal node: combine at internal ----
t_tip1 <- 0.30
t_tip2 <- 0.45
blk1 <- mitov_blocks(t_tip1)
blk2 <- mitov_blocks(t_tip2)

# simulate the two tips given the same internal parent state x_int
x_int <- c(0.4, -0.2)
x1_mean <- as.numeric(blk1$omega + blk1$Phi %*% x_int)
x2_mean <- as.numeric(blk2$omega + blk2$Phi %*% x_int)
x1 <- as.numeric(mvrnorm(1, x1_mean, blk1$V))
x2 <- as.numeric(mvrnorm(1, x2_mean, blk2$V))

# Tip kernels accumulated on the internal parent
resid1 <- x1 - as.numeric(blk1$omega)
resid2 <- x2 - as.numeric(blk2$omega)

C1 <- -0.5 * t(blk1$Phi) %*% blk1$invV %*% blk1$Phi
d1 <- - t(blk1$Phi) %*% blk1$invV %*% resid1
f1 <- as.numeric(-0.5 * t(resid1) %*% blk1$invV %*% resid1 -
                   (k/2) * log(2*pi) - 0.5 * logdet_num(blk1$V))

C2 <- -0.5 * t(blk2$Phi) %*% blk2$invV %*% blk2$Phi
d2 <- - t(blk2$Phi) %*% blk2$invV %*% resid2
f2 <- as.numeric(-0.5 * t(resid2) %*% blk2$invV %*% resid2 -
                   (k/2) * log(2*pi) - 0.5 * logdet_num(blk2$V))

L_child <- C1 + C2
m_child <- d1 + d2
r_child <- f1 + f2

# Branch from INTERNAL -> ROOT
t_int_to_root <- 0.25
blki <- mitov_blocks(t_int_to_root)
Ai <- blki$Ai; bi <- blki$bi; Phi_i <- blki$Phi

Q <- Ai + L_child                       # should be negative definite
inv_Q <- - solve(-Q)                    # numerically safer to invert -Q

z <- as.numeric(bi + m_child)

C_upd <- -0.5 * t(Phi_i) %*% inv_Q %*% Phi_i
d_upd <- - t(Phi_i) %*% inv_Q %*% z
f_upd <- as.numeric(r_child - 0.25 * t(z) %*% inv_Q %*% z +
                      0.5 * log(2*pi) - 0.5 * logdet_num(-2 * Q))

cat("\n== INTERNAL COMBINE (two children) ==\n")
cat("Eigenvalues(-Q) (should be >0):", eigen(-Q, symmetric=TRUE)$values, "\n")
cat("C_upd (symmetric?) max|C - C'| =", max(abs(C_upd - t(C_upd))), "\n")
cat("d_upd =", paste(round(d_upd, 6), collapse = ", "), "\n")
cat("f_upd =", round(f_upd, 6), "\n")

cat("\nAll equalities above should be ~machine-precision (1e-10-ish).\n")
