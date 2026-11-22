# This module consists of all functions for running simulations
source("dgp.r")
source("est.r")
source("util.r")

# Set up parameters
N = 30
T_ = 50
r = 2
p = 4
phi = 0.8
K0 = 2
m0 = 1
pi_ = 0.2
L = 5
kappa = 2

# file names
LOG_DATETIME = format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE = paste0("logs/simulation_", LOG_DATETIME, ".log")

set.seed(1233)

# Step 1: generate observations X and y
# X: N x T_ x p
# y: N x T_
data = generate_data(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, m0=m0, pi_=pi_, L=L)
X = data$X
y = data$y
assert_dimensions(X, dim_=c(N, T_, p))
assert_dimensions(y, dim_=c(N, T_))

# Step 2: fine-tune gamma_1 and gamma_2
gamma_1 = 1
gamma_2 = 1

gamma_1 = NULL
gamma_2 = NULL

# Step 3: estimate the slopes B with given gamma_1 and gamma_2
M_F_hat = compute_M_F_hat(N=N, T_=T_, p=p, X=X)
B_hat = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, M_F_hat=M_F_hat, kappa=kappa)

# Step 4: identify the structure of slopes: groups and regimes
group_hat = identify_groups(N=N, T_=T_, p=p, B=B_hat)
regime_hat = identify_regimes(N=N, T_=T_, p=p, B=B_hat)
K_hat = max(group_hat)
m_hat = max(regime_hat) - 1
clog(paste0("group=", vec2str(group_hat), ", regime=", vec2str(regime_hat), ", K_hat=", K_hat, ", m_hat=", m_hat))

# Step5: post-selection estimation of slopes
# A_hat: K_hat x (m_hat+1) x p
A_hat = estimate_slopes_post_selection(
  p=p, X=X, y=y, M_F_hat=M_F_hat, group_hat=group_hat, regime_hat=regime_hat, K_hat=K_hat, m_hat=m_hat)
assert_dimensions(A_hat, dim_=c(K_hat, m_hat+1, p))



