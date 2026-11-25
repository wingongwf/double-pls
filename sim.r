# This module consists of all functions for running simulations
source("dgp.r")
source("est.r")
source("util.r")

# Set up parameters
N = 10     # TODO: small parameters test case
T_ = 20    # TODO: small parameters test case
r = 2  
p = 2      # TODO: small parameters test case
phi = 0.8
K0 = 2
m0 = 1
pi_ = 0.2
L = 5
kappa = 2

# file names
LOG_DATETIME = format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE = paste0("logs/simulation_", LOG_DATETIME, ".log")
CSV_BETA = paste0("results/Btilde_", LOG_DATETIME, ".csv")

set.seed(1233)
make_csv_Btilde(N=N, T_=T_, p=p, filename=CSV_BETA)

# Step 1: generate observations X and y
# X: N x T_ x p
# y: N x T_
data = generate_data(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, m0=m0, pi_=pi_, L=L)
X = data$X
y = data$y
B0 = data$B
assert_dimensions(X, dim_=c(N, T_, p))
assert_dimensions(y, dim_=c(N, T_))

# Step 2: tune gamma_1 and gamma_2
gamma_max = c(100, 100)  # TODO: need tuning
gamma_grid = get_gamma_grid(gamma_max=gamma_max, gamma_min_ratio=1e-4, num_gamma=5)
gamma_hat = tune_gamma(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa, gamma_grid=gamma_grid)

# Step 3: estimate slopes with given gamma_1 and gamma_2
# A_hat: K_hat x (m_hat+1) x p
fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa, gamma_1=gamma_hat[1], gamma_2=gamma_hat[2])
A_hat = fit$A_hat
group_hat = fit$group_hat
regime_hat = fit$regime_hat
clog(paste0("group_hat=", vec2str(group_hat), ", regime_hat=", vec2str(regime_hat)))

# Step 4: reshape true slopes (B0) and estimated slopes (B_tilde) to vectors
# B0_vec: (N*T_*p) x 1
# B_tilde_vec: (N*T_*p) x 1
B0_vec = B02vec(B0=B0, N=N, T_=T_, p=p)
B_tilde_vec = Btilde2vec(A_hat=A_hat, group_hat=group_hat, regime_hat=regime_hat, N=N, T_=T_, p=p)

# Export B_tilde_vec to csv
csv_line = paste0(paste(B_tilde_vec, collapse=","), "\n")
cat(csv_line, file=CSV_BETA, append=TRUE)




