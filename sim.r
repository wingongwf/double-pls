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
gamma_1 = 10  # TODO: need tuning
gamma_2 = 50  # TODO: need tuning

# Step 3: estimation of slopes
fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa, gamma_1=gamma_1, gamma_2=gamma_2)
A_hat = fit$A_hat
group_hat = fit$group_hat
regime_hat = fit$regime_hat
clog(paste0("group_hat=", vec2str(group_hat), ", regime_hat=", vec2str(regime_hat)))


