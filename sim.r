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
gamma_1_max = 100  # TODO: need tuning
gamma_2_max = 100  # TODO: need tuning

gamma_min_ratio = 1e-4
num_gamma = 5
gamma_1_min = gamma_1_max * gamma_min_ratio
gamma_2_min = gamma_2_max * gamma_min_ratio

gamma_1_values = 10 ^ seq(from=log(gamma_1_min, base=10), to=log(gamma_1_max, base=10), length.out=num_gamma)
gamma_2_values = 10 ^ seq(from=log(gamma_2_min, base=10), to=log(gamma_2_max, base=10), length.out=num_gamma)

ic_values = matrix(NA, nrow=length(gamma_1_values), ncol=length(gamma_2_values))
for (i in 1:length(gamma_1_values)){
  for (j in 1:length(gamma_2_values)){
    gamma_1 = gamma_1_values[i]
    gamma_2 = gamma_2_values[j]
    ic_ij = compute_information_criterion(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa, gamma_1=gamma_1, gamma_2=gamma_2)
    ic_values[i,j] = ic_ij
    clog(paste0("gamma_1=", gamma_1, ", gamma_2=", gamma_2, ", ic=", round(ic_ij,2)))
  }
}

pos_best = which(ic_values == min(ic_values), arr.ind=TRUE)
gamma_1_best = gamma_1_values[pos_best[1,"row"]]
gamma_2_best = gamma_2_values[pos_best[1,"col"]]
clog(paste0("gamma_1_best=", round(gamma_1_best,4), ", gamma_2_best=", round(gamma_2_best,4)))

# Step 3: estimation of slopes with given gamma_1 and gamma_2
fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa, gamma_1=gamma_1_best, gamma_2=gamma_2_best)
A_hat = fit$A_hat
group_hat = fit$group_hat
regime_hat = fit$regime_hat
clog(paste0("group_hat=", vec2str(group_hat), ", regime_hat=", vec2str(regime_hat)))


