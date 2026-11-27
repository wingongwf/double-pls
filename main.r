# In this program, we aim to study the estimation and identification of latent
# group structures and multiple breaks in panel data models through simulations.
source("util.r")
source("dgp.r")
source("est.r")
source("sim.r")

set.seed(1233)

# parameters for simulation
n_rep = 10  # TODO: small parameters test case

# parameters for data generating process
N = 10     # TODO: small parameters test case
T_ = 20    # TODO: small parameters test case
p = 2      # TODO: small parameters test case
r = 2
phi = 0.8
K0 = 2
m0 = 1
pi_ = 0.2
L = 5
kappa = 2

# parameters for tuning gamma
gamma_max_start_value = 1
gamma_min_ratio = 1e-4
num_gamma = 5

# file names
LOG_DATETIME = format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE = paste0("logs/simulation_", LOG_DATETIME, ".log")
CSV_BETA = paste0("results/Btilde_", LOG_DATETIME, ".csv")

run_simulation(n_rep=n_rep, N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, m0=m0, pi_=pi_, L=L, kappa=kappa)


