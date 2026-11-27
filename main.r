# In this program, we aim to study the estimation and identification of latent
# group structures and multiple breaks in panel data models through simulations.
source("util.r")
source("dgp.r")
source("est.r")
source("sim.r")

set.seed(1233)

# # parameters for simulation
# n_rep = 10
# 
# # parameters for data generating process
# N = 10
# T_ = 20
# p = 2
# r = 2
# phi = 0.8
# K0 = 2
# m0 = 1
# pi_ = 0.2
# L = 5
# kappa = 2
# 
# # parameters for tuning gamma
# gamma_max_start_value = 1
# gamma_min_ratio = 1e-4
# num_gamma = 5
# 
# # file names
# LOG_DATETIME = format(Sys.time(), "%Y%m%d_%H%M%S")
# LOG_FILE = paste0("logs/simulation_", LOG_DATETIME, ".log")
# CSV_BETA = paste0("results/Btilde_", LOG_DATETIME, ".csv")
# 
# run_simulation(n_rep=n_rep, N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, m0=m0, pi_=pi_, L=L, kappa=kappa)

# file names
LOG_DATETIME = format(Sys.time(), "%Y%m%d_%H%M%S")
LOG_FILE = paste0("logs/simulation_", LOG_DATETIME, ".log")
CSV_BETA = paste0("results/Btilde_", LOG_DATETIME, ".csv")

if (interactive()){
  # when running in RStudio console
  clog("running in console")
  params = list(
    # parameters for simulation
    n_rep=10,
    
    # parameters for data generating process
    N=10,   
    T_=20,
    p=2,
    r=2,
    phi=0.8,
    K0=2,
    m0=1,
    pi_=0.2,
    L=5,
    kappa=2,
    
    # parameters for tuning gamma
    gamma_max_start_value=1,
    gamma_min_ratio=1e-4,
    num_gamma=5
  )
} else{
  # when running in terminal
  # Rscript main.r --n_rep=10 --N=10 --T_=20 --p=2 --r=2 --phi=0.8 --K0=2 --m0=1 --pi_=0.2 --L=5 --kappa=2 --gamma_max_start_value=1 --gamma_min_ratio=1e-4 --num_gamma=5
  clog("running in terminal")
  args = commandArgs(trailingOnly=TRUE)
  params = read_args(args)
}

do.call(run_simulation, params)


