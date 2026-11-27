# This module consists of all functions for running simulations
source("dgp.r")
source("est.r")
source("util.r")

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

set.seed(1233)
make_csv_Btilde(N=N, T_=T_, p=p, filename=CSV_BETA)

run_replication = function(N, T_, p, r, phi, K0, m0, pi_, L, kappa){
  # Step 1: generate observations X and y
  # X: N x T_ x p
  # y: N x T_
  data = generate_data(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, m0=m0, pi_=pi_, L=L)
  X = data$X
  y = data$y
  B0 = data$B
  assert_dimensions(X, dim_=c(N, T_, p))
  assert_dimensions(y, dim_=c(N, T_))
  
  # Step 2: make the optimizer for this data set
  optimizer = make_optimizer(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa)
  
  # Step 3: tune gamma_1 and gamma_2
  gamma_max = tune_gamma_max(N=N, T_=T_, p=p, X=X, y=y, opt=optimizer, start_value=gamma_max_start_value)
  gamma_grid = get_gamma_grid(gamma_max=gamma_max, gamma_min_ratio=gamma_min_ratio, num_gamma=num_gamma)
  gamma_hat = tune_gamma(N=N, T_=T_, p=p, X=X, y=y, opt=optimizer, gamma_grid=gamma_grid)
  
  # Step 4: estimate slopes with given gamma_1 and gamma_2
  # A_hat: K_hat x (m_hat+1) x p
  fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, opt=optimizer, gamma_1=gamma_hat[1], gamma_2=gamma_hat[2])
  
  # Step 5: evaluate the performance
  metrics = eval_perf(N=N, T_=T_, p=p, B0=B0, fit=fit)
  return(metrics)
}

# Run simulation for n_rep times
metrics = list(
  correct_num_groups=rep(NA, n_rep),
  correct_membership=rep(NA, n_rep),
  correct_num_breaks=rep(NA, n_rep),
  correct_break_dates=rep(NA, n_rep),
  num_groups=rep(NA, n_rep),
  num_breaks=rep(NA, n_rep),
  sum_squared_errors=rep(NA, n_rep),
  size=rep(NA, n_rep)
)

for (rep in 1:n_rep){
  metrics_rep = run_replication(N=N, T_=T_, p=p, r=r, phi=phi, K0=K0, m0=m0, pi_=pi_, L=L, kappa=kappa)
  
  metrics$correct_num_groups[rep] = metrics_rep$correct_num_groups
  metrics$correct_membership[rep] = metrics_rep$correct_membership
  metrics$correct_num_breaks[rep] = metrics_rep$correct_num_breaks
  metrics$correct_break_dates[rep] = metrics_rep$correct_break_dates
  metrics$num_groups[rep] = metrics_rep$num_groups
  metrics$num_breaks[rep] = metrics_rep$num_breaks
  metrics$sum_squared_errors[rep] = metrics_rep$sum_squared_errors
  metrics$size[rep] = metrics_rep$size
}

correct_num_groups_and_num_breaks = metrics$correct_num_groups & metrics$correct_num_breaks
metrics_agg = list(
  pct_wrong_num_groups=1-mean(metrics$correct_num_groups),
  pct_wrong_membership=1-mean(metrics$correct_membership[metrics$correct_num_groups]),
  pct_wrong_num_breaks=1-mean(metrics$correct_num_breaks),
  pct_wrong_break_dates=1-mean(metrics$correct_break_dates[metrics$correct_num_breaks]),
  pct_wrong_classification=1-mean(metrics$correct_membership[correct_num_groups_and_num_breaks] &
                                    metrics$correct_break_dates[correct_num_groups_and_num_breaks]),
  avg_num_groups=mean(metrics$num_groups),
  avg_num_breaks=mean(metrics$num_breaks)
)
clog(paste0("performance: ", list2str(metrics_agg)))
  




