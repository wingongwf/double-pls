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
# Decision variables
# B: p x N*T_
B = CVXR::Variable(c(p, N*T_))

# Set up the objective function
M_F_hat = compute_M_F_hat(N=N, T_=T_, p=p, X=X)
residual_sum_squares = compute_residual_sum_squares(N=N, T_=T_, p=p, X=X, y=y, B=B, M_F_hat=M_F_hat)
penalty_cross_section = compute_penalty_cross_section(N=N, T_=T_, p=p, X=X, y=y, B=B, M_F_hat=M_F_hat, kappa=kappa)
penalty_time_dimension = compute_penalty_time_dimension(N=N, T_=T_, p=p, X=X, y=y, B=B, kappa=kappa)

objective = CVXR::Minimize(1/(N*T_) * residual_sum_squares +
                             gamma_1/(N^2) * penalty_cross_section +
                             gamma_2/T_ * penalty_time_dimension)

# Run optimization to estimate B: B_hat
prob = CVXR::Problem(objective)
solution = CVXR::psolve(prob)
stopifnot(CVXR::status(prob) == "optimal")  # assert convergence
num_iters = CVXR::solver_stats(prob)@num_iters
solve_time = CVXR::solver_stats(prob)@solve_time
clog(paste0("solve_time=", round(solve_time,2), ", ", "num_iters=", num_iters))

B_hat = CVXR::value(B)
assert_dimensions(B_hat, dim_=c(p, N*T_))


