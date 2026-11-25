# This module consists of all the functions for the data generating process.
# For arrays, we use the following order for dimensions:
# N > T_ > r > p
source("util.r")

# Generate factors
generate_factors = function(r, T_, phi){
  f_0 = rep(0, r)
  F_ = matrix(NA, nrow=T_+1, ncol=r)
  F_[1,] = f_0    # initialize f_0 for the loop
  for (t in 1:T_){
    eta_t = MASS::mvrnorm(n=1, mu=rep(0,r), Sigma=diag(r))
    F_[t+1,] = (1-phi)*rep(1,r) + phi*F_[t,] + eta_t
  }
  
  F_ = F_[-1, ]    # remove f_0
  return(F_)
}

# Generate loadings for unit i
generate_loadings_i = function(r, p){
  Lambda_i = matrix(rnorm(n=r*p, mean=2, sd=1), nrow=r, ncol=p)
  lambda_i = 0.2 * Lambda_i[,1]
  
  return(list(Lambda_i=Lambda_i, lambda_i=lambda_i))
}

# Generate loadings for all units
generate_loadings = function(N, r, p){
  Lambda = array(NA, c(N, r, p))
  lambda = matrix(NA, N, r)
  for (i in 1:N){
    loadings = generate_loadings_i(r=r, p=p)
    Lambda_i = loadings$Lambda_i
    lambda_i = loadings$lambda_i
    
    assert_dimensions(Lambda[i,,], Lambda_i)
    assert_dimensions(lambda[i,], lambda_i)
    
    Lambda[i,,] = Lambda_i
    lambda[i,] = lambda_i
  }
  
  return(list(Lambda=Lambda, lambda=lambda))
}

# Generate slopes for unit i at time t
generate_slopes_it = function(i, t, N, T_, p, K0, m0){
  # some checks
  if (i>N){
    stop("i must be samller than or equal to N")
  }
  if (K0!=2){
    stop("K0 must be 2")
  }
  if (m0!=1){
    stop("m0 must be 1")
  }
  
  if (i<=0.4*N){
    if (t<=0.4*T_){
      B_it = rep(0, p)
    } else{
      B_it = rep(1, p)
    }
  } else{
    if (t<=0.4*T_){
      B_it = rep(2, p)
    } else{
      B_it = rep(3, p)
    }
  }
  
  return(B_it)
}

# Generate slopes for all units at all times
generate_slopes = function(N, T_, p, K0, m0){
  B = array(NA, c(N, T_, p))
  for (i in 1:N){
    for (t in 1:T_){
      B_it = generate_slopes_it(i=i, t=t, N=N, T_=T_, p=p, K0=K0, m0=m0)
      assert_dimensions(B[i,t,], B_it)
      B[i,t,] = B_it
    }
  }
  
  return(B)
}

# Generate error standard deviations
generate_errors_sd = function(N){
  sigma_values = runif(n=N, min=0.5, max=1)
  return(sigma_values)
}

# Generate errors
# Note that the idiosyncratic errors are weakly serially and cross-sectionally
# correlated, and so the errors need to be generated together. There is no
# correlation between variables though.
generate_errors_k = function(pi_, sigma_values, N, T_, L){
  theta = matrix(NA, nrow=N, ncol=T_)
  for (i in 1:N){
    theta[i,] = rnorm(n=T_, mean=0, sd=sigma_values[i])
  }
  
  e = matrix(NA, nrow=N, ncol=T_)
  for (i in 1:N){
    e_i = rep(NA, T_+1)
    e_0 = 0
    e_i[1] = e_0   # initialize e_0 for the loop
    
    for (t in 1:T_){
      e_it = pi_ * e_i[t] + theta[i,t]
      for (j in 1:L){
        index_back = get_circular_index(pos=i, step=-j, size=N)
        index_front = get_circular_index(pos=i, step=j, size=N)
        e_it = e_it + pi_ * (theta[index_back,t] + theta[index_front,t])
      }
      e_i[t+1] = e_it
    }
    
    e_i = e_i[-1]  # remove e_0
    e[i,] = e_i
  }
  
  return(e)
}

# Generate errors in X
generate_errors_V = function(N, T_, p, pi_, sigma_values, L){
  V = array(NA, c(N,T_,p))
  for (k in 1:p){
    V_k = generate_errors_k(pi_=pi_, sigma_values=sigma_values, N=N, T_=T_, L=L)
    assert_dimensions(V[,,k], V_k)
    V[,,k] = V_k
  }
  
  return(V)
}

# Generate errors in y
generate_errors_eps = function(N, T_, pi_, sigma_values, L){
  eps = generate_errors_k(pi_=pi_, sigma_values=sigma_values, N=N, T_=T_, L=L)
  return(eps)
}

# Compute X
compute_X = function(F_, Lambda, V){
  num_row = dim(Lambda)[1]
  num_col = dim(F_)[1]
  num_slice = dim(V)[3]
  
  X = array(NA, dim=c(num_row, num_col, num_slice))
  for (i in 1:num_row){
    X_i = F_ %*% Lambda[i,,] + V[i,,]
    assert_dimensions(X[i,,], X_i)
    X[i,,] = X_i
  }
  
  return(X)
}

# Compute y
compute_y = function(X, B, F_, lambda, eps){
  N = dim(eps)[1]
  T_ = dim(eps)[2]
  
  y = matrix(NA, N, T_)
  for (i in 1:N){
    y_i = (X[i,,] * B[i,,]) %*% rep(1, p) + F_ %*% lambda[i,] + eps[i,]
    assert_dimensions(y_i, dim_=c(T_, 1))
    y_i = y_i[,1]
    assert_dimensions(y[i,], y_i)
    y[i,] = y_i
  }
  
  return(y)
}

# Generate observations X and y
generate_data = function(N, T_, p, r, phi, K0, m0, pi_, L){
  # error standard deviations
  # sigma_values: N x 1
  sigma_values = generate_errors_sd(N)
  
  # factors
  # F_: T_ x r
  F_ = generate_factors(r=r, T_=T_, phi=phi)
  assert_dimensions(F_, dim_=c(T_, r))
  
  # loadings
  # Lambda: N x r x p
  # lambda: N x r
  loadings = generate_loadings(N=N, r=r, p=p)
  Lambda = loadings$Lambda
  lambda = loadings$lambda
  assert_dimensions(Lambda, dim_=c(N, r, p))
  assert_dimensions(lambda, dim_=c(N, r))
  
  # errors
  # V: N x T_ x p
  # eps: N x T_
  V = generate_errors_V(N=N, T_=T_, p=p, pi_=pi_, sigma_values=sigma_values, L=L)
  assert_dimensions(V, dim_=c(N, T_, p))
  
  eps = generate_errors_eps(N=N, T_=T_, pi_=pi_, sigma_values=sigma_values, L=L)
  assert_dimensions(eps, dim_=c(N, T_))
  
  # slopes
  # B: N x p
  B = generate_slopes(N=N, T_=T_, p=p, K0=K0, m0=m0)
  assert_dimensions(B, dim_=c(N, T_, p))
  
  # regressors
  # X: N x T_ x p
  X = compute_X(F_=F_, Lambda=Lambda, V=V)
  assert_dimensions(X, dim_=c(N, T_, p))
  
  # target
  # y: N x T_
  y = compute_y(X=X, B=B, F_=F_, lambda=lambda, eps=eps)
  assert_dimensions(y, dim_=c(N, T_))
  
  return(list(X=X, y=y))
}


