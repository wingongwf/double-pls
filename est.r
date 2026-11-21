# This module consists of all functions for estimation.

compute_M_F_hat = function(N, T_, p, X){
  # Compute the matrix that projects onto the space orthogonal to the
  # the col space of factors. We only use the cross-sectional means of X (as
  # opposed to both X and y) as a proxy of space spanned by F.
  # F_hat = Z_bar
  # M_F_hat: T_ x T_
  Z = array(NA, c(N, T_, p))
  for (i in 1:N){
    assert_dimensions(Z[i,,], X[i,,])
    Z[i,,] = X[i,,]
  }
  
  F_hat = apply(Z, c(2,3), mean)  # compute cross-sectional averages
  assert_dimensions(F_hat, dim_=c(T_, p))
  
  P_F_hat = F_hat %*% MASS::ginv(t(F_hat)%*%F_hat) %*% t(F_hat)
  assert_dimensions(base::diag(T_), P_F_hat)
  M_F_hat = base::diag(T_) - P_F_hat
  
  assert_dimensions(M_F_hat, dim_=c(T_, T_))
  return(M_F_hat)
}

compute_residual_sum_squares = function(N, T_, p, X, y, B, M_F_hat){
  residual_sum_squares = 0
  for (i in 1:N){
    B_i = t(get_B_i(B=B, i=i, T_=T_, p=p))
    assert_dimensions(X[i,,], B_i)
    
    term_fit = (X[i,,] * B_i) %*% rep(1, p)
    residual_sum_squares_i = CVXR::sum_squares(M_F_hat %*% (y[i,] - term_fit))
    assert_dimensions(residual_sum_squares_i, dim_=c(1, 1))
    
    residual_sum_squares = residual_sum_squares + residual_sum_squares_i
  }
  
  assert_dimensions(residual_sum_squares, dim_=c(1, 1))
  return(residual_sum_squares)
}

estimate_slopes_by_cross_section = function(N, p, X, y, M_F_hat){
  B_cs = matrix(NA, p, N)
  for (i in 1:N){
    B_cs_i = solve(t(X[i,,]) %*% M_F_hat %*% X[i,,]) %*% t(X[i,,]) %*% M_F_hat %*% y[i,]
    assert_dimensions(B_cs_i, dim_=c(p, 1))
    B_cs_i = B_cs_i[,1]
    assert_dimensions(B_cs[,i], B_cs_i)
    B_cs[,i] = B_cs_i
  }
  
  return(B_cs)
}

compute_omega_cross_section = function(N, p, X, y, M_F_hat, kappa){
  B_cs = estimate_slopes_by_cross_section(N=N, p=p, X=X, y=y, M_F_hat=M_F_hat)
  assert_dimensions(B_cs, dim_=c(p, N))
  
  omega = matrix(NA, N, N)
  for (i in 1:N){
    for (j in 1:N){
      if (i < j){
        omega_ij = sqrt(sum((B_cs[,i] - B_cs[,j])^2)) ^ (-kappa)
        stopifnot(length(omega_ij)==1)
        omega[i,j] = omega_ij
      }
    }
  }
  
  stopifnot(sum(!is.na(omega)) == (N-1)*N/2)
  return(omega)
}

compute_penalty_cross_section = function(N, T_, p, X, y, B, M_F_hat, kappa){
  penalty_cross_section = 0
  for (i in 1:N){
    for (j in 1:N){
      if (i < j){
        B_i = get_B_i(B=B, i=i, T_=T_, p=p)
        B_j = get_B_i(B=B, i=j, T_=T_, p=p)
        
        difference = B_i - B_j
        assert_dimensions(difference, dim_=c(p, T_))
        
        difference = CVXR::vec(difference)
        assert_dimensions(difference, dim_=c(p*T_, 1))
        
        omega_cs = compute_omega_cross_section(N=N, p=p, X=X, y=y, M_F_hat=M_F_hat, kappa=kappa)
        penalty_ij = omega_cs[i,j] * CVXR::p_norm(difference, 2)
        assert_dimensions(penalty_ij, dim_=c(1, 1))
        penalty_cross_section = penalty_cross_section + penalty_ij
      }
    }
  }
  
  assert_dimensions(penalty_cross_section, dim_=c(1, 1))
  return(penalty_cross_section)
}

estimate_slopes_by_period = function(N, T_, p, X, y){
  B_period = matrix(NA, p, T_)
  X_cross_section_mean = apply(X, c(2,3), mean)
  y_cross_section_mean = apply(y, 2, mean)
  assert_dimensions(X_cross_section_mean, dim_=c(T_, p))
  assert_dimensions(y_cross_section_mean, dim_=T_)
  
  for (t in 1:T_){
    sum_1 = 0
    sum_2 = 0
    
    for (i in 1:N){
      assert_dimensions(X[i,t,], dim_=p)
      X_tilde_it = X[i,t,] - X_cross_section_mean[t,]
      
      assert_dimensions(y[i,t], dim_=1)
      y_tilde_it = y[i,t] - y_cross_section_mean[t]
      
      XX = X_tilde_it %*% t(X_tilde_it)
      assert_dimensions(XX, dim_=c(p, p))
      sum_1 = sum_1 + XX
      
      Xy = X_tilde_it * y_tilde_it
      assert_dimensions(Xy, dim_=p)
      sum_2 = sum_2 + Xy
    }
    
    B_period_t = solve(sum_1) %*% sum_2
    assert_dimensions(B_period_t, dim_=c(p, 1))
    B_period_t = B_period_t[,1]
    
    assert_dimensions(B_period[,t], B_period_t)
    B_period[,t] = B_period_t
  }
  
  return(B_period)
}

compute_omega_period = function(N, T_, p, X, y, kappa){
  B_period = estimate_slopes_by_period(N=N, T_=T_, p=p, X=X, y=y)
  assert_dimensions(B_period, dim_=c(p, T_))
  
  omega = rep(NA, T_)
  for (t in 2:T_){
    omega_t = sqrt(sum((B_period[,t] - B_period[,t-1])^2)) ^ (-kappa)
    assert_dimensions(omega_t, dim_=1)
    omega[t] = omega_t
  }
  
  stopifnot(sum(!is.na(omega)) == T_-1)
  return(omega)
}

compute_penalty_time_dimension = function(N, T_, p, X, y, B, kappa){
  penalty_time_dimension = 0
  for (t in 2:T_){
    B_t = get_B_t(B=B, t=t, T_=T_, N=N, p=p)
    B_tminus1 = get_B_t(B=B, t=t-1, T_=T_, N=N, p=p)
    
    difference = B_t - B_tminus1
    assert_dimensions(difference, dim_=c(p, N))
    
    difference = CVXR::vec(difference)
    assert_dimensions(difference, dim_=c(p*N, 1))
    
    omega_period = compute_omega_period(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa)
    penalty_t = omega_period[t] * CVXR::p_norm(difference, 2)
    assert_dimensions(penalty_t, dim_=c(1, 1))
    penalty_time_dimension = penalty_time_dimension + penalty_t
  }
  
  assert_dimensions(penalty_time_dimension, dim_=c(1,1))
  return(penalty_time_dimension)
}


