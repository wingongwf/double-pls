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

compute_penalty_cross_section = function(N, T_, p, X, y, B, M_F_hat, kappa){
  omega_cs = compute_omega_cross_section(N=N, T_=T_, p=p, X=X, y=y, M_F_hat=M_F_hat, kappa=kappa)
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
        
        penalty_ij = omega_cs[i,j] * CVXR::p_norm(difference, 2)
        assert_dimensions(penalty_ij, dim_=c(1, 1))
        penalty_cross_section = penalty_cross_section + penalty_ij
      }
    }
  }
  
  assert_dimensions(penalty_cross_section, dim_=c(1, 1))
  return(penalty_cross_section)
}

compute_penalty_time_dimension = function(N, T_, p, X, y, B, kappa){
  omega_period = compute_omega_period(N=N, T_=T_, p=p, X=X, y=y, kappa=kappa)
  penalty_time_dimension = 0
  for (t in 2:T_){
    B_t = get_B_t(B=B, t=t, T_=T_, N=N, p=p)
    B_tminus1 = get_B_t(B=B, t=t-1, T_=T_, N=N, p=p)
    
    difference = B_t - B_tminus1
    assert_dimensions(difference, dim_=c(p, N))
    
    difference = CVXR::vec(difference)
    assert_dimensions(difference, dim_=c(p*N, 1))
    
    penalty_t = omega_period[t] * CVXR::p_norm(difference, 2)
    assert_dimensions(penalty_t, dim_=c(1, 1))
    penalty_time_dimension = penalty_time_dimension + penalty_t
  }
  
  assert_dimensions(penalty_time_dimension, dim_=c(1,1))
  return(penalty_time_dimension)
}

compute_omega_cross_section = function(N, T_, p, X, y, M_F_hat, kappa){
  B_cs = estimate_slopes_by_cross_section(N=N, T_=T_, p=p, X=X, y=y)
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

compute_information_criterion = function(N, T_, p, X, y, opt, gamma_1, gamma_2, verbose=FALSE){
  fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, opt=opt, gamma_1=gamma_1, gamma_2=gamma_2, verbose=verbose)
  A_hat = fit$A_hat
  M_F_hat = fit$M_F_hat
  group_hat = fit$group_hat
  regime_hat = fit$regime_hat
  K_hat = fit$K_hat
  m_hat = fit$m_hat
  
  if (verbose) clog(paste0("gamma_1=", round(gamma_1,2), ", gamma_2=", round(gamma_2,2), ", K_hat=", K_hat, ", m_hat=", m_hat))
  
  residual_sum_squares = 0
  for (l in 1:K_hat){
    for (u in 1:(m_hat+1)){
      periods = which(regime_hat==u)
      num_periods = length(periods)
      
      for (i in 1:N){
        if (group_hat[i] == l){
          X_iu = X[i,periods,]
          y_iu = y[i,periods]
          M_F_hat_iu = M_F_hat[periods,periods]
          A_hat_lu = A_hat[l,u,]
          
          # reshape vectors to conformable matrices when num_periods=1
          if (num_periods == 1){
            X_iu = matrix(data=X_iu, nrow=num_periods, ncol=p)
            M_F_hat_iu = matrix(data=M_F_hat_iu, nrow=num_periods, ncol=num_periods)
          }
          
          # check dimensions
          assert_dimensions(X_iu, dim_=c(num_periods, p))
          assert_dimensions(y_iu, dim_=num_periods)
          assert_dimensions(M_F_hat_iu, dim_=c(num_periods, num_periods))
          assert_dimensions(A_hat_lu, dim_=p)
          
          residual_sum_squares_iu = sum((M_F_hat_iu %*% (y_iu - X_iu %*% A_hat_lu))^2)
          assert_dimensions(residual_sum_squares_iu, dim_=1)
          residual_sum_squares = residual_sum_squares + residual_sum_squares_iu
        }
      }
    }
  }
  
  phi = 0.5 * log(N*T_) / sqrt(N*T_)
  ic = 1/(N*T_) * residual_sum_squares + phi * p * (K_hat + m_hat + 1)
  
  return(ic)
}

tune_gamma = function(N, T_, p, X, y, opt, gamma_grid){
  num_gamma_1 = length(gamma_grid$gamma_1)
  num_gamma_2 = length(gamma_grid$gamma_2)
  ic_values = matrix(NA, nrow=num_gamma_1, ncol=num_gamma_2)
  
  for (i in 1:num_gamma_1){
    for (j in 1:num_gamma_2){
      gamma_1 = gamma_grid$gamma_1[i]
      gamma_2 = gamma_grid$gamma_2[j]
      ic_ij = compute_information_criterion(N=N, T_=T_, p=p, X=X, y=y, opt=opt, gamma_1=gamma_1, gamma_2=gamma_2, verbose=FALSE)
      ic_values[i,j] = ic_ij
      # clog(paste0("gamma_1=", gamma_1, ", gamma_2=", gamma_2, ", ic=", round(ic_ij,2)))
    }
  }
  
  stopifnot(sum(is.na(ic_values)) == 0)
  pos_min_ic = which(ic_values == min(ic_values), arr.ind=TRUE)
  
  # Select the first row (even when there are multiple gamma pairs that give same IC)
  gamma_1_hat = gamma_grid$gamma_1[pos_min_ic[1,"row"]]
  gamma_2_hat = gamma_grid$gamma_2[pos_min_ic[1,"col"]]
  gamma_hat = c(gamma_1_hat, gamma_2_hat)
  clog(paste0("pos_min_ic=", vec2str(apply(pos_min_ic, 1, paste, collapse=",")),
              ", gamma_hat=", vec2str(round(gamma_hat,4)),
              ", gamma_hat_ratio=", round(gamma_1_hat/gamma_2_hat, 4)))
  
  return(gamma_hat)
}

tune_gamma_max = function(N, T_, p, X, y, opt, start_value, multiplier=10, lower_bound=1e-10, upper_bound=1e10){
  gamma = start_value
  fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, opt=opt, gamma_1=gamma, gamma_2=gamma, verbose=FALSE)
  
  if (fit$K_hat!=1 || fit$m_hat!=0){
    # if gamma too small
    while (fit$K_hat!=1 || fit$m_hat!=0){
      prev_gamma = gamma
      gamma = gamma * multiplier
      if (gamma > upper_bound){
        stop("upper_bound=", upper_bound, " reached, gamma_max not found")
      }
      
      clog(paste0("prev_gamma=", prev_gamma," too small, try gamma=", gamma))
      fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, opt=opt, gamma_1=gamma, gamma_2=gamma, verbose=FALSE)
    }
    
    return(c(gamma, gamma))
  } else{
    # if gamma too large
    while (fit$K_hat==1 && fit$m_hat==0){
      prev_gamma = gamma
      gamma = gamma / multiplier
      if (gamma < lower_bound){
        stop("lower_bound=", lower_bound, " reached, gamma_max not found")
      }
      
      clog(paste0("prev_gamma=", prev_gamma," too large, try gamma=", gamma))
      fit = estimate_slopes(N=N, T_=T_, p=p, X=X, y=y, opt=opt, gamma_1=gamma, gamma_2=gamma, verbose=FALSE)
    }
    
    gamma_max = gamma * multiplier
    return(c(gamma_max, gamma_max))
  }
}

make_optimizer = function(N, T_, p, X, y, M_F_hat, kappa){
  # Decision variables
  # B: p x N*T_
  B = CVXR::Variable(c(p, N*T_))
  
  # Parameters: gamma_1, gamma_2
  gamma_1 = CVXR::Parameter(nonneg=TRUE)
  gamma_2 = CVXR::Parameter(nonneg=TRUE)
  
  # Set up the objective function
  M_F_hat = compute_M_F_hat(N=N, T_=T_, p=p, X=X)
  residual_sum_squares = compute_residual_sum_squares(N=N, T_=T_, p=p, X=X, y=y, B=B, M_F_hat=M_F_hat)
  penalty_cross_section = compute_penalty_cross_section(N=N, T_=T_, p=p, X=X, y=y, B=B, M_F_hat=M_F_hat, kappa=kappa)
  penalty_time_dimension = compute_penalty_time_dimension(N=N, T_=T_, p=p, X=X, y=y, B=B, kappa=kappa)
  
  objective = CVXR::Minimize(1/(N*T_) * residual_sum_squares +
                               gamma_1/(N^2) * penalty_cross_section +
                               gamma_2/T_ * penalty_time_dimension)
  
  # Set up the problem
  prob = CVXR::Problem(objective)
  return(list(prob=prob, B=B, gamma_1=gamma_1, gamma_2=gamma_2))
}

estimate_groups = function(N, T_, p, B, verbose=FALSE, tol=1e-8){
  # Put cross-sectional units into groups
  group = rep(NA, N)
  group[1] = 1   # put i=1 into first group
  
  for (i in 2:N){
    B_hat_i = get_B_i(B=B, i=i, T_=T_, p=p)
    is_assigned = FALSE
    collected_groups = unique(group)
    collected_groups = collected_groups[!is.na(collected_groups)]   # remove NA
    
    for (l in collected_groups){
      is_group_l = TRUE
      members = which(group==l)  # all members of group l
      
      # if B_i is not equal to ANY member, then not group l
      for (mem in members){
        B_hat_mem = get_B_i(B=B, i=mem, T_=T_, p=p)
        if (max(abs(as.vector(B_hat_i - B_hat_mem))) > tol){
          if (verbose) clog(paste0("i=", i, " different than i=", mem, " in group=", l))
          is_group_l = FALSE
          break
        } else{
          if (verbose) clog(paste0("i=", i, " same as i=", mem, " in group=", l))
        }
      }
      
      # if B_i is equal to ALL members, then assign to group l
      if (is_group_l){
        if (verbose) clog(paste0("i=", i, " same as i=", vec2str(members), ", assign to group=", l))
        group[i] = l
        is_assigned = TRUE
        break
      }
    }
    
    # if B_i does not belong to any collected groups, then assign to a new group
    if(!is_assigned){
      new_group = max(collected_groups) + 1
      if (verbose) clog(paste0("i=", i, " assign to new group=", new_group))
      group[i] = new_group
    }
  }
  
  return(group)
}

estimate_regimes = function(N, T_, p, B, verbose=FALSE, tol=1e-8){
  # Put periods into regimes
  regime = rep(NA, T_)
  regime[1] = 1   # put t=1 into first regime
  
  for (t in 2:T_){
    B_t = get_B_t(B=B, t=t, T_=T_, N=N, p=p)
    B_tminus1 = get_B_t(B=B, t=t-1, T_=T_, N=N, p=p)
    
    if (max(abs(as.vector(B_t - B_tminus1))) > tol){
      if (verbose) clog(paste0("t=", t, " different regime than t=", t-1, ", assign to regime=", regime[t-1] + 1))
      regime[t] = regime[t-1] + 1 
    } else{
      if (verbose) clog(paste0("t=", t, " same regime as t=", t-1, ", assign to regime=", regime[t-1]))
      regime[t] = regime[t-1]
    }
  }
  
  return(regime)
}

estimate_slopes = function(N, T_, p, X, y, opt, gamma_1, gamma_2, verbose=FALSE){
  # Step 1: estimate the slopes B with given gamma_1 and gamma_2
  M_F_hat = compute_M_F_hat(N=N, T_=T_, p=p, X=X)
  B_hat = estimate_slopes_pre_selection(opt=opt, gamma_1=gamma_1, gamma_2=gamma_2, verbose=verbose)
  assert_dimensions(B_hat, dim_=c(p, N*T_))
  
  # Step 2: estimate groups and regimes in slopes
  group_hat = estimate_groups(N=N, T_=T_, p=p, B=B_hat)
  regime_hat = estimate_regimes(N=N, T_=T_, p=p, B=B_hat)
  K_hat = max(group_hat)
  m_hat = max(regime_hat) - 1
  if (verbose) clog(paste0("group=", vec2str(group_hat), ", regime=", vec2str(regime_hat), ", K_hat=", K_hat, ", m_hat=", m_hat))
  
  # Step 3: post-selection estimation of slopes
  # A_hat: K_hat x (m_hat+1) x p
  A_hat = estimate_slopes_post_selection(
    N=N, p=p, X=X, y=y, M_F_hat=M_F_hat, group_hat=group_hat, regime_hat=regime_hat, K_hat=K_hat, m_hat=m_hat, verbose=verbose)
  assert_dimensions(A_hat, dim_=c(K_hat, m_hat+1, p))
  
  return(list(A_hat=A_hat, M_F_hat=M_F_hat, group_hat=group_hat, regime_hat=regime_hat, K_hat=K_hat, m_hat=m_hat))
}

estimate_slopes_pre_selection = function(opt, gamma_1, gamma_2, verbose=TRUE){
  # Update parameters
  CVXR::value(opt$gamma_1) = gamma_1
  CVXR::value(opt$gamma_2) = gamma_2
  
  # Run optimization to estimate B: B_hat
  solution = CVXR::psolve(opt$prob)
  
  # Log the optimization stats
  stats = CVXR::solver_stats(opt$prob)
  status = CVXR::status(opt$prob)
  if (verbose) clog(paste0("gamma_1=", gamma_1, ", gamma_2=", gamma_2, ", status=", status,
              ", solve_time=", round(stats@solve_time,2), ", ", "num_iters=", stats@num_iters))
  stopifnot(status == "optimal")  # assert convergence

  B_hat = CVXR::value(opt$B)
  
  # Revert parameters back to NULL for safety
  CVXR::value(opt$gamma_1) = NULL
  CVXR::value(opt$gamma_2) = NULL
  
  return(B_hat)
}

estimate_slopes_by_cross_section = function(N, T_, p, X, y){
  B_cs = matrix(NA, p, N)
  X_time_mean = apply(X, c(1,3), mean)
  y_time_mean = apply(y, 1, mean)
  assert_dimensions(X_time_mean, dim_=c(N,p))
  assert_dimensions(y_time_mean, dim_=N)

  for (i in 1:N){
    sum_XX = 0
    sum_Xy = 0

    for (t in 1:T_){
      assert_dimensions(X[i,t,], dim_=p)
      X_ddot_it = X[i,t,] - X_time_mean[i,]

      assert_dimensions(y[i,t], dim_=1)
      y_ddot_it = y[i,t] - y_time_mean[i]

      XX = X_ddot_it %*% t(X_ddot_it)
      assert_dimensions(XX, dim_=c(p, p))
      sum_XX = sum_XX + XX

      Xy = X_ddot_it * y_ddot_it
      assert_dimensions(Xy, dim_=p)
      sum_Xy = sum_Xy + Xy
    }

    B_cs_i = solve(sum_XX) %*% sum_Xy
    assert_dimensions(B_cs_i, dim_=c(p, 1))
    B_cs_i = B_cs_i[,1]

    assert_dimensions(B_cs[,i], B_cs_i)
    B_cs[,i] = B_cs_i
  }

  return(B_cs)
}

estimate_slopes_by_cross_section_CCE = function(N, p, X, y, M_F_hat){
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

estimate_slopes_by_period = function(N, T_, p, X, y){
  B_period = matrix(NA, p, T_)
  X_cross_section_mean = apply(X, c(2,3), mean)
  y_cross_section_mean = apply(y, 2, mean)
  assert_dimensions(X_cross_section_mean, dim_=c(T_, p))
  assert_dimensions(y_cross_section_mean, dim_=T_)
  
  for (t in 1:T_){
    sum_XX = 0
    sum_Xy = 0
    
    for (i in 1:N){
      assert_dimensions(X[i,t,], dim_=p)
      X_tilde_it = X[i,t,] - X_cross_section_mean[t,]
      
      assert_dimensions(y[i,t], dim_=1)
      y_tilde_it = y[i,t] - y_cross_section_mean[t]
      
      XX = X_tilde_it %*% t(X_tilde_it)
      assert_dimensions(XX, dim_=c(p, p))
      sum_XX = sum_XX + XX
      
      Xy = X_tilde_it * y_tilde_it
      assert_dimensions(Xy, dim_=p)
      sum_Xy = sum_Xy + Xy
    }
    
    B_period_t = solve(sum_XX) %*% sum_Xy
    assert_dimensions(B_period_t, dim_=c(p, 1))
    B_period_t = B_period_t[,1]
    
    assert_dimensions(B_period[,t], B_period_t)
    B_period[,t] = B_period_t
  }
  
  return(B_period)
}

estimate_slopes_post_selection = function(N, p, X, y, M_F_hat, group_hat, regime_hat, K_hat, m_hat, verbose=FALSE){
  A_hat = array(NA, dim=c(K_hat, m_hat+1, p))
  
  for (l in 1:K_hat){
    for (u in 1:(m_hat+1)){
      periods = which(regime_hat==u)
      num_periods = length(periods)
      
      sum_XMX = 0
      sum_XMy = 0
      num_obs = 0
      for (i in 1:N){
        if (group_hat[i] == l){
          X_iu = X[i,periods,]
          y_iu = y[i,periods]
          M_F_hat_iu = M_F_hat[periods,periods]
          
          # reshape vectors to conformable matrices when num_periods=1
          if (num_periods == 1){
            X_iu = matrix(data=X_iu, nrow=num_periods, ncol=p)
            M_F_hat_iu = matrix(data=M_F_hat_iu, nrow=num_periods, ncol=num_periods)
          }
          
          # check dimensions
          assert_dimensions(X_iu, dim_=c(num_periods, p))
          assert_dimensions(y_iu, dim_=num_periods)
          assert_dimensions(M_F_hat_iu, dim_=c(num_periods, num_periods))        

          sum_XMX_i = t(X_iu) %*% M_F_hat_iu %*% X_iu
          assert_dimensions(sum_XMX_i, dim_=c(p, p))
          sum_XMX = sum_XMX + sum_XMX_i
          
          sum_XMy_i = t(X_iu) %*% M_F_hat_iu %*% y_iu
          assert_dimensions(sum_XMy_i, dim_=c(p, 1))
          sum_XMy = sum_XMy + sum_XMy_i
          
          num_obs = num_obs + num_periods
        }
      }
      
      # Use the Moore-Penrose inverse to handle the cases when sum_XMX is singular
      A_hat_lu = MASS::ginv(sum_XMX) %*% sum_XMy
      assert_dimensions(A_hat_lu, dim_=c(p, 1))
      A_hat_lu = A_hat_lu[,1]

      assert_dimensions(A_hat[l,u,], A_hat_lu)
      A_hat[l,u,] = A_hat_lu
      if (verbose) clog(paste0("l=", l, ", u=", u, ", A_hat_lu=", vec2str(round(A_hat_lu,2)), ", num_obs=", num_obs))
    }
  }
  
  return(A_hat)
}

# Time the functions
tune_gamma = timer(tune_gamma)
tune_gamma_max = timer(tune_gamma_max)
estimate_slopes = timer(estimate_slopes)
