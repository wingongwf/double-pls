# This module consists of all functions for running simulations

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

run_simulation = function(n_rep, N, T_, p, r, phi, K0, m0, pi_, L, kappa){
  make_csv_Btilde(N=N, T_=T_, p=p, filename=CSV_BETA)
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
  
  # run n_rep times
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
  return(metrics_agg)
}

# Time the functions
run_replication = timer(run_replication)
run_simulation = timer(run_simulation)


