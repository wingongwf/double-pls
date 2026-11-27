# This module consists of all utility functions.

get_circular_index = function(pos, step, size){
  new = pos + step
  if (new<1){
    new = new + size
  } else if (new>size){
    new = new - size
  }
  
  return(new)
}

get_gamma_grid = function(gamma_max, gamma_min_ratio, num_gamma){
  gamma_1_max = gamma_max[1]
  gamma_2_max = gamma_max[2]
  gamma_1_min = gamma_1_max * gamma_min_ratio
  gamma_2_min = gamma_2_max * gamma_min_ratio
  
  gamma_1_grid = 10 ^ seq(from=log10(gamma_1_min), to=log10(gamma_1_max), length.out=num_gamma)
  gamma_2_grid = 10 ^ seq(from=log10(gamma_2_min), to=log10(gamma_2_max), length.out=num_gamma)
  
  gamma_grid = list(gamma_1=gamma_1_grid, gamma_2=gamma_2_grid)
  return(gamma_grid)
}

# Get slopes for unit i at all times
get_B_i = function(B, i, T_, p){
  index_start = (i-1) * T_ + 1
  index_end = i * T_
  B_i = B[,index_start:index_end]
  
  assert_dimensions(B_i, dim_=c(p, T_))
  return(B_i)
}

# Get slopes for all units at time t
get_B_t = function(B, t, T_, N, p){
  cols = seq(from=t, to=t+(N-1)*T_, by=T_)
  B_t = B[,cols]
  
  assert_dimensions(B_t, dim_=c(p, N))
  return(B_t)
}

is_invertible = function(A){
  # Check if a square matrix A is invertible or not
  num_row = dim(A)[1]
  num_col = dim(A)[2]
  stopifnot(num_row==num_col)
  
  return(qr(A)$rank == nrow(A))
}

B02vec = function(B0, N, T_, p){
  B0_vec = rep(NA, N*T_*p)
  for (i in 1:N){
    for (t in 1:T_){
      B0_it = B0[i,t,]
      assert_dimensions(B0_it, dim_=p)
      
      start_index = (i-1)*T_*p + (t-1)*p + 1
      end_index = start_index + p - 1
      assert_dimensions(B0_vec[start_index:end_index], B0_it)
      B0_vec[start_index:end_index] = B0_it
    }
  }
  
  stopifnot(sum(is.na(B0_vec)) == 0)
  return(B0_vec)
}

B02matrix = function(B0, N, T_, p){
  B0_matrix = matrix(NA, p, N*T_)
  for (i in 1:N){
    B0_i = t(B0[i,,])
    assert_dimensions(B0_i, dim_=c(p, T_))
    
    start_index = (i-1) * T_ + 1
    end_index = start_index + T_ - 1 
    assert_dimensions(B0_matrix[,start_index:end_index], B0_i)
    B0_matrix[,start_index:end_index] = B0_i
  }
  
  stopifnot(sum(is.na(B0_matrix)) == 0)
  return(B0_matrix)
}

Btilde2vec = function(A_hat, group_hat, regime_hat, N, T_, p){
  B_tilde_vec = rep(NA, N*T_*p)
  for (i in 1:N){
    for (t in 1:T_){
      group_i = group_hat[i]
      regime_t = regime_hat[t]
      B_tilde_it = A_hat[group_i,regime_t,]
      assert_dimensions(B_tilde_it, dim_=p)
      
      start_index = (i-1)*T_*p + (t-1)*p + 1
      end_index = start_index + p - 1
      assert_dimensions(B_tilde_vec[start_index:end_index], B_tilde_it)
      B_tilde_vec[start_index:end_index] = B_tilde_it
    }
  }
  
  stopifnot(sum(is.na(B_tilde_vec)) == 0)
  return(B_tilde_vec)
}

vec2str = function(vec){
  string = paste(vec, collapse = " ")
  return(string)
}

list2str = function(x, digits=NULL){
  out = ""
  for (item in names(x)){
    val = ifelse(is.null(digits), x[[item]], round(x[[item]], digits))
    out = paste0(out, item, "=", val, ", ")
  }
  
  out = substring(out, 1, nchar(out)-2) # remove trailing comma
  return(out)
}

eval_perf = function(N, T_, p, B0, K0, m0, fit, export=TRUE){
  A_hat = fit$A_hat
  group_hat = fit$group_hat
  regime_hat = fit$regime_hat
  K_hat = fit$K_hat
  m_hat = fit$m_hat
  
  # Reshape true slopes (B0) and estimated slopes (B_tilde) to vectors
  # B0_vec: (N*T_*p) x 1
  # B_tilde_vec: (N*T_*p) x 1
  B0_vec = B02vec(B0=B0, N=N, T_=T_, p=p)
  B_tilde_vec = Btilde2vec(A_hat=A_hat, group_hat=group_hat, regime_hat=regime_hat, N=N, T_=T_, p=p)
  assert_dimensions(B0_vec, B_tilde_vec)
  
  # Export B_tilde_vec to csv
  if (export){
    csvline = paste0(paste(B_tilde_vec, collapse=","), "\n")
    cat(csvline, file=CSV_BETA, append=TRUE)    
  }

  B0_matrix = B02matrix(B0=B0, N=N, T_=T_, p=p)
  sum_squared_errors = sum((B_tilde_vec - B0_vec)^2)
  group0 = estimate_groups(N=N, T_=T_, p=p, B=B0_matrix)
  regime0 = estimate_regimes(N=N, T_=T_, p=p, B=B0_matrix)
  stopifnot(max(group0) == K0)
  stopifnot(max(regime0) == m0+1)
  
  metrics = list(
    correct_num_groups=(K_hat == K0),
    correct_membership=all(group_hat == group0),
    correct_num_breaks=(m_hat == m0),
    correct_break_dates=all(regime_hat == regime0),
    num_groups=K_hat,
    num_breaks=m_hat,
    sum_squared_errors=sum_squared_errors,
    size=length(B0_vec)
  )
  clog(paste0("performance: ", list2str(metrics,4)))
  return(metrics)
}

make_csv_Btilde = function(N, T_, p, filename){
  header = ""
  for (i in 1:N){
    for (t in 1:T_){
      for (k in 1:p){
        header = paste0(header, "Btilde_", i, "_", t, "_", k, ",")
      }
    }
  }
  
  header = paste0(substr(header, 1, nchar(header)-1), "\n")
  cat(header, file=filename, append=TRUE)
}

read_args = function(args){
  args = sub("^--", "", args)
  params = strsplit(args, "=", fixed=TRUE)
  params = setNames(lapply(params, `[`, 2), sapply(params, `[`, 1))
  
  out = list(
    # parameters for simulation
    n_rep=as.integer(params$n_rep),
    
    # parameters for data generating process
    N=as.integer(params$N),
    T_=as.integer(params$T_),
    p=as.integer(params$p),
    r=as.integer(params$r),
    phi=as.numeric(params$phi),
    K0=as.integer(params$K0),
    m0=as.integer(params$m0),
    pi_=as.numeric(params$pi_),
    L=as.integer(params$L),
    kappa=as.integer(params$kappa),
    
    # parameters for tuning gamma
    gamma_max_start_value=as.numeric(params$gamma_max_start_value),
    gamma_min_ratio=as.numeric(params$gamma_min_ratio),
    num_gamma=as.integer(params$num_gamma)
  )
  
  return(out)
}

timer = function(f){
  # Can only log scalar and named parameters
  fn_name = deparse(substitute(f))
  force(f)
  function(...){
    start = Sys.time()
    fn_params = list(...)
    clog(paste0(fn_name, ": status=start, ", list2str(fn_params)))
    output = f(...)
    
    elapsed = difftime(Sys.time(), start, units="secs")
    clog(paste0(fn_name, ": status=end, elapsed=", round(elapsed, 2), ", ", list2str(fn_params)))
    return(output)
  }
}

clog = function(msg, file=LOG_FILE){
  logline = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%OS6 "), msg, "\n")
  cat(logline)   # print to console
  cat(logline, file=file, append=TRUE)    # print to log
}

assert_dimensions = function(A, B=NULL, dim_=NULL){
  dim_A = if (is.array(A) || inherits(A, "CVXR::Expression")) dim(A) else length(A)
  
  if (!is.null(B) && is.null(dim_)){
    dim_B = if (is.array(B) || inherits(B, "CVXR::Expression")) dim(B) else length(B)
    stopifnot(identical(dim_A, dim_B))
    
  } else if (is.null(B) && (!is.null(dim_))){
    dim_ = as.integer(dim_)
    stopifnot(identical(dim_A, dim_))
    
  } else{
    stop("One and only one of B or dim_ should be supplied. Or some other error.")
  }
}

# # Some test cases
# # Good
# assert_dimensions(matrix(NA,2,3), matrix(NA,2,3))
# assert_dimensions(matrix(NA,2,3), dim_=c(2,3))
# assert_dimensions(array(NA, c(1,2,3)), array(NA, c(1,2,3)))
# assert_dimensions(array(NA, c(1,2,3)), dim_=c(1,2,3))
# assert_dimensions(c(NA,NA), c(3,5))
# assert_dimensions(c(NA,NA), dim_=2)
# 
# # Bad
# assert_dimensions(matrix(NA,2,3), matrix(NA,3,3))
# assert_dimensions(matrix(NA,2,3), dim_=c(3,3))
# assert_dimensions(matrix(NA,2,3), matrix(NA,3,3), c(1,2))
# assert_dimensions(matrix(NA,2,3))
# assert_dimensions(array(NA, c(1,2,3)), array(NA, c(1,2,4)))
# assert_dimensions(array(NA, c(1,2,3)), dim_=c(1,2,4))
# assert_dimensions(c(1,2), c(3,4,5))
# assert_dimensions(c(NA,NA), dim_=3)


