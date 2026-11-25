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


