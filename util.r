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

vec2str = function(vec){
  string = paste(vec, collapse = " ")
  return(string)
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


