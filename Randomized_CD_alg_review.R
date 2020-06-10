# This .R code file consists of:
  # Algorithm 1: Coordinate Descent method with randomized / cyclic rules and fixed step size 
  # for solving quadratic form objective function

# Arthurs: STA 243 Final Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun


##### Coordinate Descent method #####
RCDM = function(A, b, xs, xk = NULL, cr = NULL, iter_k = NULL, 
                alpha=0.001, tol=10^-2, maxIter=10^7, rule="random") {
  # Params: 
  # columns of A
  n = ncol(A)
  # set k as the counter
  k = 1
  
  # denote cr = norm(xk-xs)/norm(xs) as the a criterion
  cr = c(1)
  
  # initialize xk as the estimates' vector 
  if (is.null(xk)){
    xk = zeros(n, 1)
  }
  
  # Define the objective function
  quadratic_obj = function(y, yhat){
    fun_val = 0.5*norm((y - yhat), "2")^2
    return(fun_val)
  }
  
  # initialize objective function vector
  fstar = quadratic_obj(A%*%xs, b)  # true function value 
  fx = c(quadratic_obj(A%*%xk, b) - fstar)
  error = c()  # initialize error vector 
  
  if (rule == "random") {
    iter_k = sample(ncol(A), 1)
  } else if (rule == "cyclic") {
    iter_k = 1
  } else {
    print(paste("Need to specify variants of CD method."))
    break
  }
  
  # main loop 
  while (fx[k] >= tol) {
    # update the gradient
    u = A%*%xk - b
    A1 = t(A)
    gd_k = A1[iter_k, ]%*%u
    
    # update xk
    xk[iter_k] = xk[iter_k] - alpha*gd_k
    
    # update criterion cr
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    # if (mod(k, 1000) == 0) {
    #   print(c(paste("step", k), paste("error", cr[k+1])))
    # }
    
    # update k
    k = k+1
    
    # update estimated function value and error 
    fx = c(fx, quadratic_obj(A%*%xk, b) - fstar)
    error = c(error, norm((xk - xs), "2"))
    
    # update iter_k based on random / cyclic rules
    if (rule == "cyclic") {
      iter_k = mod(iter_k, n) + 1
    } else if (rule == "random") {
      iter_k = sample(ncol(A), 1)
    }
    
    # set algorithm iter bound 
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k = k, cr = cr, error = error, fx = fx))
}




#### Gradient Descent ####

GD = function(A, b, xs, xk = NULL, cr = NULL, iter_k = NULL, 
                alpha=0.001, tol=10^-2, maxIter=10^7, rule="random") {
  # Params: 
  # columns of A
  n = ncol(A)
  # set k as the counter
  k = 1
  # denote cr = norm(xk-xs)/norm(xs) as the a criterion
  cr = c(1)
  
  # initialize xk as the estimates' vector 
  if (is.null(xk)){
    xk = zeros(n, 1)
  }
  
  # Define the objective function
  quadratic_obj = function(y, yhat){
    fun_val = 0.5*norm((y - yhat), "2")^2
    return(fun_val)
  }
  
  # initialize objective function vector
  fstar = quadratic_obj(A%*%xs, b)  # true function value 
  fx = c(quadratic_obj(A%*%xk, b) - fstar)
  error = c()  # initialize error vector 
  
  if (rule == "random") {
    iter_k = sample(ncol(A), 1)
  } else if (rule == "cyclic") {
    iter_k = 1
  } else {
    print(paste("Need to specify variants of CD method."))
    break
  }
  
  # main loop 
  while ( fx[k] >= tol) {
    # update the gradient
    u = A%*%xk - b
    gd_k = t(A)%*%u
    
    # update xk
    xk = xk - alpha * gd_k

    # update criterion cr
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    # if (mod(k, 1000) == 0) {
    #   print(c(paste("step", k), paste("error", cr[k+1])))
    # }
    
    # update k
    k = k+1
    
    # update estimated function value and error 
    fx = c(fx, quadratic_obj(A%*%xk, b) - fstar)
    error = c(error, norm((xk - xs), "2"))
    

    
    # set algorithm iter bound 
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k = k, cr = cr, error = error, fx = fx))
}



