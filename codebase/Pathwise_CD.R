# This .R code file consists of:
# Algorithm 4: Pathwise Coordinate descent for the lasso
# for solving quadratic form objective function with L1 penalty

# Arthurs: STA 243 Final Project Group Members:
# Han Chen, Ninghui Li, Chenghan Sun


Path_RCDM = function(A, b, xs, xk = NULL, lamb=1, alpha=0.001, 
                     tol=10^-3, maxIter=10^7, rule="cyclic") {
  ### Solve quadratic form functions min_x f(x) = (Ax - b) ^ 2 + lambda|x| ###
  ### Algorithm : Pathwise Coordinate descent for the lasso "PATHWISE COORDINATE OPTIMIZATION"
  ### A the input matrix, b vector, xs the true parameters 
  ### lamb the tuning parameter
  ### alpha : stepsize for each step coordinate descent.
  ### xk initial value of the optimization problem 
  ### stopping criterion f(xk) - fstar < tol, where fstar = f(xs), we stop the function if the iteration exceed maxIter..
  
  # initialize k 
  k = 1
  #initialize cr
  cr = c(1)
  # iter_k = ik
  iter_k = 1
  #fxout 
  fxopt = 1/2*(norm((A%*%xs - b), "2"))^2 + lamb*norm(xs, "1")
  
  
  # initialize x 
  if(is.null(xk)){
    xk = zeros(n, 1)
  }
  # initialize gradient vector
  gd = zeros(n, 1)
  
  fx = c(1/2*(norm(A%*%xk-b, "2"))^2 + lamb*norm(xk, "1"))
  error = c(fx[k] - fxopt)
  while (error[k] >= tol) {
    # update the gradient
    u = A%*%xk - b
    A1 = t(A)
    gd1 = A1[iter_k, ]%*%u
    # gd2 = shrink_operator(xk[iter_k], lamb)
    gd2 = sign(xk[iter_k])%*%lamb
    gd[iter_k] = gd1 + gd2
    
    # update xk
    xk[iter_k] = xk[iter_k] - alpha*gd[iter_k]
    
    # update stopping criterion 
    cr[k+1] = norm((xk-xs), "2") / norm(xs, "2")
    
    if (mod(k, 1000) == 0) {
      print(error[k])
    }
    
    # update k
    k = k+1
    
    # update error 
    # fx = c(fx, lasso_obj(A%*%xk, b, xk, lamb))
    fx = c(fx, 1/2*(norm(A%*%xk-b, "2"))^2 + lamb*norm(xk, "1"))
    error = c(error, fx[k] - fxopt)
    
    # update iter_k
    iter_k = mod(iter_k, n)+1
    
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k = k, cr = cr, error = error))
}




