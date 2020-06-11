# This .R code file consists of:
  # Algorithm 3: Accelerated Randomized Coordinate Descent (Nesterov 2012)
  # for solving quadratic form objective function

# Arthurs: STA 243 Final Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun



# Coordinate Descent method
A_RCDM = function(A, b, xs, alpha = 1, Sigma = NULL, xk = NULL, iter_k = NULL, 
                   gamma_k=0, tol=10^-2, maxIter=10^7) {
  ### Solve quadratic form functions min_x f(x) = (Ax - b) ^ 2  ###
  ### Algorithm : Accelerated Coordinate Descent method "Nesterov, Y.: Effciency of coordinate descent methods on huge-scale optimization problems."
  ### A the input matrix, b vector, xs the true parameters 
  ### Sigma eigenvalues of A^TA, we use this to control L_i in the algorithm, 
  ### alpha fixed stepsize 
  ### xk initial value of the optimization problem 
  ### stopping criterion f(xk) - fstar < tol, where fstar = f(xs), we stop the function if the iteration exceed maxIter.
  
  # Params: 
  #dim of xs 
  n = ncol(A)
  
  # set k as the counter
  k = 1
  
  # calculate Sigma 
  if(is.null(Sigma)){
    Sigma = alpha * eigen(t(A) %*% A)$values
  }
  
  # denote cr = norm(xk-xs)/norm(xs) as the a criterion
  cr = c(1)
  
  # initialize xk as the estimates' vector 
  if (is.null(xk)){
    xk = zeros(n, 1)
  }
  
  # initialize mu_k as a vector of sub-direction 
  mu_k = xk
  
  # Define the objective function
  quadratic_obj = function(y, yhat){
    fun_val = 0.5*norm((y - yhat), "2")^2
    return(fun_val)
  }
  
  # initialize objective function vector
  fstar = quadratic_obj(A%*%xs, b)  # true function value 
  fx = c(quadratic_obj(A%*%xk, b) - fstar)
  error = c()  # initialize error vector 
  
  large_root_finder = function(a, b, c) {
    sol1 = (-b - sqrt(b^2 - 4*a*c)) / (2*a)
    sol2 = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
    return(max(sol1, sol2))
  }
  
  # main loop 
  while (fx[k] >= tol) {
    # choose gamma_k to be the larger root 
    para.a = 1
    para.b = (min(Sigma)*gamma_k^2/n) - (1/n)
    para.c = -gamma_k^2
    gamma_k_new = large_root_finder(para.a, para.b, para.c)
    
    # set alpha_k and beta_k
    alpha_k = (n - gamma_k_new*min(Sigma)) / (gamma_k_new*(n^2 - min(Sigma)))
    beta_k  =  1 - gamma_k_new*min(Sigma)  / n
    
    # set y_k
    y_k = alpha_k*mu_k + (1-alpha_k)*xk
    
    # choose index iter_k with uniform probability
    iter_k = sample(n, 1)
    
    # set dk
    # u1 = alpha_k*A%*%mu_k + (1 - alpha_k)*A%*%xk - b
    # A1 = t(A)
    # dk = (1 - alpha_k)*A1[iter_k, ]%*%u1
    
    dk = t(A[, iter_k]) %*% ( A %*% y_k - b ) 
    
    # update xk
    #xk[iter_k] = y_k[iter_k] - (1/sigma_max^2)*dk
    xk = y_k
    xk[iter_k] = xk[iter_k] - (1 / Sigma[iter_k]) * dk
    
    # update criterion cr
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    # if (mod(k, 1000) == 0) {
    #   print(c(paste("step", k), paste("error", cr[k+1])))
    #   print(paste("value of objective function", fx[k]))
    #   
    # }
    
    # update mu_k
    #print(beta_k*mu_k)
    #print((1 - beta_k)*y_k)
    #print((gamma_k_new/sigma_max^2)*dk)
    mu_k = beta_k * mu_k + (1 - beta_k) * y_k 
    mu_k[iter_k] = mu_k[iter_k] - (gamma_k_new/Sigma[iter_k]^2) * dk
    #mu_k[iter_k] = beta_k*mu_k[iter_k] + (1 - beta_k)*y_k[iter_k] - (gamma_k_new/Sigma[iter_k]^2)*dk
    
    # update gamma
    gamma_k = gamma_k_new
    
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
  return(list(k = k, cr = cr, error = error, fx = fx, xk = xk))
}




