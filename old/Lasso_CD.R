# This .R code file consists of:
  # Algorithm: Coordinate Descent method with cyclic rule and fixed step size 
  # for solving Lasso Problem 
# Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun

library(pracma)
# Simulation dataset
# input data points, here xs is the true solution that we want to find
m = 300
n = 500
s = 2
A = randn(m, n)
xs = zeros(n, 1)
picks = randperm(n)
xs[picks[1:s]] = 100*randn(s, 1)
b = A%*%xs

# Define the objective function: Lasso
lasso_obj = function(y, yhat, xk, lamb){
  # 1 / (2 * nobs) RSS + lambda * penalty
  # fun_val = 0.5 * sum((y-yhat)^2) + lam*sum(abs(betat))
  fun_val = 0.5*norm((y - yhat), "2")^2 + lamb*norm(xk, "1")
  return(fun_val)
}

# Closed form solution of lasso problem: Soft-Threshold
shrink_operator = function(b, lamb) {
  flag = max(0, b-lamb) - max(0, -b-lamb)
  return(flag)
}

# Experiment Setup

# initialize fixed CGD step size alpha
alpha = 0.001

# CGD method terminates when norm(xk-xs)/norm(xs) smaller than given epsi = 10^-3
# denote norm(xk-xs)/norm(xs) = cr (criterion)
# set k as the counter
k = 1
cr = c(1)
lamb=1
fxopt = 1/2*(norm((A%*%xs - b), "2"))^2 + lamb*norm(xs, "1")

# Initialize index iter_k (apply for both cyclic and randomized rules)
iter_k = 1 

# Coordinate Descent method
RCDM = function(A, b, xk, xs, cr, fxopt, iter_k, 
                lamb=1, alpha=0.001, tol=10^-3, maxIter=10^7, rule="cyclic") {
  # iter_k = ik
  
  # initialize x 
  xk = zeros(n, 1)
  
  # initialize gradient vector
  gd = zeros(n, 1)
  
  fx = c()
  error = c()
  while (cr[k] >= tol) {
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
      print(cr[k+1])
    }
    
    # update k
    k = k+1
    
    # update error 
    # fx = c(fx, lasso_obj(A%*%xk, b, xk, lamb))
    fx = c(fx, 1/2*(norm(A%*%xk-b, "2"))^2 + lamb*norm(xk, "1"))
    error = c(error, abs(fx[k] - fxopt))
    
    # update iter_k
    iter_k = mod(iter_k, n)+1
    
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k, cr, error))
}

RCDM_results = RCDM(A, b, xk, xs, cr, fxopt, iter_k)

