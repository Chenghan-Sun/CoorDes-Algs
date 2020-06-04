# This .R code file consists of:
  # Algorithmc1: Coordinate Descent method with randomized / cyclic rules and fixed step size 
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

# Define the objective function
quadratic_obj = function(y, yhat){
  fun_val = 0.5*norm((y - yhat), "2")^2
  return(fun_val)
}

# Experiment Setup

# initialize fixed CGD step size alpha
alpha = 0.001

# CGD method terminates when norm(xk-xs)/norm(xs) smaller than given epsi = 10^-3
# denote norm(xk-xs)/norm(xs) = cr (criterion)
# set k as the counter
k = 1
cr = c(1)

# Initialize index iter_k (apply for both cyclic and randomized rules)
iter_k = 1 

# Coordinate Descent method
RCDM = function(A, b, xk, xs, cr, iter_k, 
                alpha=0.001, tol=10^-2, maxIter=10^7, rule="cyclic") {
  
  # initialize x 
  xk = zeros(n, 1)
  
  # initialize gradient vector
  gd =  zeros(n, 1)
  
  fx = c()
  error = c()
  while (cr[k] >= tol) {
    # update the gradient
    u = A%*%xk - b
    A1 = t(A)
    gd_k = A1[iter_k, ]%*%u
    gd[iter_k] = gd_k
    
    # update xk
    xk[iter_k] = xk[iter_k] - alpha%*%gd[iter_k]
    # print(xk[iter_k])
    
    # update stopping criterion 
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    if (mod(k, 1000) == 0) {
      print(cr[k+1])
    }
    
    # update error 
    fx = c(fx, quadratic_obj(A%*%xk, b))
    error = c(error, norm((xk - xs), "2"))
    #print(error[k])
    
    # update k
    k = k+1
    
    # update iter_k
    iter_k = mod(iter_k, n)+1
    
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k, cr, error))
}

RCDM_results = RCDM(A, b, xk, xs, cr, iter_k)






