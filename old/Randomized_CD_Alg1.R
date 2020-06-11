# This .R code file consists of:
  # Algorithmc1: Coordinate Descent method with randomized / cyclic rules and fixed step size 
  # for solving Lasso Problem 
# Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun

library(pracma)

### Coordinate Descent method ###
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
    xk[iter_k] = xk[iter_k] - alpha*gd[iter_k]
    # print(xk[iter_k])
    
    # update stopping criterion 
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    if (mod(k, 1000) == 0) {
      print(cr[k+1])
    }
    
    # update k
    k = k+1
    
    # update error 
    fx = c(fx, quadratic_obj(A%*%xk, b))
    # fx = c(fx, 1/2*(norm(A%*%xk-b, "2"))^2)
    error = c(error, norm((xk - xs), "2"))
    #print(error[k])
    
    # update iter_k
    iter_k = mod(iter_k, n)+1
    
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k = k, cr = cr, error = error, xk = xk))
}



### Experiment ### 
# input data points, here xs is the true solution that we want to find
m = 300
n = 100
k = 5


u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
v = randortho(n)

s_c_diag = runif(min(m, n), min= 1 / sqrt(k), max = 1)
s_c = diag(s_c_diag, nrow=m, ncol=n)  # for convexity assumption 

# sigular value decomposition
A = u%*%s_c%*%v  # for convexity assumption

xs = rnorm(n)
b = A%*%xs + 1 / (k) * rnorm(m)

xs = zeros(n, 1)
picks = randperm(n)
xs[picks[1:s]] = 100*randn(s, 1)
b = A%*%xs + 5 * rnorm(m)



RCDM_results = RCDM(A, b, xk, xs, cr, iter_k, alpha = 1 / 1000)

beta_hat = RCDM_results$xk
norm(beta_hat - xs, "2")



