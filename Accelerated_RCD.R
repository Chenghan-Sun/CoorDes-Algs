# This .R code file consists of:
  # Algorithm 2: Accelerated Randomized Coordinate Descent (Nesterov 2012)
  # for solving quadratic form objective function

# Arthurs: STA 243 Final Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun

library(pracma)


# Coordinate Descent method
A_RCDM = function(A, b, xs, sigma_max, xk = NULL, cr = NULL, iter_k = NULL, 
                  mu_k=NULL, gamma_k=0, tol=10^-2, maxIter=10^7) {
  # Params: 
  # 
  
  n = ncol(A)
  
  # set k as the counter
  k = 1
  
  # denote cr = norm(xk-xs)/norm(xs) as the a criterion
  cr = c(1)
  
  # initialize xk as the estimates' vector 
  if (is.null(xk)){
    xk = zeros(n, 1)
  }
  
  # initialize mu_k as a vector of sub-direction 
  if (is.null(mu_k)){
    mu_k = zeros(n, 1)
  }
  
  # Define the objective function
  quadratic_obj = function(y, yhat){
    fun_val = 0.5*norm((y - yhat), "2")^2
    return(fun_val)
  }
  
  # initialize objective function vector
  fx = c(quadratic_obj(A%*%xk, b))
  fstar = quadratic_obj(A%*%xs, b)  # true function value 
  error = c()  # initialize error vector 
  
  large_root_finder = function(a, b, c) {
    sol1 = (-b - sqrt(b^2 - 4*a*c)) / (2*a)
    sol2 = (-b + sqrt(b^2 - 4*a*c)) / (2*a)
    return(max(sol1, sol2))
  }
  
  # main loop 
  while (abs(fx[k] - fstar) >= tol) {
    # choose gamma_k to be the larger root 
    a = 1
    b = (sigma_max*gamma_k^2/n) - (1/n)
    c = -gamma_k^2
    gamma_k_new = large_root_finder(a, b, c)
    
    # set alpha_k and beta_k
    alpha_k = (n - gamma_k_new*sigma_max) / (gamma_k_new*(n^2 - sigma_max))
    beta_k = 1 - gamma_k_new*sigma_max/n
    
    # set y_k
    y_k = alpha_k*mu_k + (1-alpha_k)*xk
    
    # choose index iter_k with uniform probability
    iter_k = sample(n, 1)
    
    # set dk
    u1 = alpha_k*A%*%mu_k + (1 - alpha_k)*A%*%xk - b
    A1 = t(A)
    dk = (1 - alpha_k)*A1[iter_k, ]%*%u1
    
    # update xk
    xk[iter_k] = y_k[iter_k] - (1/sigma_max^2)*dk
    
    # update criterion cr
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    if (mod(k, 1000) == 0) {
      print(c(paste("step", k), paste("error", cr[k+1])))
    }
    
    # update mu_k
    #print(beta_k*mu_k)
    #print((1 - beta_k)*y_k)
    #print((gamma_k_new/sigma_max^2)*dk)
    mu_k[iter_k] = beta_k*mu_k[iter_k] + (1 - beta_k)*y_k[iter_k] - (gamma_k_new/sigma_max^2)*dk
    
    # update gamma
    gamma_k = gamma_k_new
    
    # update k
    k = k+1
    
    # update estimated function value and error 
    fx = c(fx, quadratic_obj(A%*%xk, b))
    error = c(error, norm((xk - xs), "2"))
    
    # set algorithm iter bound 
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  return(list(k = k, cr = cr, error = error, fx = fx))
}

### Experiment ### 
# input data points, here xs is the true solution that we want to find
m = 100
n = 50
k = 30


u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
v = randortho(n)

#s_c_diag = runif(min(m, n), min= 1 / sqrt(k), max = 1)
s_c_diag  = seq(from = 1 / sqrt(k), to = 1, length.out = min(m, n))
s_c = diag(s_c_diag, nrow=m, ncol=n)  # for convexity assumption 

# sigular value decomposition
A = u%*%s_c%*%v  # for convexity assumption

#xs = rnorm(n)
xs  = ones(n, 1)
b = A%*%xs + 1 / (1 * 1000) * rnorm(m)
# solve(t(A) %*% A, t(A) %*% b)
# t(A) %*% A 

A_RCDM_results = A_RCDM(A, b, xs, sigma_max=max(s_c_diag), tol = 0.0001)
print(paste("The total number of iteration for ARCD algorithm = ", A_RCDM_results$k))

#gap vs iteration 
plot(A_RCDM_results$cr)



