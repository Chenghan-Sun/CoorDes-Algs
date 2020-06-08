# This .R code file consists of:
  # Algorithm 1: Coordinate Descent method with randomized / cyclic rules and fixed step size 
  # for solving quadratic form objective function

# Arthurs: STA 243 Final Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun

library(pracma)

# Coordinate Descent method
RCDM = function(A, b, xs, xk = NULL, cr = NULL, iter_k = NULL, 
                alpha=0.001, tol=10^-2, maxIter=10^7, rule="random") {
  # Params: 
    # 
  
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
  fx = c(quadratic_obj(A%*%xk, b))
  fstar = quadratic_obj(A%*%xs, b)  # true function value 
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
  while (abs(fx[k] - fstar) >= tol) {
    # update the gradient
    u = A%*%xk - b
    A1 = t(A)
    gd_k = A1[iter_k, ]%*%u
    
    # update xk
    xk[iter_k] = xk[iter_k] - alpha*gd_k
    
    # update criterion cr
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    if (mod(k, 1000) == 0) {
      print(c(paste("step", k), paste("error", cr[k+1])))
    }
    
    # update k
    k = k+1
    
    # update estimated function value and error 
    fx = c(fx, quadratic_obj(A%*%xk, b))
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

RCDM_results = RCDM(A, b, xs, alpha = 1, tol = 0.005)
RCDM_results$k

#gap vs iteration 
plot(RCDM_results$cr)

#eps vs nums of iteration 
set.seed(100)
N = 30
eps = seq(0.005, 0.1, length.out = N)
num_iter = sapply(eps, function(eps){
  RCDM_results = RCDM(A, b, xs, alpha = 1, tol = eps)
  RCDM_results$k
})
plot(log(1 / eps), num_iter)

#sigma vs nums of iteration 
set.seed(100)
N = 100
kappa = seq(1, 30, length.out = N)
num_iter = sapply(kappa, function(kappa){
  m = 100
  n = 50
  k = kappa

  u = randortho(m)  
  v = randortho(n)
  s_c_diag  = seq(from = 1 / sqrt(k), to = 1, length.out = min(m, n))
  s_c = diag(s_c_diag, nrow=m, ncol=n)  # for convexity assumption 
  
  A = u%*%s_c%*%v  
  
  xs  = ones(n, 1)
  b = A%*%xs + 1 / (1 * 1000) * rnorm(m)
  
  RCDM_results = RCDM(A, b, xs, alpha = 1, tol = 0.01)
  RCDM_results$k
})
plot(kappa, num_iter)



### seperable CD ###
SpCD <- function(A, b, xs, lambda = 1, iter_k = 1, xk = NULL, cr = NULL, 
                 alpha = 0.001, tol = 1e-2, maxIter = 1e7){
  # set k as the counter
  # CGD method terminates when norm(xk-xs)/norm(xs) smaller than given epsi = 10^-3
  # denote norm(xk-xs)/norm(xs) = cr (criterion)
  k = 1
  cr = c(1)
  
  # initialize x 
  if (is.null(xk)){
    xk = zeros(n, 1)
  }
  #gradient 
  gd_k = 0 
  
  
  # Define the objective function
  quadratic_obj = function(xk, y){
    fun_val = 0.5*norm((y - A%*%xk), "2")^2 + lambda * sum(abs(xk))
    return(fun_val)
  }
  
  fx = c(quadratic_obj(xk, b))
  fstar  =  quadratic_obj(xs, b)
  error = c()
  
  while (abs(fx[k] - fstar) >= tol) {
    # update the gradient
    u = A%*%xk - b
    A1 = t(A)
    gd_k = A1[iter_k, ]%*%u
    
    #### update xk
    sub.optim <- function(x){
      gd_k * (x - xk[iter_k]) + 1 / (2 * alpha) * (x - xk[iter_k]) ^ 2 + lambda * abs(x)
    }
    sub.optim.result = optim(0, sub.optim, method = "BFGS")
    z_k = sub.optim.result$par
    xk[iter_k] = z_k
    #print(z_k)
    
    # update stopping criterion 
    cr[k+1] = norm(xk-xs, "2") / norm(xs, "2")
    
    if (mod(k, 1000) == 0) {
      #print(c(paste("step", k),paste("error", cr[k+1]) ))
      print(c(paste("step", k), paste("error", fx[k] - fstar) ))
      #print(gd[iter_k])
      #print(xk)
      #print(z_k)
      #print(xk)
    }
    
    
    # update error 
    fx = c(fx, quadratic_obj(xk, b))
    error = c(error, norm((xk - xs), "2"))
    
    # update k
    k = k+1
    
    # update iter_k
    iter_k = mod(iter_k, n) + 1
    
    if (k > maxIter) {
      print(paste("Algorithm unfinished by reaching the maximum iterations."))
      break
    }
  }
  print(xk)
  return(list(k = k, cr = cr, error = error, fx = fx ))
}


### Experiment ### 
# input sparse data points, here xs is the true solution that we want to find
m = 100
n = 200
k = 50
s = 30


u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
v = randortho(n)

#s_c_diag = runif(min(m, n), min= 1 / sqrt(k), max = 1)
s_c_diag  = seq(from = 1 / sqrt(k), to = 1, length.out = min(m, n))
s_c = diag(s_c_diag, nrow=m, ncol=n)  # for convexity assumption 

# sigular value decomposition
A = u%*%s_c%*%v  # for convexity assumption

#xs = rnorm(n)
xs  = zeros(n, 1)
xs[1:s,] = 1
b = A %*% xs 
#solve(t(A) %*% A, t(A) %*% b)
# t(A) %*% A 

SpCD_results = SpCD(A, b, xs, lambda = 0.01, alpha = 1, tol = 0.0001)
RCDM_results = RCDM(A, b, xs, alpha = 1, tol = 0.005)
SpCD_results$k
RCDM_results$k

plot(SpCD_results$fx)
plot(RCDM_results$fx)


