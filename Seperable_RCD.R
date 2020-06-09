

##### seperable CD #####

SpCD <- function(A, b, xs, lambda = 1, iter_k = 1, xk = NULL, cr = NULL, 
                 alpha = 0.001, tol = 1e-2, maxIter = 1e7){
  ### Solve quadratic form functions with L1 penalty min_x f(x) = (Ax - b) ^ 2 + lambda|x| ###
  ### Algorithm: Seperable Regularised version Coordinate Descent "Richtarik, P., Takac, M.: Iteration complexity of a randomized block-coordinate descent
  ### methods for minimizing a composite function"
  ### A the input matrix, b vector, xs the true parameters 
  ### lambda the tuning parameter
  ### alpha : usually set as 1 / L_max, where L_max is the maximum of component Lipschitz constants.
  ### xk initial value of the optimization problem 
  ### stopping criterion f(xk) - fstar < tol, where fstar = f(xs), we stop the function if the iteration exceed maxIter.
  
  
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
  
  while (fx[k] - fstar >= tol) {
    # update the gradient
    u = A%*%xk - b
    A1 = t(A)
    gd_k = A1[iter_k, ]%*%u
    
    # update xk
    ## Here we use the soft_threshold to solve the suboptimization problem 
    a_k = 1 / alpha 
    c_k = 1 / alpha * xk[iter_k] - gd_k 
    if(c_k < -1 * lambda){
      z_k = (c_k + lambda) / a_k 
    }else if (c_k > lambda){
      z_k = (c_k - lambda) / a_k 
    }else{
      z_k = 0
    }
    
    
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


#####  Experiment  #####

# input sparse data points, here xs is the true solution that we want to find
m = 100
n = 50
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
#RCDM_results = RCDM(A, b, xs, alpha = 1, tol = 0.005)
SpCD_results$k
#RCDM_results$k

plot(SpCD_results$fx)
#plot(RCDM_results$fx)


### gap vs iteration ###
plot(SpCD_results$cr)




### eps vs nums of iteration  (strong convex) ###
set.seed(300)
m = 100
n = 50
k = 50
s = 30

u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
v = randortho(n)

#s_c_diag = runif(min(m, n), min= 1 / sqrt(k), max = 1)
s_c_diag  = seq(from = 1 / sqrt(k), to = 1, length.out = min(m, n))
s_c = diag(s_c_diag, nrow=m, ncol=n)  # for convexity assumption 

# sigular value decomposition
A = u%*%s_c%*%v  # for convexity assumption

xs  = zeros(n, 1)
xs[1:s,] = 1
b = A %*% xs 

#produce eps with length N
N = 100
eps = seq(0.005, 0.1, length.out = N)
num_iter = sapply(eps, function(eps){
  SpCD_results = SpCD(A, b, xs, lambda = 0.01, alpha = 1, tol = eps)
  SpCD_results$k
})
plot(log(1 / eps), num_iter)





### sigma vs nums of iteration (strong convex) ### 
set.seed(100)
N = 100
kappa = seq(1, 30, length.out = N)
num_iter = sapply(kappa, function(kappa){
  m = 100
  n = 50
  k = kappa
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
  
  SpCD_results = SpCD(A, b, xs, lambda = 0.01, alpha = 1, tol = eps)
  SpCD_results$k
})

plot(kappa, num_iter)








### eps vs nums of iteration  (strong convex) ###
set.seed(300)
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

xs  = zeros(n, 1)
xs[1:s,] = 1
b = A %*% xs 

#produce eps with length N
N = 100
eps = seq(0.005, 0.1, length.out = N)
num_iter = sapply(eps, function(eps){
  SpCD_results = SpCD(A, b, xs, lambda = 0.01, alpha = 1, tol = eps)
  SpCD_results$k
})
plot(log(1 / eps), num_iter)







