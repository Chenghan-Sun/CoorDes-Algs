# This .R code file consists of:
  # 1. Algorithms Convergence Analysis and plots
  # 2. Algorithms Complexity Analysis and plots
# For: 
  # 1. Algorithm 1: Randomized Coordinate Descent method
  # 2. Algorithm 2: Accelerated Randomized Coordinate Descent (Nesterov 2012)
  # 3. Algorithm 3: Seperable Coordinate Descent Algorithm
# Under:
  # 1. Strong convexity assumption 
  # 2. Convexity assumption 

# Arthurs: STA 243 Final Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun


library(pracma)
library(ggplot2)

# load all algorithms from codebase 
setwd("/Users/furinkazan/Box/STA_243/CoorDes-Algs/codebase/")  # please set your current directory path 
source("RCD.R")
source("Seperable_RCD.R")
source("Accelerated_RCD.R")
Quad_generator <- function(m = 100, n = 50, k = 30){
  ### problem set up    min_x |Ax - b|^2
  ### if m >= n, we generate data matrix A with A^TA has conditional number k, where the largest singualr value of A is 1
  ### if m < n , we generate data matrix A where the largest singular value of A is 1
  
  # input data points, here xs is the true solution that we want to find
  m = m
  n = n
  k = k
  
  
  u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
  v = randortho(n)
  
  #generate singular value of matrix A; A^TA has conditional number k
  s_c_diag  = seq(from = 1 / sqrt(k), to = 1, length.out = min(m, n))
  s_c = diag(s_c_diag, nrow=m, ncol=n)  
  
  # sigular value decomposition
  A = u%*%s_c%*%v  # for convexity assumption
  
  #xs = rnorm(n)
  xs  = ones(n, 1)
  #b = A%*%xs + 1 / (1 * 500) * rnorm(m)
  b= A%*%xs 
  # solve(t(A) %*% A, t(A) %*% b)
  # t(A) %*% A 
  return(list(A = A, xs = xs, b = b, n = n, m = m ))
}

Quad_sparse_generator <- function(m = 100, n = 50, k = 30, s = 30){
  ### problem set up    min_x |Ax - b|^2 
  ### here we assume that the true value x is sparse with |x|_0 = s
  ### if m >= n, we generate data matrix A with A^TA has conditional number k, where the largest singualr value of A is 1
  ### if m < n , we generate data matrix A where the largest singular value of A is 1
  
  # input data points, here xs is the true solution that we want to find
  m = m
  n = n
  k = k
  
  
  u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
  v = randortho(n)
  
  #generate singular value of matrix A; A^TA has conditional number k
  s_c_diag  = seq(from = 1 / sqrt(k), to = 1, length.out = min(m, n))
  s_c = diag(s_c_diag, nrow=m, ncol=n)  
  
  # sigular value decomposition
  A = u%*%s_c%*%v  # for convexity assumption
  
  xs  = zeros(n, 1)
  xs[1:s,] = 1
  #b = A%*%xs + 1 / (1 * 5000) * rnorm(m)
  b= A%*%xs
  return(list(A = A, xs = xs, b = b, n = n, m = m ))
  
}

########## Experiments on Convergence Analysis #############
######### general comparision ##########
m = 100 
n = 50
k = 1000
tol = -5
maxIter = 100000

data1 = Quad_generator(m = m, n = n, k = k)
t1 = Sys.time()
RCDM_results = RCDM(data1$A, data1$b, data1$xs, alpha = 1, tol = tol, maxIter = maxIter)
t2 = Sys.time()
RCDM.t = t2 - t1
RCDM_results$k
RCDM.gap = RCDM_results$fx
#plot(RCDM.gap)

t1 = Sys.time()
A_RCDM_results = A_RCDM(data1$A, data1$b, data1$xs, alpha = 1, Sigma = rep(1, data1$n), tol = tol, maxIter =maxIter)
t2 = Sys.time()
A_RCDM.t = t2 - t1 
A_RCDM_results$k
A_RCDM.gap = A_RCDM_results$fx
#plot(A_RCDM.gap)

t1 = Sys.time()
GD_results = GD(data1$A, data1$b, data1$xs, alpha = 0.5, tol = tol, maxIter = maxIter)
t2 = Sys.time()
GD.t = t2 - t1
GD_results$k
GD.gap = GD_results$fx 
#plot(GD.gap)

## ggplot

color <- c("RCDM" = "red", "A_RCDM" = "black", "GD" = "blue")
gg <- ggplot() + 
  geom_line(aes(x = 1:length(RCDM.gap), y = log(abs(RCDM.gap)), color = "RCDM"), size = 0.5)  + 
  geom_line(aes(x = 1:length(A_RCDM.gap),y  = log(abs(A_RCDM.gap)), color = "A_RCDM"), size = 0.5) +
  geom_line(aes(x = 1:length(GD.gap),y  = log(abs(GD.gap)), color = "GD"), size = 0.5) +
  xlab("nums of iteration") + 
  ylab("log of optimality gap") + 
  labs(title = "Optimality gap VS iteration (k = 1000)", color = "Legend") +
  theme(plot.title = element_text(hjust = 0.5, size = 15)) + 
  scale_color_manual(values = color)

show(gg)




########## Experiments on Complexity Analysis #############
########## randomized_CD #############
#strongly convex setting 
data1 = Quad_generator(m = 100, n = 50, k = 30)
RCDM_results = RCDM(data1$A, data1$b, data1$xs, alpha = 1, tol = 0.005)
RCDM_results$k

### gap vs iteration ###
plot(RCDM_results$cr)

### log(1/eps) vs nums of iteration (strongly convex)###
set.seed(100)
m = 100
n = 50
k = 30
data1 = Quad_generator(m = m, n = n, k = k)

N = 500
#eps = seq(0.005, 0.1, length.out = N)
eps = exp(-seq(1, 3, length.out = N))
num_iter = sapply(eps, function(eps){
  # print(eps)
  RCDM_results = RCDM(data1$A, data1$b, data1$xs, alpha = 1, tol = eps)
  RCDM_results$k
})
# plot(log(1 / eps), num_iter)

ggplot() + 
  geom_point(aes(x = log(1 / eps), y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = log(1 / eps), y = num_iter), show.legend = TRUE)  + 
  xlab(expression(log(1 / epsilon))) + 
  ylab("numbers of iteration") +
  labs(title = "Epsilon complexity for RCDM under strong convexity assumption") + 
  theme(plot.title = element_text(hjust = 0.5))



### (1/sigma) vs nums of iteration (strong convex)###
set.seed(100)
N = 500
kappa = seq(1, 30, length.out = N)
num_iter = sapply(kappa, function(kappa){
  m = 100
  n = 50
  k = kappa
  data1 = Quad_generator(m = m, n = n, k = k)
  RCDM_results = RCDM(data1$A, data1$b, data1$xs, alpha = 1, tol = 0.01)
  RCDM_results$k
})

ggplot() + 
  geom_point(aes(x = kappa, y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = kappa, y = num_iter), show.legend = TRUE)  + 
  xlab(expression(frac(1,sigma))) + 
  ylab("numbers of iteration") + 
  labs(title = "Sigma complexity for RCDM under strong convexity assumption") +
  theme(plot.title = element_text(hjust = 0.5))



### (1/eps) vs nums of iteration(convex) ###
set.seed(100)
m = 100
n = 50
k = 100000
data1 = Quad_generator(m = m, n = n, k = k)

N = 500
#eps = seq(0.005, 0.1, length.out = N)
eps = 1 / seq(10, 1000, length.out = N)
num_iter = sapply(eps, function(eps){
  RCDM_results = RCDM(data1$A, data1$b, data1$xs, alpha = 1, tol = eps)
  RCDM_results$k
})
plot(1 / eps, num_iter)

ggplot() + 
  geom_point(aes(x = 1 / eps, y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = 1 / eps, y = num_iter), show.legend = TRUE)  + 
  xlab(expression(frac(1,epsilon))) + 
  ylab("numbers of iteration") + 
  labs(title = "Epsilon complexity for RCDM under convexity assumption") + 
  theme(plot.title = element_text(hjust = 0.5))




########## Seperable_RCD ###########

###  Experiment  ###
# input sparse data points, here xs is the true solution that we want to find
m = 100
n = 50
k = 50
s = 30
data2 = Quad_sparse_generator(m = m, n = n, k = k, s=s)

SpCD_results = SpCD(data2$A, data2$b, data2$xs, lambda = 0.01, alpha = 1, tol = 0.0001)
SpCD_results$k

plot(SpCD_results$fx)

### gap vs iteration ###
plot(SpCD_results$cr)


### eps vs nums of iteration  (strong convex) ###
set.seed(300)
m = 100
n = 50
k = 50
s = 30
data2 = Quad_sparse_generator(m = m, n = n, k = k, s=s)

#produce eps with length N
N = 500
#eps = seq(0.005, 0.1, length.out = N)
eps = exp(-seq(1, 3, length.out = N))
num_iter = sapply(eps, function(eps){
  SpCD_results = SpCD(data2$A, data2$b, data2$xs, lambda = 0.01, alpha = 1, tol = eps)
  SpCD_results$k
})
# plot(log(1 / eps), num_iter)

ggplot() + 
  geom_point(aes(x = log(1 / eps), y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = log(1 / eps), y = num_iter), show.legend = TRUE)  + 
  xlab(expression(log(1 / epsilon))) + 
  ylab("numbers of iteration") + 
  labs(title = "Epsilon complexity for SCDM under strong convexity assumption") + 
  theme(plot.title = element_text(hjust = 0.5))



### sigma vs nums of iteration (strong convex) ### 
set.seed(300)
N = 500
kappa = seq(1, 300, length.out = N)
num_iter = sapply(kappa, function(kappa){
  m = 100
  n = 50
  k = kappa
  s = 10
  tol = 0.001
  data2 = Quad_sparse_generator(m = m, n = n, k = k, s=s)
  SpCD_results = SpCD(data2$A, data2$b, data2$xs, lambda = 0.01, alpha = 1, tol = tol)
  SpCD_results$k
})

ggplot() + 
  geom_point(aes(x = kappa, y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = kappa, y = num_iter), show.legend = TRUE)  + 
  xlab(expression(frac(1,sigma))) + 
  ylab("numbers of iteration") + 
  labs(title = "Sigma complexity for SCDM under strong convexity assumption") +
  theme(plot.title = element_text(hjust = 0.5))


### eps vs nums of iteration  (convex) ###
set.seed(300)
m = 100
n = 200
k = 50
s = 30
data2 = Quad_sparse_generator(m = m, n = n, k = k, s=s)

#produce eps with length N
N = 100
eps = seq(0.005, 0.1, length.out = N)
num_iter = sapply(eps, function(eps){
  SpCD_results = SpCD(data2$A, data2$b, data2$xs, lambda = 0.01, alpha = 1, tol = eps)
  SpCD_results$k
})
# plot(log(1 / eps), num_iter)

ggplot() + 
  geom_point(aes(x = log(1 / eps), y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = log(1 / eps), y = num_iter), show.legend = TRUE)  + 
  xlab(expression(log(frac(1,epsilon)))) + 
  ylab("numbers of iteration") +
  labs(title = "Epsilon complexity for SCDM under convexity assumption") + 
  theme(plot.title = element_text(hjust = 0.5))





######### Accelerated_RCD ##########
##### Experiment #####
# input data points, here xs is the true solution that we want to find
m = 100
n = 50
k = 30

data1 = Quad_generator(m = m, n = n, k = k)

A_RCDM_results = A_RCDM(data1$A, data1$b, data1$xs, alpha = 15, tol = 0.01)
A_RCDM_results = A_RCDM(data1$A, data1$b, data1$xs, alpha = 1, Sigma = rep(1, n), tol = 0.01)
print(paste("The total number of iteration for ARCD algorithm = ", A_RCDM_results$k))

### gap vs iteration ###
plot(A_RCDM_results$cr)


### eps vs nums of iteration (strong convex)###
set.seed(100)
m = 100
n = 50
k = 30
data1 = Quad_generator(m = m, n = n, k = k)

N = 100
# eps = seq(0.005, 0.1, length.out = N)
eps = exp(-seq(1, 3, length.out = N))
num_iter = sapply(eps, function(eps){
  A_RCDM_results = A_RCDM(data1$A, data1$b, data1$xs, alpha = 1, Sigma = rep(1, n), tol = eps)
  A_RCDM_results$k
})
# plot(log(1 / eps), num_iter)

ggplot() + 
  geom_point(aes(x = log(1 / eps), y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = log(1 / eps), y = num_iter), show.legend = TRUE)  + 
  xlab(expression(log(1 / epsilon))) + 
  ylab("numbers of iteration") + 
  labs(title = "Epsilon complexity for ACDM under strong convexity assumption") + 
  theme(plot.title = element_text(hjust = 0.5))



### sigma vs nums of iteration (strong convex) ###
set.seed(100)
N = 100
kappa = seq(1, 30, length.out = N)
num_iter = sapply(kappa, function(kappa){
  m = 100
  n = 50
  k = kappa
  data1 = Quad_generator(m = m, n = n, k = k)
  A_RCDM_results = A_RCDM(data1$A, data1$b, data1$xs, alpha = 1, Sigma = rep(1, n), tol = 0.01)
  A_RCDM_results$k
})

# plot(kappa, num_iter)
ggplot() + 
  geom_point(aes(x = sqrt(kappa)*log(1/kappa), y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = sqrt(kappa)*log(1/kappa), y = num_iter), show.legend = TRUE)  + 
  xlab(expression(sqrt(frac(1, sigma))*log(sigma))) + 
  ylab("numbers of iteration") + 
  labs(title = "Sigma complexity for ACDM under strong convexity assumption") +
  theme(plot.title = element_text(hjust = 0.5))




### eps vs nums of iteration (convex)###
set.seed(100)
m = 100
n = 150
k = 10000
data1 = Quad_generator(m = m, n = n, k = k)

N = 100
# eps = seq(10, 0.1, length.out = N)

eps = 1 / seq(10, 1000, length.out = N)
num_iter = sapply(eps, function(eps){
  A_RCDM_results = A_RCDM(data1$A, data1$b, data1$xs, alpha = 1, Sigma = rep(1, n), tol = eps)
  A_RCDM_results$k
})
# plot(1 / sqrt(eps), num_iter)

ggplot() + 
  geom_point(aes(x = sqrt(1/eps), y = num_iter)) + 
  geom_smooth( method = "lm", aes(x = sqrt(1/eps), y = num_iter), show.legend = TRUE)  + 
  xlab(expression(sqrt(frac(1,epsilon)))) + 
  ylab("numbers of iteration") +
  labs(title = "Epsilon complexity for ACDM under convexity assumption") + 
  theme(plot.title = element_text(hjust = 0.5))

