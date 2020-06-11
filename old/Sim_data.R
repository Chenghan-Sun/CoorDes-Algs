# This .R code file consists of:
  # Simualtion data for Coordinate Descent methods (Experiments)
# Project Group Members:
  # Han Chen, Ninghui Li, Chenghan Sun


library(pracma)

m = 300
n = 500

u = randortho(m)  # Generates random orthonormal or unitary matrix of size m
v = randortho(n)

s_c_diag = runif(min(m, n), min=30, max=50)
s_c_diag[271:300] = 0
s_c = diag(s_c_diag, nrow=m, ncol=n)  # for convexity assumption 
s_sc = diag(runif(min(m, n), min=30, max=50), nrow=m, ncol=n)  # for sigma-strong convexity assumption 

# s -> dim(m X n)
# if (m > n) {
#  s = vbind(s, zeros(m-n, n))
#} else {
#  s = cbind(s, zeros(m, n-m))
#}

# sigular value decomposition
A_c = u%*%s_c%*%v  # for convexity assumption
A_sc = u%*%s_sc%*%v  # for sigma-strong convexity assumption 

# test 
t_scd = svd(A_c)
t_scd$d



