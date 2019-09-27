
rm(list=ls())

library(lars)
library(MASS) # for mvrnorm function

# read source file
source("joint_selection_Github.R")


####################################################
# 2. Data generation
####################################################

n = 100
k = 30 # the number of transcription factors (TFs)
pk = 5 # pk-1 = the number of genes which are regulated by each TF
p = k*pk
TF.ind = ((1:k) - 1)*pk + 1 # indices for TFs


# generation of covariate matrix X (based on AR representation)
A0 = matrix(0, p,p)
for(i in TF.ind){
  A0[i+(1:(pk-1)), i] = 1 
}
D0 = diag(runif(n = p, min = 3, max = 5)) 
X = matrix(0, nrow=n, ncol=p)
X[,1] = rnorm(n, sd = sqrt(D0[1,1]))
for(j in 2:p){
  mean.vec.j = X[, 1:(j-1)]%*%as.matrix(A0[j, 1:(j-1)])
  X[,j] = rnorm(n, mean = mean.vec.j, sd = sqrt(D0[j,j]))
}


# generation of coefficient vector \beta
beta0 = matrix(0, nrow=p, ncol=1)
gam0.loc = 1:(pk*4)
gam0.len = pk*4
beta0[gam0.loc,1] = runif(n = gam0.len, min = 0.5, max = 1)
gamma0 = c(rep(1,gam0.len), rep(0, p-gam0.len))


# Generate a data vector Y : n X 1
sigma2 = sum(beta0^2)/4
Y = X %*% beta0 + rnorm(n = n, mean = rep(0,n), sd = sqrt(sigma2))


####################################################
# hyperparameters
####################################################
a.hp = 2.75 # a in MRF prior
b.hp = 0.5 # b in MRF prior
qn = 0.005 # q for indep. Bernoulli prior

a0 = 0.1 # sigma2 ~ IG(a0, b0)
b0 = 0.01 
c0 = 1 # beta_gamma | gamma, sigma2 ~ N(0, c0*sigma2*I) for some c0 > 0

alpha.plus = 10 # alpha.plus = \alpha_i(D) - \nu_i(D) in DAG-Wishart prior
U = diag(p) # U in DAG-Wishart prior

niter = 5000
nburn = 5000
nadap = 0
init.gam = rep(0, p)

res = DAGWishart(Y, X, qn, alpha.plus, U, init.A=NULL, init.gam=init.gam, Rj=NULL, a.hp, b.hp, tau2, sigma2, a0, b0, c0, niter, nburn, nadap)


post.ind = as.numeric(colMeans(res$gamma.mat)>0.5) # median probability model
post.ind
