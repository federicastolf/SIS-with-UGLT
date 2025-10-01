rm(list=ls())

source("gibbs.R")

n = 30
p = 10
y = matrix(rnorm(n*p),n,p)
xf = as.factor(c(rep(1,5), rep(2,5)))
wB = model.matrix(~xf) # meta-covariate matrix
# gibbs parameters
nrun = 50
burn = round(nrun/4)
thin = 1
mseed = 434
# adaptation parameters
b0 = 1
b1 = 5 * 10^(-4)
start_adapt = 5
# hyperparameters 
alpha = 5 # param CUSP
a_sigma = 1
b_sigma = 1
a_theta = 1
b_theta = 1
sd_gammaB = 1 # Standard deviation for \eqn{\gamma_{hB}} prior distribution.
p_constant = 0.5
scale_factor_MH = 1

fit1 = gibbs_adaptive(y, wB, nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH)








