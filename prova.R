rm(list=ls())

library(Rcpp)
library(RcppArmadillo)

sourceCpp("Cwrapper.cpp")

# # hyperparameters and adaptation
# alpha = 5 # param CUSP
# a_sigma = 1 
# b_sigma = 1
# a_theta = 1
# b_theta = 1
# sd_gammaB = 1 # Standard deviation for \eqn{\gamma_{hB}} prior distribution.
# p_constant = 0.5
# start_adapt = 20
# 
# # gibbs parameters
# nrun = 100
# burn = round(nrun/4)
# thin = 1
# # adaptation parameters
# b0 = 1
# b1 = 5 * 10^(-4)
# 
# mseed=343
# verbose=T
# y # data
# wB # meta-covariate

