library(pheatmap)
library(gridExtra)
library(sparvaride)
library(tidyverse)

rm(list=ls())

source("gibbs.R")
source("helper_functions.R")

#----------------------# 3 block covariance structures #-----------------------#

# simulate data
n = 50
p = 30
mseed = 435
out <- simulate_factor_cov_block(p, K = 10,  cross_block_sd = 0.1, 
                                 n_blocks = 3, psi = 5, seed = mseed, 
                                 loading_sd = 2)
Sigma <- out$Sigma
Lambda <- out$Lambda
pheatmap(Lambda, cluster_rows = F, cluster_cols = F)
pheatmap(Sigma, cluster_rows = F, cluster_cols = F)

Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
covariate = model.matrix(~out$xlab-1) # 

# gibbs parameters
nrun = 8000
burn = round(nrun/2)
thin = 1
# adaptation parameters
b0 = 1
b1 = 5 * 10^(-4)
start_adapt = 100
# hyperparameters 
alpha = 5 # param CUSP
a_sigma = 1
b_sigma = 1
a_theta = 1
b_theta = 1
sd_gammaB = 1 
p_constant = 0.2
scale_factor_MH = 1
cMH = 0.001

fit = gibbs_adaptive(Y, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = F)

# ---- Covariance ----
Lambda_outer = lapply(fit$lambda, function(A) A %*% t(A))
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

p3 = pheatmap(cov2cor(Sigma), cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance")

p1 = pheatmap(cov2cor(cov_mean), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Posterior covariance")

pblock3 = grid.arrange(p3[[4]], p1[[4]], nrow = 1)


################################################################################
#-----------------------# different data generating #--------------------------#

out <- simulate_factor_block_with_cross(p = p, block_sizes = c(10,10,10),
                                        block_rho = c(0.4, 0.9, 0.6),
                                        cross_rho = 0.3, loading_noise_sd = 0.1, 
                                        seed = mseed)

Sigma <- out$Sigma
Lambda <- out$Lambda
pheatmap(Lambda, cluster_rows = F, cluster_cols = F)
pheatmap(Sigma, cluster_rows = F, cluster_cols = F)

Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
covariate = model.matrix(~out$xlab-1) # 

fit = gibbs_adaptive(Y, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = F)

# ---- Covariance ----
Lambda_outer = lapply(fit$lambda, function(A) A %*% t(A))
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

max_val = max(Sigma)
p3 = pheatmap(cov2cor(Sigma), cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance")

p1 = pheatmap(cov2cor(cov_mean), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Posterior covariance")

pblock3 = grid.arrange(p3[[4]], p1[[4]], nrow = 1)
