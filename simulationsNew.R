library(pheatmap)
library(gridExtra)
library(sparvaride)
library(tidyverse)

rm(list=ls())

source("gibbs.R")
source("helper_functions.R")

#----------------------# 3 block covariance structures #-----------------------#
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
p_constant = 0.8
scale_factor_MH = 1
cMH = 0.001

# simulate data
n = 50
p = 30
mseed = 435
out <- simulate_factor_block_with_cross(p = p, block_sizes = c(10,10,10),
                                        block_rho = c(0.7, 0.5, 0.6),
                                        cross_rho = 0.1, loading_noise_sd = 0.1, 
                                        seed = mseed, diag_var = 2)
Sigma <- out$Sigma
Lambda <- out$Lambda
pheatmap(Lambda, cluster_rows = F, cluster_cols = F)
pheatmap(Sigma, cluster_rows = F, cluster_cols = F)

Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
corrplot::corrplot(cov2cor(cor(Y)))
covariate = model.matrix(~out$xlab-1) # 

#--- model non order-dependent
fit = gibbs_adaptive(Y, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = F)

Lambda_outer <- mapply(function(A, s) {A %*% t(A) + diag(1/s)},
                       fit$lambda, fit$sigmacol, SIMPLIFY = FALSE)
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

p3 = pheatmap(cov2cor(Sigma), cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance", 
              color = colorRampPalette(c("white","grey50", "grey10"))(100))
p1 = pheatmap(cov2cor(cov_mean), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "No order",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))


#--- model order-dependent
fit2 = gibbs_adaptive(Y, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = T)

Lambda_outer2 <- mapply(function(A, s) {A %*% t(A) + diag(1/s)},
                       fit2$lambda, fit2$sigmacol, SIMPLIFY = FALSE)
cov_mean2 = apply(simplify2array(Lambda_outer2), c(1, 2), mean)
p2 = pheatmap(cov2cor(cov_mean2), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Order dependent",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))
# correlation matrix
pblock3 = grid.arrange(p3[[4]], p1[[4]], p2[[4]], nrow = 1)


# covariance matrix
p3 = pheatmap(Sigma, cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))
p1 = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "No order",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))
p2 = pheatmap(cov_mean2, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Order dependent",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))       
pblock3 = grid.arrange(p3[[4]], p1[[4]], p2[[4]], nrow = 1)



p1_bis = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
                  cluster_cols = F, border_color ="NA", legend=T,
                  main = "Posterior covariance",
                  color = colorRampPalette(c("white","grey50", "grey10"))(100))

#----------------------# Reordered #-----------------------#

new_order = c(21:30,1:10,11:20)
Y2 = Y[,new_order]

fit = gibbs_adaptive(Y2, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = F)

Lambda_outer <- mapply(function(A, s) {A %*% t(A) + diag(1/s)},
                       fit$lambda, fit$sigmacol, SIMPLIFY = FALSE)
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)


fit2 = gibbs_adaptive(Y2, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = T)

Lambda_outer2 <- mapply(function(A, s) {A %*% t(A) + diag(1/s)},
                        fit2$lambda, fit2$sigmacol, SIMPLIFY = FALSE)
cov_mean2 = apply(simplify2array(Lambda_outer2), c(1, 2), mean)



# covariance matrix
p6 = pheatmap(Sigma[new_order,new_order], cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))
p4 = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "No order",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))
p5 = pheatmap(cov_mean2, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Order dependent",
              color = colorRampPalette(c("white","grey50", "grey10"))(100))       
pblock3_n = grid.arrange(p3[[4]], p1[[4]], p2[[4]],
  p6[[4]], p4[[4]], p5[[4]], nrow = 2)


p4_bis = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
                  cluster_cols = F, border_color ="NA", legend=T,
                  main = "Posterior covariance",
                  color = colorRampPalette(c("white","grey50", "grey10"))(100))

grid.arrange(p3[[4]], p1_bis[[4]], p6[[4]], p4_bis[[4]], nrow = 1)




#----------------------# Spatio-temporal covariance structures #-----------------------#

library(geoR)

?varcov.spatial


p = 36
n = 50
mseed = 435

# monthly seasonality
month_cov = simulate_factor_block_with_cross(p = p, block_sizes = rep(3,12),
                                 block_rho = runif(12, 0.5,0.9),
                                 cross_rho = 0.1, loading_noise_sd = 0.1, 
                                 seed = mseed, diag_var = 2)$Sigma

month_cov = month_cov[c(seq(1,36,3),seq(2,36,3),seq(3,36,3) ),
          c(seq(1,36,3),seq(2,36,3),seq(3,36,3) )]

# timepoint
t = seq(1, p)

# matern covariance structure
temp_cov = varcov.spatial(coords = cbind(t,t), cov.model = "matern",
                     cov.pars = c(2, 10))$varcov
pheatmap(temp_cov, cluster_rows = F, cluster_cols = F, 
         border_color ="NA",   main = "True covariance")

Sigma = month_cov + temp_cov

Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
corrplot::corrplot(cov2cor(cor(Y)))


covariate = model.matrix( ~ ., 
                          data = data.frame(t, as.factor(rep(1:3, each = 12)),
                         as.factor(rep(1:12, 3))) )


#--- model non order-dependent
p_constant = 0.8
#b_theta = 5
#a_sigma = 5
alpha = 10
fit = gibbs_adaptive(Y, cbind(1,scale(t), ), #as.matrix(rep(1,p)),
                     nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = F)

Lambda_outer <- mapply(function(A, s) {A %*% t(A) + diag(1/s)},
                       fit$lambda, fit$sigmacol, SIMPLIFY = FALSE)
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

p3 = pheatmap(cov2cor(Sigma), cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance", colorRampPalette(c("white","grey", "grey10"))(100))
p1 = pheatmap(cov2cor(cov_mean), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "No order", colorRampPalette(c("white","grey", "grey10"))(100))


#--- model order-dependent
fit2 = gibbs_adaptive(Y, cbind(1,scale(t)), nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH, cMH, order_dependent = T)

Lambda_outer2 <- mapply(function(A, s) {A %*% t(A) + diag(1/s)},
                        fit2$lambda, fit2$sigmacol, SIMPLIFY = FALSE)
cov_mean2 = apply(simplify2array(Lambda_outer2), c(1, 2), mean)
p2 = pheatmap(cov2cor(cov_mean2), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Order dependent", colorRampPalette(c("white","grey", "grey10"))(100))
# correlation matrix
pblock3 = grid.arrange(p3[[4]], p1[[4]], p2[[4]], nrow = 1)
