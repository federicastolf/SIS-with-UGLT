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
data_synt = simulate_block_mvn(n, p, n_blocks = 3, within_cor = c(0.4, 0.9, 0.6),
                               between_cor = 0,variances = 1, mseed)

Y = data_synt$data
covariate = model.matrix(~data_synt$xlab) # X matrix

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
p_constant = 0.05
scale_factor_MH = 1
cMH = 0.001

fit = gibbs_adaptive(Y, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                     b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                     b_theta, sd_gammaB, scale_factor_MH, cMH)

# ---- Covariance ----
Lambda_outer = lapply(fit$lambda, function(A) A %*% t(A))
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

max_val = max(data_synt$Sigma)
p3 = pheatmap(data_synt$Sigma, cluster_rows = F, cluster_cols = F, 
              border_color ="NA",   main = "True covariance",
              breaks = seq(0, max_val, length.out = 100),
              color = colorRampPalette(c("white","yellow2", "orange","darkred"))(100))

p1 = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              main = "Posterior covariance",
              breaks = seq(0, max_val, length.out = 100),
              color = colorRampPalette(c("white","yellow2", "orange","darkred"))(100))

pblock3 = grid.arrange(p3[[4]], p1[[4]], nrow = 1)
# ggsave("combined_heatmaps.pdf", pblock3, width = 12, height = 6)

#---- Inference on loadings ----

ppLambda = identify_lambda(fit, p)
Lambda_mean = apply(simplify2array(ppLambda$lambda), c(1, 2), mean)

pp1 = pheatmap(Lambda_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
               cluster_cols = F, border_color ="NA", legend=T,
               main="Lambda",
               breaks = seq(min(Lambda_mean), max(Lambda_mean), length.out = 100),
               color = colorRampPalette(c("white","gold2", "orange","darkred"))(100))

#-------------# reorder blocks #---------#

data_reordered = reorder_blocks(Sigma = data_synt$Sigma, Y = Y,
                                block_membership = data_synt$xlab,
                                new_block_order = c(3, 1, 2))
Yr = data_reordered$Y
covariater = model.matrix(~data_reordered$xlab) # X matrix
fit_r = gibbs_adaptive(Yr, covariater, nrun, burn, thin, mseed, verbose = T, p_constant,
                       b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                       b_theta, sd_gammaB, scale_factor_MH, cMH)

#--# Covariance 
Lambda_outerr = lapply(fit_r$lambda, function(A) A %*% t(A))
cov_meanr = apply(simplify2array(Lambda_outerr), c(1, 2), mean)


pr3 = pheatmap(data_reordered$Sigma, cluster_rows = F, cluster_cols = F, 
               border_color ="NA",   main = "True covariance- reordered",
               breaks = seq(0, max_val, length.out = 100),
               color = colorRampPalette(c("white","yellow2", "orange","darkred"))(100))
pr1 = pheatmap(cov_meanr, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
               cluster_cols = F, border_color ="NA", legend=T,
               main = "Posterior covariance-reordered",
               breaks = seq(0, max_val, length.out = 100),
               color = colorRampPalette(c("white","yellow2", "orange","darkred"))(100))

prblock3 = grid.arrange(pr3[[4]], pr1[[4]], nrow = 1)
# ggsave("reordere_heatmaps.pdf", prblock3, width = 12, height = 6)


#----# Lambda
ppLambdar = identify_lambda(fit_r, p)
Lambda_meanr = apply(simplify2array(ppLambdar$lambda), c(1, 2), mean)

ppr1 = pheatmap(Lambda_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
                cluster_cols = F, border_color ="NA", legend=T,
                main="Lambda",
                breaks = seq(min(Lambda_mean), max(Lambda_mean), length.out = 100),
                color = colorRampPalette(c("white","gold2", "orange","darkred"))(100))