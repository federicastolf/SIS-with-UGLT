library(pheatmap)
library(gridExtra)
library(sparvaride)

rm(list=ls())

source("gibbs.R")
source("helper_functions.R")
source("my_color_palettes.R")
my_palette = my_color_palettes()

#----------------------# scenario 1 : block covariance #------------------------#

# simulate data
n = 50
p = 30
mseed = 435
data_synt <- simulate_block_mvn(n, p, n_blocks = 3, within_cor = c(0.5, 0.7, 0.6), 
                                between_cor = 0, mseed)
Y <- data_synt$data
Sigma <- data_synt$Sigma
pheatmap(data_synt$Sigma, cluster_rows = F, cluster_cols = F)
# library(corrplot)
# corrplot(result$Sigma, method = "color", tl.cex = 0.5)

xf = as.factor(c(rep(1,p/3), rep(2,p/3), rep(3,p/3)))
covariate = model.matrix(~xf) # X matrix
covariate1 = matrix(1, p, 3)

# gibbs parameters
nrun = 5000
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
sd_gammaB = 1 # Standard deviation for \eqn{\gamma_{hB}} prior distribution.
p_constant = log(p)/p*2*exp(1)
scale_factor_MH = 1
cMH = 0.001

fit = gibbs_adaptive(Y, covariate, nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH, cMH)
fit_nocov = gibbs_adaptive(Y, covariate1, nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH, cMH)

# ---- Covariance ----

Lambda_outer <- lapply(fit$lambda, function(A) A %*% t(A))
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

p1 = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(cov_mean)),-0.0001,length=50),0,
                         seq(0.0001,max(abs(cov_mean)),length=50)))

Lambda_outer_nocov <- lapply(fit_nocov$lambda, function(A) A %*% t(A))
cov_mean_nocov = apply(simplify2array(Lambda_outer_nocov), c(1, 2), mean)

p3 = pheatmap(cov_mean_nocov, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(cov_mean_nocov)),-0.0001,length=50),0,
                         seq(0.0001,max(abs(cov_mean_nocov)),length=50)))

grid.arrange(p1[[4]], p3[[4]],  nrow = 1)


#---- Inference on loadings ----

ppLambda <- identify_lambda(fit, p)
Lambda_mean = apply(simplify2array(ppLambda$lambda), c(1, 2), mean)

pp1 = pheatmap(Lambda_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(Lambda_mean)),-0.0001,length=50),0,
                         seq(0.0001,max(abs(Lambda_mean)),length=50)))

ppLambda_nocov <- identify_lambda(fit_nocov, p)
Lambda_mean_nocov = apply(simplify2array(ppLambda_nocov$lambda), c(1, 2), mean)

pp3 = pheatmap(Lambda_mean_nocov, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
               cluster_cols = F, border_color ="NA", legend=T,
               color = my_palette$sparse(100), 
               breaks = c(seq(-max(abs(Lambda_mean_nocov)),-0.0001,length=50),0,
                          seq(0.0001,max(abs(Lambda_mean_nocov)),length=50)))

grid.arrange(pp1[[4]],  pp3[[4]], nrow = 1)


#----------------------# scenario 2 : linear dependence #------------------------#

data_synt2 <- simulate_decreasing_cov(p, max_var = 15, min_var = 1, max_cor = 0.5, 
                                      min_cor = 0, mseed = mseed)
Y2 <- data_synt2$data
pheatmap(data_synt2$Sigma, cluster_rows = F, cluster_cols = F)

covariate_lin = cbind(rep(1,p), data_synt2$x_norm)
covariate_lin2 = matrix(1, p, 2)

fitLin = gibbs_adaptive(Y2, covariate_lin, nrun, burn, thin, mseed, verbose = T, 
                        p_constant, b0, b1, start_adapt, alpha, a_sigma, b_sigma, 
                        a_theta, b_theta, sd_gammaB, scale_factor_MH, cMH)
fitLin_nocov = gibbs_adaptive(Y2, covariate_lin2, nrun, burn, thin, mseed, verbose = T, p_constant,
                           b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                           b_theta, sd_gammaB, scale_factor_MH, cMH)

# ---- Covariance ----

Lambda_outer_lin <- lapply(fitLin$lambda, function(A) A %*% t(A))
cov_mean_lin = apply(simplify2array(Lambda_outer_lin), c(1, 2), mean)

pl1 = pheatmap(cov_mean_lin, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(cov_mean_lin)),-0.0001,length=50),0,
                         seq(0.0001,max(abs(cov_mean_lin)),length=50)))

Lambda_outer_linnocov <- lapply(fitLin_nocov$lambda, function(A) A %*% t(A))
cov_meanlin_nocov = apply(simplify2array(Lambda_outer_linnocov), c(1, 2), mean)

pl3 = pheatmap(cov_meanlin_nocov, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(cov_meanlin_nocov)),-0.0001,length=50),
                         0,seq(0.0001,max(abs(cov_meanlin_nocov)),length=50)))

grid.arrange(pl1[[4]], pl3[[4]],  nrow = 1)

#---- Inference on loadings ----

ppLambda_lin <- identify_lambda(fitLin, p)
Lambda_mean_lin = apply(simplify2array(ppLambda_lin$lambda), c(1, 2), mean)

ppl1 = pheatmap(Lambda_mean_lin, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
               cluster_cols = F, border_color ="NA", legend=T,
               color = my_palette$sparse(100), 
               breaks = c(seq(-max(abs(Lambda_mean_lin)),-0.0001,length=50),0,
                          seq(0.0001,max(abs(Lambda_mean_lin)),length=50)))

ppLambda_nocov_lin <- identify_lambda(fitLin_nocov, p)
Lambda_mean_nocov_lin = apply(simplify2array(ppLambda_nocov_lin$lambda), c(1, 2), mean)

ppl3 = pheatmap(Lambda_mean_nocov_lin, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
               cluster_cols = F, border_color ="NA", legend=T,
               color = my_palette$sparse(100), 
               breaks = c(seq(-max(abs(Lambda_mean_nocov_lin)),-0.0001,length=50),0,
                          seq(0.0001,max(abs(Lambda_mean_nocov_lin)),length=50)))

grid.arrange(ppl1[[4]],  ppl3[[4]], nrow = 1)

