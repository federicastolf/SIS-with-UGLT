library(pheatmap)

rm(list=ls())

source("gibbs.R")
source("helper_functions.R")

# simulate data
n = 50
p = 30
k = 10
sd_gamma0 = 5
cp0 = 0.5
sigma_sq = 0.1
mseed = 434
xf = as.factor(c(rep(1,p/3), rep(2,p/3), rep(3,p/3)))
cov = model.matrix(~xf) # X matrix
data_synt = simulate_data(n, p, k, cov, sd_gamma0, cp0, sigma_sq, mseed)

y = data_synt$Y

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
p_constant = 0.2
# p_constant <-  exp(2)* log(p) / p
scale_factor_MH = 1
cMH = 0.7

fit1 = gibbs_adaptive(y, cov, nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH, cMH)

summary(fit1$accRate) 
vec = sapply(fit1$gamma, function(m) m[1, 2])
plot(vec, type="l")

Lambda_outer <- lapply(fit1$lambda, function(A) A %*% t(A))
covariance <- lapply(Lambda_outer, function(x) (x != 0) * 1)
cov_mean = apply(simplify2array(covariance), c(1, 2), mean)

p1 = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=F,
              color = colorRampPalette(c("white", "blue"))(60),
              breaks = seq(min(cov_mean),max(cov_mean), length.out = 60))
p1
