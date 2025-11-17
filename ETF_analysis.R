# ---- library ----

library(readxl)
library(pheatmap)
library(sparvaride)

rm(list = ls())

# colors
source("my_color_palettes.R")
my_palette = my_color_palettes()


# ---- Data loading ----

library(readxl)

datafilepath = "../Data/ETF_data.xlsx"
#datafilepath = "ETF_data.xlsx"

# relative returns
data = read_excel(datafilepath, sheet = "ETF RR")
data = as.data.frame(data)

# rownames
date = as.data.frame(read_excel(datafilepath, sheet = "ETF price",
                  skip = 1,
                  col_types = c("date", rep("numeric",32))))[,1]
rownames(data) = date[2:753]

# meta covariate
Info_etf =  as.data.frame(read_excel(datafilepath, sheet = "ETFs"))

x = model.matrix(~Info_etf[,1]-1) 
colnames(x) = sort(unique(Info_etf[,1]))



# ---- Exploratory Analysis ----
p = ncol(data)

# Exploratory plots
matplot(data, type = "l", ylim = c(-0.2,0.2), col = my_palette$scatter(p))

# remove missing values and outliers
apply(data, 1, function(x) sum(abs(x)>0.5))
# all in the last 4 days

data = data[1:749,]
n = nrow(data)

## put histograms on the diagonal
panel.hist <- function(x, ...){
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "gray", ...)
}

pairs(data[,x[,1]==1], diag.panel = panel.hist)
pairs(data[,x[,2]==1], diag.panel = panel.hist)
pairs(data[,x[,3]==1], diag.panel = panel.hist)
pairs(data[,x[,4]==1], diag.panel = panel.hist)
pairs(data[,x[,5]==1], diag.panel = panel.hist)
# strong within sector positive correlation
# gaussian is fine

library(corrplot)
corrplot(cor(data), col = my_palette$sparse(100))
# highly correlated, maybe too much?
# commodities and crypto uncorrelated

# let's scale the data to set more easily hyperpar
y = scale(data)

# ---- Function loading ----
source("gibbs.R")
source("helper_functions.R")


# ---- Factor model estimation ----
# Gibbs parameters
nrun = 12000
burn = round(nrun/4)
thin = 1
mseed = 123
# adaptation parameters
b0 = 1
b1 = 5 * 10^(-4)
start_adapt = 100
# hyperparameters 
alpha = 15 # param CUSP
a_sigma = 3
b_sigma = 2
a_theta = 3
b_theta = 2
sd_gammaB = 1 # Standard deviation for \eqn{\gamma_{hB}} prior distribution.
p_constant = 0.01#log(p)/p*2*exp(1)
scale_factor_MH = 1
cMH = 0.0001

fit1 = gibbs_adaptive(y, x, nrun, burn, thin, mseed, verbose = T, p_constant,
                      b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta,
                      b_theta, sd_gammaB, scale_factor_MH, cMH)

# ---- checks ----
fit1$delta[[1000]]

nsample = length(fit1$numFactors)
plot(fit1$numFactors, type = "l",  col=my_palette$scatter(5)[5])
summary(fit1$numFactors, type = "l")

plot(fit1$accRat, type = "l",  col=my_palette$scatter(5)[5])
summary(fit1$accRate) 
vec = sapply(fit1$gamma, function(m) m[1, 1])
plot(vec, type="l", col=my_palette$scatter(5)[1], ylim= c(-5,5))
lines(1:nsample, sapply(fit1$gamma, function(m) m[2, 1]),
      col=my_palette$scatter(5)[2])
lines(1:nsample, sapply(fit1$gamma, function(m) m[3, 1]),
      col=my_palette$scatter(5)[3])
lines(1:nsample, sapply(fit1$gamma, function(m) m[4, 1]),
      col=my_palette$scatter(5)[5])
colnames(x)
# there was a switch, but in general the x seems relevant


# ---- Covariance ----

matplot(t(simplify2array(fit1$sigmacol)), type = "l",  col=my_palette$scatter(18))
Sigma_mean <- diag(rowMeans(simplify2array(fit1$sigmacol)))
diag(Sigma_mean)
# sam residual variance very high: idunno

Lambda_outer <- lapply(fit1$lambda, function(A) A %*% t(A))
cov_mean = apply(simplify2array(Lambda_outer), c(1, 2), mean)

p1 = pheatmap(cov_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(cov_mean)),-0.0001,length=50),0,seq(0.0001,max(abs(cov_mean)),length=50))
)
p1 = pheatmap(cov2cor(cov_mean+ Sigma_mean), treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-1,-0.0001,length=50),0,seq(0.0001,1,length=50))
              #breaks = c(seq(-max(abs(cov_mean)),-0.0001,length=50),0,seq(0.0001,max(abs(cov_mean)),length=50))
)


#---- Loadings Post-processing ----


# rule 3579 - library sparvaride, ordering of pivots with positive sign

lambda_ordered = list()
delta_ordered = list()
pivot_location = list()
c =0

for(i in 1:nsample){
  delta = matrix(fit1$delta[[i]][,c(fit1$activeFactors[[i]]==T)],
                 nrow = p)
  lambda = matrix(fit1$lambda[[i]][,c(fit1$activeFactors[[i]]==T)],
                  nrow = p)
  
  # if zero columns are present, not identifiable
  if((sum(colSums(delta)==0)==0)){
    
    # if zero rows are present, we check counting rule only on non zero rows
    delta_tocheck = as.matrix(delta[rowSums(delta)>0,])
    
    # check 3579 rule
    if(counting_rule_holds(delta_tocheck)){
      c =c+1
      
      # order pivots decreasing and sign switch when needed
      neword = reorder_and_sign_swap(delta, lambda)
      lambda_ordered[[c]] = neword$lambda
      delta_ordered[[c]] = neword$delta
      pivot_location[[c]] = neword$pivot_location
    }
  }
  
}
c

# mode of pivot locations
modal_config = find_modal_vector(pivot_location)
modal_config$mode
modal_rank = length(modal_config$mode)

c=0
lambda_identified = delta_identified =list()
# keep only modal lambdas
for(i in 1:length(lambda_ordered)){
  if(length(pivot_location[[i]])==modal_rank){
    if(all(pivot_location[[i]]==modal_config$mode)){
      c = c+1
      lambda_identified[[c]] = lambda_ordered[[i]] 
      delta_identified[[c]] = delta_ordered[[i]]
    }
  }
}
c


#---- Inference on loadings ----

Lambda_mean = apply(simplify2array(lambda_identified), c(1, 2), mean)
rownames(Lambda_mean)=c(1:p)

p1 = pheatmap(Lambda_mean, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
              cluster_cols = F, border_color ="NA", legend=T,
              color = my_palette$sparse(100), 
              breaks = c(seq(-max(abs(Lambda_mean)),-0.0001,length=50),0,
                         seq(0.0001,max(abs(Lambda_mean)),length=50))
)
pheatmap(x, treeheight_row = 0, treeheight_col = 0, cluster_rows = F,
         cluster_cols = F, border_color ="NA", legend=T,
         color = my_palette$sparse(100), 
         breaks = c(seq(-max(abs(x)),-0.0001,length=50),0,
                    seq(0.0001,max(abs(x)),length=50))
)

