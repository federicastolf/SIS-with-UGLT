library(Rcpp)
library(RcppArmadillo)

sourceCpp("git/Cwrapper.cpp")
#sourceCpp("Cwrapper.cpp")
gibbs_adaptive = function(y, wB, nrun, burn, thin, mseed, verbose, p_constant, 
                           b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta, 
                          b_theta, Sigma_gamma, scale_factor_MH, cMH, y_max = Inf,
                          star = FALSE, kinit = NULL, kmax = NULL, 
                          order_dependent = FALSE, mu_mean0 = 0, mu_sd0 = 10,
                          column_intercept = TRUE){
  set.seed(mseed)
  p = dim(y)[2]
  n = dim(y)[1]
  qB = dim(wB)[2]
  if(is.null(kmax)) {
    kmax = p + 1
  }
  if(is.null(kinit)) {
    kinit = min(floor(log(p)*6), p)
  }
  k = kinit # number of factors to start with (active and inactive)
  kstar = k # number of active factors
  sp = floor((nrun - burn) / thin) # number of posterior samples
  # adaptive probability
  prob = 1 / exp(b0 + b1 * seq(1, nrun))
  uu = runif(nrun)
  
  # Transformation for the response
  a_j = function(j, y_max) {
    val = j
    val[j == y_max + 1] = Inf
    val
  }
  # Bounds for truncated normal
  a_y = a_yp1=- matrix(NA, nrow = n, ncol = p)
  for(j in 1:p) {
    a_y[, j] = a_j(y[, j], y_max)        # a_y = a_j(y)
    a_yp1[, j] = a_j(y[, j] + 1, y_max)  # a_yp1 = a_j(y + 1)
  }
  # Replace NA with 0/Inf in a_y/a_yp1
  a_y[is.na(y)] = 0      # log(0)=-Inf
  a_yp1[is.na(y)] = Inf  # log(Inf)=Inf
  
  #-----------# Initialization #-------------#
  mu = rep(0, p)
  ps = rgamma(p, a_sigma, b_sigma) # sigma^-2
  Lambda_star = matrix(rnorm(p*k), nrow=p, ncol=k) # loading matrix
  eta = matrix(rnorm(n*k), nrow = n, ncol = k) # latent factors
  # Initialize GammaB (qxk) and pred (pxk)
  GammaB = matrix(rnorm(qB * k), nrow = qB, ncol = k) # traits effect on local shrinkage
  pred = wB %*% GammaB # local shrinkage coefficients
  logit = plogis(pred)
  # Initialize Phi pxk
  Phi = matrix(rbinom(p * k, size = 1, prob = p_constant), nrow = p, ncol = k)
  # Initialize pi_h, h = 1, ..., k
  v = c(rbeta(k - 1, shape1 = 1, shape2 = alpha), 1)
  w = v * c(1, cumprod(1 - v[-k]))  # product up to  l - 1
  d = rep(k - 1, k)                 # augmented data
  rho = rep(1, k)                   # preallocation for Bernoulli
  # Initialize the precision matrix of lambda star
  Plam = diag(rgamma(k, a_theta, b_theta))
  # Compute Lambda (pxk)
  
  # pivots
  Lcal = c(1:p) # set of all possible values for l
  lpiv = sample(Lcal, k, replace = F) # pivots vector
  
  Delta = matrix(0, p, k)
  for(i in 1:k){
    l_h = lpiv[i]
    for(j in l_h:p){
      Delta[j, i] = Phi[j, i]
    }
    Delta[l_h, i] = 1  # force pivot to 1
  }

  Lambda = Lambda_star * Delta * matrix(rho, nrow=p, ncol=k, byrow=TRUE)
  
  # Allocate output object memory
  output = c("gamma",      # shrinkCoefSamples    : qBxk
             "eta",         # etaval               : nxk
             "lambda",      # loadSamples          : pxk
             "preccol",     # preccol (1/sigma^2) : p
             "activeFactors"     # rho : k
  )
  if (star) output = c(output, "mu")
  
  out = list("numFactors" = NA)
  if("gamma" %in% output) out["gamma"] = NA
  if("eta" %in% output) out["eta"] = NA
  if("lambda" %in% output) out["lambda"] = NA
  if("preccol" %in% output) out["preccol"] = NA
  if("activeFactors" %in% output) out["activeFactors"] = NA 
  if ("mu" %in% output) out["mu"] = NA
  # start time
  t0 = proc.time()
  # -------------------------------------------------------------------------- #
  # ADAPTIVE GIBBS SAMPLING
  # -------------------------------------------------------------------------- #
  out = Rcpp_gibbs(alpha, a_sigma, b_sigma, a_theta, b_theta, Sigma_gamma, 
                   p_constant, y, wB, burn, nrun, thin, start_adapt, kmax,  eta,
                   GammaB, Lambda, Lambda_star, d, kstar, logit, rho, Phi, Plam,
                   pred, ps, v, w, out, verbose, uu, prob, sp, lpiv, Delta,
                   scale_factor_MH, cMH, a_y, a_yp1, order_dependent, star,
                   mu, mu_mean0, mu_sd0, column_intercept)
  # -------------------------------------------------------------------------- #
  
  if ("preccol" %in% output) out[["preccol"]] <- lapply(out[["preccol"]], c)
  out[["numFactors"]] <- c(out[["numFactors"]])
  out[["time"]] <- (proc.time() - t0)[1]
  out[["y"]] <- y                 # data                       : nxp
  out[["wB"]]  <- wB              # biological meta-covariates : pxqB
  out[["hyperparameters"]] <- list(alpha = alpha, a_theta = a_theta,
                                   b_theta = b_theta, 
                                   Sigma_gamma = Sigma_gamma, a_sigma = a_sigma, 
                                   b_sigma = b_sigma, p_constant = p_constant,
                                   mu_mean0 = mu_mean0, mu_sd0 = mu_sd0,
                                   column_intercept = column_intercept)
  
  return(out)
}

