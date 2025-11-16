library(Rcpp)
library(RcppArmadillo)

sourceCpp("Cwrapper.cpp")
gibbs_adaptive = function(y, wB, nrun, burn, thin, mseed, verbose, p_constant, 
                           b0, b1, start_adapt, alpha, a_sigma, b_sigma, a_theta, 
                          b_theta, sd_gammaB, scale_factor_MH, cMH, kinit = NULL, 
                          kmax = NULL){
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
  kstar = k - 1 # number of active factors
  sp = floor((nrun - burn) / thin) # number of posterior samples
  # adaptive probability
  prob = 1 / exp(b0 + b1 * seq(1, nrun))
  uu = runif(nrun)
  
  #-----------# Initialization #-------------#
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
  Lambda = t(t(Lambda_star) * sqrt(rho)) * sqrt(Phi)
  
  # pivots
  Lcal = c(1:p) # set of all possible values for l
  lpiv = sample(Lcal, k, replace = F) # pivots vector
  Delta = matrix(0, p, k) # sparsity matrix
  for(i in 1:k){
    Delta[c(1:(lpiv[i]-1)),i] = 0
    Delta[lpiv[i],i] = 1
  }
  
  # Allocate output object memory
  output = c("gamma",      # shrinkCoefSamples    : qBxk
             "eta",         # etaval               : nxk
             "lambda",      # loadSamples          : pxk
             "sigmacol",     # sigmacol (1/sigma^2) : p
             "activeFactors"     # rho : k
  )
  
  out = list("numFactors" = NA)
  if("gamma" %in% output) out["gamma"] = NA
  if("eta" %in% output) out["eta"] = NA
  if("lambda" %in% output) out["lambda"] = NA
  if("sigmacol" %in% output) out["sigmacol"] = NA 
  if("activeFactors" %in% output) out["activeFactors"] = NA 
  
  # start time
  t0 = proc.time()
  # -------------------------------------------------------------------------- #
  # ADAPTIVE GIBBS SAMPLING
  # -------------------------------------------------------------------------- #
  out = Rcpp_gibbs(alpha, a_sigma, b_sigma, a_theta, b_theta, sd_gammaB, p_constant,
                   y, wB, burn, nrun, thin, start_adapt, kmax, 
                   eta, GammaB, Lambda,
                   Lambda_star, d, kstar, logit, rho, Phi, Plam, pred, ps, v, w, 
                   out, verbose, uu, prob, sp, lpiv, Delta, scale_factor_MH, cMH)
  # -------------------------------------------------------------------------- #
  if ("sigmacol" %in% output) out[["sigmacol"]] <- lapply(out[["sigmacol"]], c)
  out[["numFactors"]] <- c(out[["numFactors"]])
  out[["time"]] <- (proc.time() - t0)[1]
  out[["y"]] <- y                 # data                       : nxp
  out[["wB"]]  <- wB              # biological meta-covariates : pxqB
  out[["hyperparameters"]] <- list(alpha = alpha, a_theta = a_theta,
                                   b_theta = b_theta, 
                                   sd_gammaB = sd_gammaB, a_sigma = a_sigma, 
                                   b_sigma = b_sigma, p_constant = p_constant)
  
  return(out)
}
