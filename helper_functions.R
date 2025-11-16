# function to simulate synthetic data
simulate_data = function(n, p, k, covariate, sd_gamma, cp, sigma_sq, mseed){
  set.seed(mseed)
  q = dim(covariate)[2]
  gamma = matrix(rnorm(q*p, sd = sd_gamma), q,p) # regression coefficient
  pilogit = cp*plogis(covariate %*% gamma) 
  get_phi = simPhi(p, k, pilogit)
  Phi = get_phi$Phi
  theta_inv = rgamma(k,1,1)
  Lambda = matrix(0,p,k)
  for(j in 1:p){
    for(h in 1:k){
      Lambda[j,h]= Phi[j,h]/theta_inv[h]
    }
  }
  M = matrix(rnorm(n*k, 0, 1), ncol=k) # factor scores
  # data
  Y = M %*% t(Lambda) + sqrt(sigma_sq)*matrix(rnorm(n*p), nrow = n, ncol = p) 
  return(list("Y" = Y, "Lambda" = Lambda, "Phi"=Phi, "pivot"=get_phi$pivot,
              "delta"=get_phi$Delta, "gamma"=gamma))
}

# function for simulate Phi via gibbs-sampler
simPhi = function(p, k, pilogit, NiterGibbs = 200){
  # Initialization 
  phi = matrix(rbinom(p*k, 1, 0.2), p, k) 
  delta = phi
  l=rep(NA, k) # pivots
  l[1] = which(phi[,1]>0)[1]
  for(h in 2:k){
    w = which(phi[,h]>0) # counting 1's in the h column of Delta
    l[h] = w[!w %in% l[-h]][1] # take the first
    if(is.na(l[h])){
      delta[,h] = 0 # put all column to zeros
    }else{
      delta[1:l[h],h] = 0 # before pivots all zeros
      delta[l[h]:p, h] =  phi[l[h]:p,h] # from pivot and after equal to phi[,l]
    }
  }
  
  for(t in 2:NiterGibbs){
    for(h in 1:k){ 
      Lh = sort(l[-h])
      p_h = pilogit[,h]
      k = 0     
      # definisco prob di phi_h condizionatamente agli altri pivot
      for(pl in Lh){
        k = k+1
        if (pl-k>0){
          scale = 1-prod(1-pilogit[1:pl,h][-Lh[1:k]])
          p_h[pl] = min(p_h[pl]/scale, 1)
        }
        else{
          p_h[pl] = 1
        }
      }
      phi[,h] = rbinom(p,1,p_h) # generate phi
      w = which(phi[,h]>0)
      l[h] = w[!w %in% l[-h]][1]
      if(is.na(l[h])){
        delta[,h] = 0
      }else{
        delta[1:l[h],h] = 0
        delta[l[h]:p, h] =  phi[l[h]:p,h]
      }
    }
  }
  # return the last phi and delta
  return(list("Phi" = phi, "Delta"= delta, "pivot"=l))
}


# function to reorder and sign swap based on delta and lambda
reorder_and_sign_swap <- function(delta, lambda) {
  # Ensure both have same dimensions
  if (!all(dim(delta) == dim(lambda))) {
    stop("delta and lambda must have the same dimensions")
  }
  
  p <- nrow(lambda)
  k <- ncol(lambda)
  
  # Identify first nonzero position and sign for each column
  pivot_pos <- rep(NA_integer_, k)
  pivot_sign <- rep(1, k)
  
  for (h in seq_len(k)) {
    nz <- which(lambda[, h] != 0)
    if (length(nz) > 0) {
      pivot_pos[h] <- nz[1]
      pivot_sign[h] <- sign(lambda[nz[1], h])
    } else {
      pivot_pos[h] <- p + 1  # if column is all zeros, put it last
    }
  }
  
  # Order columns by first nonzero position
  order_idx <- order(pivot_pos)
  pivot_pos_sorted = sort(pivot_pos)
  
  # Apply the order
  delta_new <- delta[, order_idx, drop = FALSE]
  lambda_new <- lambda[, order_idx, drop = FALSE]
  pivot_sign <- pivot_sign[order_idx]
  
  # Flip sign where needed
  for (h in seq_len(k)) {
    if (pivot_sign[h] < 0) {
      lambda_new[, h] <- -lambda_new[, h]
    }
  }
  
  list(delta = delta_new, lambda = lambda_new, order = order_idx,
       pivot_location = pivot_pos_sorted)
}

# find modal vector
find_modal_vector <- function(pivot_list) {
  # convert each vector to a string key
  keys <- vapply(pivot_list, function(v) paste(v, collapse = ","), character(1))
  
  # count frequency
  freq <- table(keys)
  
  # find the most frequent key
  modal_key <- names(freq)[which.max(freq)]
  
  # convert back to numeric vector
  if (nzchar(modal_key)) {
    modal_vec <- as.integer(strsplit(modal_key, ",", fixed = TRUE)[[1]])
  } else {
    modal_vec <- integer(0)  # handle empty vector case
  }
  
  list(mode = modal_vec, count = max(freq), freq_table = freq)
}
