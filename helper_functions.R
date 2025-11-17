library(MASS)

# simulate data with bloc-structure covariance
simulate_block_mvn <- function(n, p, n_blocks, within_cor, between_cor, mseed) {
  set.seed(mseed)
  block_size <- p / n_blocks
  # Allow different correlations for each block
  if (length(within_cor) == 1) {
    within_cor <- rep(within_cor, n_blocks)
  } else if (length(within_cor) != n_blocks) {
    stop("within_cor must be length 1 or n_blocks")
  }
  
  # Initialize covariance matrix
  Sigma <- matrix(between_cor, p, p)
  
  # Fill in blocks
  for (i in 1:n_blocks) {
    start_idx <- (i - 1) * block_size + 1
    end_idx <- i * block_size
    # Set within-block correlation
    Sigma[start_idx:end_idx, start_idx:end_idx] <- within_cor[i]
    # Set diagonal to 1
    diag(Sigma[start_idx:end_idx, start_idx:end_idx]) <- 1
  }
  
  # Simulate data
  Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  return(list(data = Y, Sigma = Sigma))
}

# Covariance decreases with index
simulate_decreasing_cov <- function(p, x = NULL, max_var, min_var, max_cor, 
                                    min_cor, x_range = c(0, 10), mseed) {
  set.seed(mseed)
  # Generate covariate if not provided
  if (is.null(x)) {
    x <- seq(x_range[1], x_range[2], length.out = p)
  } else {
    if (length(x) != p) {
      stop("Length of x must equal p")
    }
  }
  # Normalize x to [0, 1]
  x_min <- min(x)
  x_max <- max(x)
  x_norm <- (x - x_min) / (x_max - x_min)
  
  # Variances DECREASE linearly with x
  # At x_min: max_var, at x_max: min_var
  variances <- max_var - (max_var - min_var) * x_norm
  sds <- sqrt(variances)
  
  # Correlation matrix: correlation DECREASES with x
  R <- matrix(0, p, p)
  
  for (i in 1:p) {
    for (j in 1:p) {
      if (i == j) {
        R[i, j] <- 1
      } else {
        # Correlation based on average of normalized x values
        avg_x_norm <- (x_norm[i] + x_norm[j]) / 2
        R[i, j] <- max_cor - (max_cor - min_cor) * avg_x_norm
      }
    }
  }
  
  # Ensure positive definiteness
  eigenvalues <- eigen(R, only.values = TRUE)$values
  if (any(eigenvalues <= 0)) {
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE)$mat)
    warning("Matrix adjusted for positive definiteness")
  }
  
  # Create covariance matrix: Sigma = D * R * D
  D <- diag(sds)
  Sigma <- D %*% R %*% D
  
  Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  return(list(data = Y, Sigma = Sigma, x = x, x_norm = x_norm))
}



identify_lambda <- function(fit, p, nsample = NULL, 
                            return_delta = TRUE,
                            return_pivot_info = FALSE) {
  if (is.null(nsample)) {
    nsample <- length(fit$lambda)
  }
  
  # Initialize lists
  lambda_ordered <- list()
  delta_ordered <- list()
  pivot_location <- list()
  c <- 0
  
  # Step 1: Filter and order samples based on 3579 rule
  for (i in 1:nsample) {
    delta <- matrix(fit$delta[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                    nrow = p)
    lambda <- matrix(fit$lambda[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                     nrow = p)
    
    # Check if zero columns are present (not identifiable)
    if (sum(colSums(delta) == 0) == 0) {
      
      # If zero rows are present, check counting rule only on non-zero rows
      delta_tocheck <- as.matrix(delta[rowSums(delta) > 0, ])
      
      # Check 3579 rule
      if (counting_rule_holds(delta_tocheck)) {
        c <- c + 1
        
        # Order pivots decreasing and sign switch when needed
        neword <- reorder_and_sign_swap(delta, lambda)
        lambda_ordered[[c]] <- neword$lambda
        delta_ordered[[c]] <- neword$delta
        pivot_location[[c]] <- neword$pivot_location
      }
    }
  }
  
  # Step 2: Find modal pivot configuration
  modal_config <- find_modal_vector(pivot_location)
  modal_rank <- length(modal_config$mode)
  
  # Step 3: Keep only lambdas with modal configuration
  c <- 0
  lambda_identified <- list()
  delta_identified <- list()
  
  for (i in 1:length(lambda_ordered)) {
    if (length(pivot_location[[i]]) == modal_rank) {
      if (all(pivot_location[[i]] == modal_config$mode)) {
        c <- c + 1
        lambda_identified[[c]] <- lambda_ordered[[i]]
        delta_identified[[c]] <- delta_ordered[[i]]
      }
    }
  }
  
  # Prepare output
  result <- list(
    lambda = lambda_identified,
    modal_config = modal_config$mode,
    modal_rank = modal_rank,
    n_identified = length(lambda_identified)
  )
  
  if (return_delta) {
    result$delta <- delta_identified
  }
  
  if (return_pivot_info) {
    result$all_pivot_locations <- pivot_location
    result$lambda_ordered <- lambda_ordered
    result$delta_ordered <- delta_ordered
  }
  
  return(result)
}



#------------------------------------------------------------------------------#
# # function to simulate synthetic data
# simulate_data = function(n, p, k, covariate, sd_gamma, cp, sigma_sq, mseed){
#   set.seed(mseed)
#   q = dim(covariate)[2]
#   gamma = matrix(rnorm(q*p, sd = sd_gamma), q,p) # regression coefficient
#   pilogit = cp*plogis(covariate %*% gamma) 
#   get_phi = simPhi(p, k, pilogit)
#   Phi = get_phi$Phi
#   theta_inv = rgamma(k,1,1)
#   Lambda = matrix(0,p,k)
#   for(j in 1:p){
#     for(h in 1:k){
#       Lambda[j,h]= Phi[j,h]/theta_inv[h]
#     }
#   }
#   M = matrix(rnorm(n*k, 0, 1), ncol=k) # factor scores
#   # data
#   Y = M %*% t(Lambda) + sqrt(sigma_sq)*matrix(rnorm(n*p), nrow = n, ncol = p) 
#   return(list("Y" = Y, "Lambda" = Lambda, "Phi"=Phi, "pivot"=get_phi$pivot,
#               "delta"=get_phi$Delta, "gamma"=gamma))
# }
# 
# # function for simulate Phi via gibbs-sampler
# simPhi = function(p, k, pilogit, NiterGibbs = 200){
#   # Initialization 
#   phi = matrix(rbinom(p*k, 1, 0.2), p, k) 
#   delta = phi
#   l=rep(NA, k) # pivots
#   l[1] = which(phi[,1]>0)[1]
#   for(h in 2:k){
#     w = which(phi[,h]>0) # counting 1's in the h column of Delta
#     l[h] = w[!w %in% l[-h]][1] # take the first
#     if(is.na(l[h])){
#       delta[,h] = 0 # put all column to zeros
#     }else{
#       delta[1:l[h],h] = 0 # before pivots all zeros
#       delta[l[h]:p, h] =  phi[l[h]:p,h] # from pivot and after equal to phi[,l]
#     }
#   }
#   
#   for(t in 2:NiterGibbs){
#     for(h in 1:k){ 
#       Lh = sort(l[-h])
#       p_h = pilogit[,h]
#       k = 0     
#       # definisco prob di phi_h condizionatamente agli altri pivot
#       for(pl in Lh){
#         k = k+1
#         if (pl-k>0){
#           scale = 1-prod(1-pilogit[1:pl,h][-Lh[1:k]])
#           p_h[pl] = min(p_h[pl]/scale, 1)
#         }
#         else{
#           p_h[pl] = 1
#         }
#       }
#       phi[,h] = rbinom(p,1,p_h) # generate phi
#       w = which(phi[,h]>0)
#       l[h] = w[!w %in% l[-h]][1]
#       if(is.na(l[h])){
#         delta[,h] = 0
#       }else{
#         delta[1:l[h],h] = 0
#         delta[l[h]:p, h] =  phi[l[h]:p,h]
#       }
#     }
#   }
#   # return the last phi and delta
#   return(list("Phi" = phi, "Delta"= delta, "pivot"=l))
# }


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
