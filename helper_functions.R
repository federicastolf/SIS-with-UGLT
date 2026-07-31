library(MASS)

# simulate data with block structure covariance
simulate_block_mvn = function(n, p, n_blocks, within_cor, between_cor, variances = 1,
                              mseed, count = FALSE, lambda = 5) {
  set.seed(mseed)
  block_size = p / n_blocks
  
  if (length(within_cor) == 1) within_cor = rep(within_cor, n_blocks)
  if (length(within_cor) != n_blocks) stop("within_cor must be length 1 or n_blocks")
  
  if (length(variances) == 1) variances = rep(variances, p)
  if (length(variances) != p) stop("variances must be length 1 or p")
  
  if (length(lambda) == 1) lambda = rep(lambda, p)
  if (length(lambda) != p) stop("lambda must be length 1 or p")
  
  # matrice di correlazione a blocchi
  R = matrix(between_cor, p, p)
  for (i in 1:n_blocks) {
    idx = ((i - 1) * block_size + 1):(i * block_size)
    R[idx, idx] = within_cor[i]
  }
  diag(R) = 1
  
  D = diag(sqrt(variances))
  Sigma = D %*% R %*% D
  
  if (!is.null(mseed)) set.seed(mseed)
  Z = MASS::mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  
  if (count) {
    sds = sqrt(variances)
    U = pnorm(Z, mean = 0, sd = rep(sds, each = n))
    Y = sapply(seq_len(p), function(j) qpois(U[, j], lambda = lambda[j]))
  } else {
    Y = Z
  }
  xlab = as.factor(rep(1:n_blocks, each = block_size))
  return(list(data = Y, Sigma = Sigma, xlab = xlab))
}


reorder_blocks <- function(Sigma, Y, block_membership, new_block_order) {
  p <- ncol(Y)
  n_blocks <- length(unique(block_membership))
  row_order <- unlist(lapply(new_block_order, function(b) {
    which(block_membership == b)
  }))
  Sigma_reordered <- Sigma[row_order, row_order]
  Y_reordered <- Y[, row_order]
  block_membership_reordered <- block_membership[row_order]
  
  return(list(Sigma = Sigma_reordered,
              Y = Y_reordered, xlab = block_membership_reordered, 
              row_order = row_order, new_block_order = new_block_order
  ))
}


identify_lambda_eta <- function(fit, p, nsample = NULL, 
                            return_delta = TRUE,
                            return_eta = TRUE,
                            return_pivot_info = FALSE) {
  if (is.null(nsample)) {
    nsample <- length(fit$lambda)
  }
  
  # Initialize lists
  n = nrow(fit$eta[[1]])
  lambda_ordered <- list()
  eta_ordered <- list()
  delta_ordered <- list()
  pivot_location <- list()
  to_keep = rep(0, nsample)
  c <- 0
  
  # Step 1: Filter and order samples based on 3579 rule
  for (i in 1:nsample) {
    delta <- matrix(fit$delta[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                    nrow = p)
    lambda <- matrix(fit$lambda[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                     nrow = p)
    eta <- matrix(fit$eta[[i]][, c(fit$activeFactors[[i]] == TRUE)],
                     nrow = n)
    
    # Check if zero columns are present (not identifiable)
    if (sum(colSums(delta) == 0) == 0) {
      
      # If zero rows are present, check counting rule only on non-zero rows
      delta_tocheck <- as.matrix(delta[rowSums(delta) > 0, ])
      
      # Check 3579 rule
      if (counting_rule_holds(delta_tocheck)) {
        c <- c + 1
        
        # Order pivots decreasing and sign switch when needed
        neword <- reorder_and_sign_swap(delta, lambda, eta)
        lambda_ordered[[c]] <- neword$lambda
        delta_ordered[[c]] <- neword$delta
        eta_ordered[[c]] <- neword$eta
        pivot_location[[c]] <- neword$pivot_location
        to_keep[i] = 1
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
  eta_identified <- list()
  id_tmp = which(to_keep>0)
  sample_index = NULL

  for (i in 1:length(lambda_ordered)) {
    if (length(pivot_location[[i]]) == modal_rank) {
      if (all(pivot_location[[i]] == modal_config$mode)) {
        c <- c + 1
        lambda_identified[[c]] <- lambda_ordered[[i]]
        delta_identified[[c]] <- delta_ordered[[i]]
        eta_identified[[c]] <- eta_ordered[[i]]
        sample_index = c(sample_index, id_tmp[i])
      }
    }
  }
  
  # Prepare output
  result <- list(
    lambda = lambda_identified,
    modal_config = modal_config$mode,
    modal_rank = modal_rank,
    n_identified = length(lambda_identified),
    sample_index = sample_index
  )
  
  if (return_delta) {
    result$delta <- delta_identified
  }
  
  if (return_eta) {
    result$eta <- eta_identified
  }
  
  if (return_pivot_info) {
    result$all_pivot_locations <- pivot_location
    result$lambda_ordered <- lambda_ordered
    result$delta_ordered <- delta_ordered
  }
  
  return(result)
}


# function to reorder and sign swap based on delta and lambda
reorder_and_sign_swap <- function(delta, lambda, eta) {
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
  eta_new <- eta[, order_idx, drop = FALSE]
  pivot_sign <- pivot_sign[order_idx]
  
  # Flip sign where needed
  for (h in seq_len(k)) {
    if (pivot_sign[h] < 0) {
      lambda_new[, h] <- -lambda_new[, h]
      eta_new[, h] <- -eta_new[, h]
    }
  }
  
  list(delta = delta_new, lambda = lambda_new, eta = eta_new,
       order = order_idx,
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


# Create pairly negative correlated beta covariance
# rho < 0: negative correlation
# rho > -1/(q-1) for positive definiteness
make_beta_cov <- function(sd_beta, rho) {
  q <- length(sd_beta)
  
  if (q > 1L && (rho <= -1 / (q - 1) || rho >= 1)) {
    stop("rho must satisfy -1/(q-1) < rho < 1")
  }
  
  R <- matrix(rho, q, q)
  diag(R) <- 1
  
  D <- diag(sd_beta, q)
  Sigma_beta <- D %*% R %*% D
  
  if (min(eigen(Sigma_beta, symmetric = TRUE,
                only.values = TRUE)$values) <= 0) {
    stop("Sigma_beta is not positive definite")
  }
  
  Sigma_beta
}
# Example
# q = 3 # number meta-covariate.rho <0 and rho> -1/(3-1) -> rho>-0.5 
# Sigma_gamma = make_beta_cov(sd_beta = rep(1,3), rho = -1/2+0.01)