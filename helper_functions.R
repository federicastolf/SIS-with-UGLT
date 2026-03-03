library(MASS)

# simulate data with block structure covariance
simulate_block_mvn <- function(n, p, n_blocks, within_cor, between_cor, variances, 
                               mseed) {
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
  
  D <- diag(rep(sqrt(variances), p))
  Sigma1 <- D %*% Sigma %*% D
  # Simulate data
  Y <- mvrnorm(n = n, mu = rep(0, p), Sigma = Sigma)
  xlab <- as.factor(rep(1:n_blocks, each = block_size))
  
  return(list(data = Y, Sigma = Sigma1, xlab = xlab))
}

# new
simulate_factor_cov_block <- function(p,
                                      K,
                                      n_blocks = 3,
                                      block_sizes = NULL,      # length n_blocks; defaults to equal split of p
                                      K_block = NULL,          # length n_blocks; defaults to split K across blocks
                                      loading_sd = 1,
                                      psi = NULL,              # NULL -> no diagonal; scalar or length p vector allowed
                                      cross_block_sd = 0,      # >0 gives small cross-block loadings
                                      seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # block sizes over variables
  if (is.null(block_sizes)) {
    base <- p %/% n_blocks
    rem  <- p %% n_blocks
    block_sizes <- rep(base, n_blocks)
    if (rem > 0) block_sizes[seq_len(rem)] <- block_sizes[seq_len(rem)] + 1
  }
  stopifnot(sum(block_sizes) == p)
  
  # factors per block
  if (is.null(K_block)) {
    base <- K %/% n_blocks
    rem  <- K %% n_blocks
    K_block <- rep(base, n_blocks)
    if (rem > 0) K_block[seq_len(rem)] <- K_block[seq_len(rem)] + 1
  }
  stopifnot(sum(K_block) == K)
  
  # variable indices by block
  idx_blocks <- split(seq_len(p), rep(seq_len(n_blocks), times = block_sizes))
  
  # factor indices by block
  fac_blocks <- split(seq_len(K), rep(seq_len(n_blocks), times = K_block))
  
  # build Lambda (p x K)
  Lambda <- matrix(0, nrow = p, ncol = K)
  
  for (b in seq_len(n_blocks)) {
    rows <- idx_blocks[[b]]
    cols <- fac_blocks[[b]]
    
    # main within-block loadings
    Lambda[rows, cols] <- matrix(
      rnorm(length(rows) * length(cols), mean = 0, sd = loading_sd),
      nrow = length(rows), ncol = length(cols)
    )
    
    # optional small cross-block loadings (adds weak correlation across blocks)
    if (cross_block_sd > 0) {
      other_cols <- setdiff(seq_len(K), cols)
      if (length(other_cols) > 0) {
        Lambda[rows, other_cols] <- matrix(
          rnorm(length(rows) * length(other_cols), mean = 0, sd = cross_block_sd),
          nrow = length(rows), ncol = length(other_cols)
        )
      }
    }
  }
  
  # covariance from factor part
  Sigma <- Lambda %*% t(Lambda)
  
  # optional diagonal noise Psi
  if (!is.null(psi)) {
    if (length(psi) == 1) psi <- rep(psi, p)
    stopifnot(length(psi) == p)
    Sigma <- Sigma + diag(psi)
  }
  
  # --- categorical block label ---
  xlab <- factor(rep(seq_len(n_blocks), times = block_sizes))
  
  list(
    Sigma = Sigma,
    Lambda = Lambda,
    xlab = xlab,
    block_sizes = block_sizes,
    K_block = K_block
  )
}


simulate_factor_block_with_cross <- function(p,
                                             block_sizes,
                                             block_rho,
                                             cross_rho = 0.1,
                                             diag_var = 1,
                                             loading_noise_sd = 0,
                                             psi_jitter_sd = 0,
                                             seed = NULL) {
  
  if (!is.null(seed)) set.seed(seed)
  
  n_blocks <- length(block_sizes)
  stopifnot(sum(block_sizes) == p)
  stopifnot(length(block_rho) == n_blocks)
  stopifnot(cross_rho < min(block_rho))
  
  # One factor per block + 1 global factor
  K <- n_blocks + 1
  Lambda <- matrix(0, nrow = p, ncol = K)
  
  # ---- xlab (block membership for variables) ----
  xlab <- factor(rep(seq_len(n_blocks), times = block_sizes))
  
  # Global loading controls cross-block correlation
  l_global <- sqrt(cross_rho)
  
  # Block-specific loadings
  l_block <- sqrt(pmax(block_rho - cross_rho, 0))
  
  for (b in seq_len(n_blocks)) {
    rows <- which(xlab == b)
    
    # block factor
    Lambda[rows, b] <- l_block[b] +
      rnorm(length(rows), 0, loading_noise_sd)
    
    # global factor (last column)
    Lambda[rows, K] <- l_global +
      rnorm(length(rows), 0, loading_noise_sd)
  }
  
  Sigma_factor <- Lambda %*% t(Lambda)
  
  psi <- diag_var - rowSums(Lambda^2)
  psi <- pmax(psi, 1e-8)
  
  if (psi_jitter_sd > 0) {
    psi <- pmax(psi + rnorm(p, 0, psi_jitter_sd), 1e-8)
  }
  
  Sigma <- Sigma_factor + diag(psi)
  
  list(
    Sigma = Sigma,
    Lambda = Lambda,
    psi = psi,
    xlab = xlab,
    block_sizes = block_sizes
  )
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
