// ============================================================================
// FUNCTION 1: Calculate prior probability varpi_jh
// ============================================================================
// Calcola la probabilità prior che il pivot sia prima della posizione j
// varpi_jh = pr(l_h < j | l_{r,-h})
// = 1 - prod_{m in L(l_{r,-h}), m<j} {1 - logit^{-1}(w_m^T gamma_h) * c_p}
//
// Arguments:
//   j: row index
//   h: column index
//   Delta: indicator matrix for pivots in other columns
//   logit: matrix of logit probabilities
//   p_constant: constant c_p
//
// Returns: varpi_jh probability
double calculate_varpi(int j, int h, const arma::mat& Delta, 
                       const arma::mat& logit, double p_constant) {
  double prod = 1.0;
  
  // Loop over all m < j that are in L(l_{r,-h})
  for (int m = 0; m < j; m++) {
    // Check if m is in L(l_{r,-h}): no other column has pivot at m
    bool in_L = true;
    for (int hh = 0; hh < Delta.n_cols; hh++) {
      if (hh != h && Delta(m, hh) == 1) {
        in_L = false;
        break;
      }
    }
    
    if (in_L) {
      prod *= (1.0 - logit(m, h) * p_constant);
    }
  }
  
  return 1.0 - prod;
}


// Update phi_jh for j in L(l_{r,-h}) - potential pivot positions
// These are positions that are NOT pivots in other columns
//
// Arguments:
//   h: column index to update
//   Phi: matrix of local scales (modified in place)
//   Delta: indicator matrix (current pivot locations)
//   rho: vector indicating active/inactive factors
//   logit: matrix of logit probabilities
//   p_constant: constant c_p
//   epsilon: residual matrix
//   eta: factor scores
//   Lambda_star: loadings matrix (tilde lambda)
//   ps: precision vector
//   n: number of observations
//   p: number of variables
//
void update_phi_potential_pivots(int h, arma::mat& Phi, const arma::mat& Delta,
                                  const arma::vec& rho, const arma::mat& logit,
                                  double p_constant, const arma::mat& epsilon,
                                  const arma::mat& eta, const arma::mat& Lambda_star,
                                  const arma::vec& ps, int n, int p) {
  
  // Precompute which rows are in L(l_{r,-h})
  arma::vec in_L_vec = arma::sum(Delta, 1) - Delta.col(h);
  
  for (int j = 0; j < p; j++) {
    // Skip if j is NOT in L(l_{r,-h})
    if (in_L_vec(j) > 0) continue;
    
    double prob_phi_1;  // Probability that phi_jh = 1
    
    // Case 1: Inactive factor (rho_h = 0) - sample from prior
    if (rho(h) == 0) {
      prob_phi_1 = logit(j, h) * p_constant;
    }
    // Case 2: Active factor - use full conditional
    else {
      double prior_prob = logit(j, h) * p_constant;
      arma::vec resid_without_h = epsilon.row(j).t();
      for (int r = 0; r < Phi.n_cols; r++) {
        if (r != h) {
          resid_without_h -= eta.col(r) * Lambda_star(j, r) * Phi(j, r);
        }
      }
      
      // Log-likelihood for phi_jh = 0
      double sse_0 = arma::accu(arma::pow(resid_without_h, 2));
      double log_lik_0 = -0.5 * ps(j) * sse_0;
      
      // Log-likelihood for phi_jh = 1
      arma::vec resid_with_h = resid_without_h - eta.col(h) * Lambda_star(j, h);
      double sse_1 = arma::accu(arma::pow(resid_with_h, 2));
      double log_lik_1 = -0.5 * ps(j) * sse_1;
      
      // Posterior probability using Bayes rule
      double log_prob_0 = std::log(1.0 - prior_prob) + log_lik_0;
      double log_prob_1 = std::log(prior_prob) + log_lik_1;
      
      // Normalize (log-sum-exp trick)
      double max_log = std::max(log_prob_0, log_prob_1);
      double sum_exp = std::exp(log_prob_0 - max_log) + std::exp(log_prob_1 - max_log);
      prob_phi_1 = std::exp(log_prob_1 - max_log) / sum_exp;
    }
    
    // Sample phi_jh
    Phi(j, h) = (R::runif(0, 1) < prob_phi_1) ? 1.0 : 0.0;
  }
}

// ============================================================================
// FUNCTION 3: Update remaining phi_jh (non-pivot positions)
// ============================================================================
// Update phi_jh for j NOT in L(l_{r,-h}) - positions that CANNOT be pivots
// These are positions that ARE pivots in other columns
//
// Arguments:
//   h: column index to update
//   l_h: current pivot location for column h
//   Phi: matrix of local scales (modified in place)
//   Delta: indicator matrix (current pivot locations)
//   rho: vector indicating active/inactive factors
//   logit: matrix of logit probabilities
//   p_constant: constant c_p
//   epsilon: residual matrix
//   eta: factor scores
//   Lambda_star: loadings matrix (tilde lambda)
//   ps: precision vector
//   n: number of observations
//   p: number of variables
//
void update_phi_remaining(int h, int l_h, arma::mat& Phi, const arma::mat& Delta,
                          const arma::vec& rho, const arma::mat& logit,
                          double p_constant, const arma::mat& epsilon,
                          const arma::mat& eta, const arma::mat& Lambda_star,
                          const arma::vec& ps, int n, int p) {
  
  // Precompute which rows are NOT in L(l_{r,-h})
  arma::vec not_in_L_vec = arma::sum(Delta, 1) - Delta.col(h);
  
  for (int j = 0; j < p; j++) {
    // Skip if j IS in L(l_{r,-h}) (already updated in previous step)
    if (not_in_L_vec(j) == 0) continue;
    
    // Calculate varpi_jh
    double varpi_jh = calculate_varpi(j, h, Delta, logit, p_constant);
    
    // Prior probability with varpi adjustment
    double prior_prob = std::min(p_constant * logit(j, h) / varpi_jh, 1.0);
    
    double prob_phi_1;  // Probability that phi_jh = 1
    
    // Case 1: Inactive factor (rho_h = 0) - sample from prior
    if (rho(h) == 0) {
      prob_phi_1 = prior_prob;
    }
    // Case 2: Active factor and l_h > j - sample from prior
    // (because delta_jh = 0 always when l_h > j)
    else if (l_h > j) {
      prob_phi_1 = prior_prob;
    }
    // Case 3: Active factor and j > l_h - use full conditional
    else {
      // Compute residual without factor h
      arma::vec resid_without_h = epsilon.row(j).t();
      for (int r = 0; r < Phi.n_cols; r++) {
        if (r != h) {
          resid_without_h -= eta.col(r) * Lambda_star(j, r) * Phi(j, r);
        }
      }
      
      // Log-likelihood for phi_jh = 0
      double sse_0 = arma::accu(arma::pow(resid_without_h, 2));
      double log_lik_0 = -0.5 * ps(j) * sse_0;
      
      // Log-likelihood for phi_jh = 1
      arma::vec resid_with_h = resid_without_h - eta.col(h) * Lambda_star(j, h);
      double sse_1 = arma::accu(arma::pow(resid_with_h, 2));
      double log_lik_1 = -0.5 * ps(j) * sse_1;
      
      // Posterior probability using Bayes rule
      double log_prob_0 = std::log(1.0 - prior_prob) + log_lik_0;
      double log_prob_1 = std::log(prior_prob) + log_lik_1;
      
      // Normalize (log-sum-exp trick)
      double max_log = std::max(log_prob_0, log_prob_1);
      double sum_exp = std::exp(log_prob_0 - max_log) + std::exp(log_prob_1 - max_log);
      prob_phi_1 = std::exp(log_prob_1 - max_log) / sum_exp;
    }
    
    // Sample phi_jh
    Phi(j, h) = (R::runif(0, 1) < prob_phi_1) ? 1.0 : 0.0;
  }
}

// ============================================================================
// FUNCTION 4: Update Delta matrix based on Phi and pivot location
// ============================================================================
// Update delta_jh: 
//   - delta_jh = 0 for all j < l_h
//   - delta_jh = phi_jh for all j >= l_h
//
// Arguments:
//   h: column index to update
//   l_h: pivot location for column h (-1 if no pivot)
//   Phi: matrix of local scales
//   Delta: indicator matrix (modified in place)
//   p: number of variables
//
void update_delta_column(int h, int l_h, const arma::mat& Phi, 
                         arma::mat& Delta, int p) {
  
  // If no pivot found, set entire column to 0
  if (l_h < 0) {
    Delta.col(h).zeros();
    return;
  }
  
  // Set delta_jh = 0 for all j < l_h
  for (int j = 0; j < l_h; j++) {
    Delta(j, h) = 0.0;
  }
  
  // Set delta_jh = phi_jh for all j >= l_h
  for (int j = l_h; j < p; j++) {
    Delta(j, h) = Phi(j, h);
  }
}