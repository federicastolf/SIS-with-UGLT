#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]


/* -------------------------------------------------------------------------- */

arma::vec join_elem(const arma::vec& v1, const int& v2) {
  return arma::join_cols(v1, v2 * arma::ones(1));
}

arma::vec join_elem(const arma::vec& v1, const double& v2) {
  return arma::join_cols(v1, v2 * arma::ones(1));
}

arma::vec join_elem(const int& v1, const arma::vec& v2) {
  return arma::join_cols(v1 * arma::ones(1), v2);
}

arma::vec join_elem(const double& v1, const arma::vec& v2) {
  return arma::join_cols(v1 * arma::ones(1), v2);
}

/* -------------------------------------------------------------------------- */

arma::mat rbinom_vec(const int& len, const int& size, const double& prob) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::rbinom(size, prob);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat rnorm_mat(const int& rows, const int& cols, const double& mean, const double& sd) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::rnorm(mean, sd);
    }
  }
  return M;
}

arma::mat rnorm_vec(const int& len, const double& mean, const double& sd) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::rnorm(mean, sd);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat runif_mat(const int& rows, const int& cols, const double& minVal, const double& maxVal) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::runif(minVal, maxVal);
    }
  }
  return M;
}

arma::mat runif_vec(const int& len, const double& minVal, const double& maxVal) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::runif(minVal, maxVal);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat mvrnormArma(int n, const arma::vec& mu, const arma::mat& Sigma) {
  int d = mu.n_elem;
  arma::mat L = arma::chol(Sigma, "lower");
  arma::mat Z = arma::randn(n, d);
  arma::mat X = arma::repmat(mu.t(), n, 1) + Z * L.t();
  return X;
}

/* -------------------------------------------------------------------------- */

double dmvnorm_log(const arma::vec& x, const arma::vec& mu, const arma::mat& Sigma) {
  int k = mu.n_elem;
  arma::vec diff = x - mu;
  double log_det_val;
  double sign;
  arma::log_det(log_det_val, sign, Sigma);
  double quad_form = arma::as_scalar(diff.t() * arma::inv(Sigma) * diff);
  double log_dens = -0.5 * (k * std::log(2.0 * M_PI) + log_det_val + quad_form);
  return log_dens;
}

/* -------------------------------------------------------------------------- */
// Functions in this file:
//  - is_pivot_row
//  - update_bmu
//  - update_mu
//  - update_eta
//  - update_beta
//  - update_Lambda_star
//  - update_d
//  - update_Phi
/* -------------------------------------------------------------------------- */

// Update eta in the Adaptive Gibbs Sampler
//
// @param Lambda A pxk matrix.
// @param ps A p-dimensional vector.
// @param Z A nxp matrix.
//
// @return A nxk matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_eta(const arma::mat& Lambda, const arma::vec& ps, const arma::mat& Z) {
  int k = Lambda.n_cols;
  int n = Z.n_rows;
  arma::mat Lmsg = Lambda.each_col() % ps;
  arma::mat Veta1(k, k);
  if (k > 1) {
    Veta1 = arma::diagmat(arma::ones(k)) + trans(Lmsg) * Lambda;
  } else {
    Veta1 = 1 + trans(Lmsg) * Lambda;
  }
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Veta1);
  eigval = arma::reverse(eigval);
  eigvec = arma::fliplr(eigvec);
  arma::mat Tmat(k, k);
  if (arma::min(eigval) > 0.000001) {
    Tmat = trans(eigvec.each_row() % trans(sqrt(eigval)));
  } else {
    Tmat = arma::chol(Veta1);
  }
  arma::mat Q;
  arma::mat R;
  arma::qr(Q, R, Tmat);
  arma::mat S = arma::inv(R);
  arma::mat eta = Z * Lmsg * S * trans(S) + rnorm_mat(n, k, 0, 1) * trans(S);
  return eta;
}

/* -------------------------------------------------------------------------- */

// @param j An integer number.
// @param hh An integer number.
// @param Delta a binary matrix.
// A row j is the pivot of column hh iff it's the first 1 in Delta(:,hh)
inline bool is_pivot_row(int j, int hh, const arma::mat& Delta) {
  if (Delta(j, hh) == 0) return false;
  if (j == 0) return true;
  return Delta(j - 1, hh) == 0;
}


/* -------------------------------------------------------------------------- */

// Update the jth row of Lambda_star in the Adaptive Gibbs Sampler
//
// @param j An integer number.
// @param etarho A kxn matrix.
// @param Phi A pxk matrix.
// @param Plam A kxk matrix.
// @param ps A p-dimensional vector;
// @param Z A nxp matrix.
//
// @return A 1xk matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_Lambda_star(const int& j, const arma::mat& etarho, const arma::mat& Phi,
                             const arma::mat& Plam, const arma::vec& ps, const arma::mat& Z) {
  int k = Phi.n_cols;
  arma::mat etaj = trans(etarho.each_col() % trans(Phi.row(j)));
  arma::mat Qlam = Plam + ps(j) * trans(etaj) * etaj;
  arma::mat Llam = trans(arma::chol(Qlam));
  arma::vec blam = ps(j) * trans(etaj) * Z.col(j);
  // mean
  arma::mat vlam = solve(trimatl(Llam), blam);
  arma::mat mlam = solve(trimatu(trans(Llam)), vlam);
  // var
  arma::vec zlam = rnorm_vec(k, 0, 1);
  arma::mat ylam = solve(trimatu(trans(Llam)), zlam);
  arma::mat lambda_starj = trans(ylam + mlam);
  return lambda_starj;
}

/* -------------------------------------------------------------------------- */

// Update the hth element of d in the Adaptive Gibbs Sampler
//
// @param h An integer.
// @param Phi A pxk matrix.
// @param rho A k-dimensional vector.
// @param eta A nxk matrix.
// @param lambdastar A pxk matrix.
// @param Z A nxp matrix.
// @param ps A p-dimensional vector.
// @param w A k-dimensional vector.
//
// @return An integer in 1, ..., k.
//
// @note This function uses \code{Rcpp} for computational efficiency.
int update_d(const int& h, const arma::mat& Phi, const arma::vec& rho,
             const arma::mat& eta, const arma::mat& lambdastar,
             const arma::mat& Z, const arma::vec& ps, const arma::vec& w) {
  int k = Phi.n_cols;
  int n = eta.n_rows;
  int p = Phi.n_rows;
  int i, j, l;
  double lnorm0 = 0.0, lnorm1 = 0.0;
  arma::vec sdy = sqrt(1 / ps);
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      // initialize mean
      double muijh = -1 * sqrt(rho(h)) * Phi(j, h) * lambdastar(j, h) * eta(i, h);
      for (l = 0; l < k; l++) {
        muijh += sqrt(rho(l)) * Phi(j, l) * lambdastar(j, l) * eta(i, l);
      }
      lnorm0 += R::dnorm(Z(i, j), muijh, sdy(j), 1);
      // update mean
      muijh += Phi(j, h) * lambdastar(j, h) * eta(i, h);
      lnorm1 += R::dnorm(Z(i, j), muijh, sdy(j), 1);
    }
  }
  // adjust the scale
  double mlnorm = std::max(lnorm0, lnorm1);
  lnorm0 -= mlnorm;
  lnorm1 -= mlnorm;
  Rcpp::NumericVector prob_h = as<NumericVector>(wrap(log(w)));
  for (i = 0; i <= h; i++) {
    prob_h(i) += lnorm0;
  }
  for (i = h + 1; i < k; i++) {
    prob_h(i) += lnorm1;
  }
  prob_h = exp(prob_h);
  if (sum(prob_h) == 0.0) {
    prob_h = rep(0, k);
    prob_h(k - 1) = 1.0;
  }
  prob_h = prob_h / sum(prob_h);
  // draw d from a multinomial distribution
  Rcpp::IntegerVector d(k);
  R::rmultinom(1, prob_h.begin(), k, d.begin());
  return Rcpp::which_max(d);
}

/* -------------------------------------------------------------------------- */
// Calculate  varpi_jh = pr(l_h < j | l_{r,-h})
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
      //if (hh != h && Delta(m, hh) == 1) {
      if (hh != h && is_pivot_row(m, hh, Delta)){
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


/* -------------------------------------------------------------------------- */
// Identify the pivot l_h as the first non-zero element in column h of Phi
// that is not a pivot in other columns
// Returns: pivot location l_h (or -1 if no pivot found)
int identify_pivot(int h, const arma::mat& Phi, const arma::mat& Delta, int p) {
  for (int j = 0; j < p; j++) {
    // Check if j is in L(l_{r,-h}): not a pivot in other columns
    bool in_L = true;
    for (int hh = 0; hh < Delta.n_cols; hh++) {
      //if (hh != h && Delta(j, hh) == 1) {
      if (hh != h && is_pivot_row(j, hh, Delta)) {
        in_L = false;
        break;
      }
    }
    // If j is in L and phi_jh = 1, it's the pivot
    if (in_L && Phi(j, h) == 1) {
      return j;
    }
  }
  return -1;  // No pivot found
}

/* -------------------------------------------------------------------------- */
// Update Delta based on pivot location
//   - delta_jh = 0 for all j < l_h
//   - delta_jh = phi_jh for all j >= l_h
void update_delta_column(int h, int l_h, const arma::mat& Phi, arma::mat& Delta) {
  int p = Phi.n_rows;
  // Copy from Phi
  for (int j = 0; j < p; j++) {
    Delta(j, h) = Phi(j, h);
  }
  // Set delta_jh = 0 for all j < l_h
  if (l_h >= 0) {
    for (int j = 0; j < l_h; j++) {
      Delta(j, h) = 0.0;
    }
  }
}

/* -------------------------------------------------------------------------- */
// Compute L(l_{r,-h}) membership for all columns (used for MH step)
// Returns a matrix where in_L(j,h) = 1 if j is in L(l_{r,-h}), 0 otherwise
arma::mat compute_L_membership(const arma::mat& Delta, int p, int k) {
  arma::mat in_L = arma::ones(p, k);
  for (int h = 0; h < k; h++) {
    for (int j = 0; j < p; j++) {
      // Check if j is a pivot in any other column
      for (int hh = 0; hh < k; hh++) {
        
        //if (hh != h && Delta(j, hh) == 1) {
        if (hh != h && is_pivot_row(j, hh, Delta)) {
          in_L(j, h) = 0;
          break;
        }
        
      }
    }
  }
  return in_L;
}

/* -------------------------------------------------------------------------- */
// Compute L(l_{r,-h}) membership for a single column h
// Returns a vector where in_L_h(j) = 1 if j is in L(l_{r,-h}), 0 otherwise
arma::vec compute_L_membership_single(int h, const arma::mat& Delta, int p) {
  arma::vec in_L_h = arma::ones(p);
  
  for (int j = 0; j < p; j++) {
    // Check if j is a pivot in any other column (excluding h)
    for (int hh = 0; hh < Delta.n_cols; hh++) {
      
      //if (hh != h && Delta(j, hh) == 1) {
      if (hh != h && is_pivot_row(j, hh, Delta)) {
        in_L_h(j) = 0;
        break;
      }
   
    }
  }
  return in_L_h;
}

/* -------------------------------------------------------------------------- */
// Update Phi for a single column h - potential pivots only (j in L(l_{r,-h}))
void update_phi_potential_pivots_single(int h, arma::mat& Phi, const arma::vec& in_L_h,
                                         const arma::vec& rho, const arma::mat& logit, 
                                         double p_constant, const arma::mat& eta, 
                                         const arma::mat& lambdastar, const arma::mat& Z, 
                                         const arma::vec& ps) {
  int k = Phi.n_cols;
  int n = eta.n_rows;
  int p = Phi.n_rows;
  int i, j, l;
  arma::vec p_phi0(p), p_phi1(p), p_phi_sum(p), lnorm0(p), lnorm1(p);
  arma::vec sdy = sqrt(1 / ps);
  
  // Check if factor h is active or inactive
  if (rho(h) == 0) {
    // Inactive factor
    p_phi0 = 1 - logit.col(h) * p_constant;
    p_phi1 = logit.col(h) * p_constant;
    p_phi_sum = p_phi0 + p_phi1;
    
    for (j = 0; j < p; j++) {
      if (in_L_h(j) == 1) {
        Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
      }
    }
  } else {
    // Active factor
    lnorm0 = arma::zeros(p);
    lnorm1 = arma::zeros(p);
    for (j = 0; j < p; j++) {
      if (in_L_h(j) == 0) continue; // Skip if j not in L(l_{r,-h})
      for (i = 0; i < n; i++) {
        double muijh0 = 0.0;
        for (l = 0; l < k; l++) {
          if (l != h) {
            muijh0 += rho(l) * Phi(j, l) * lambdastar(j, l) * eta(i, l);
          }
        }
        lnorm0(j) += R::dnorm(Z(i, j), muijh0, sdy(j), 1);
        double muijh1 = muijh0 + rho(h) * lambdastar(j, h) * eta(i, h);
        lnorm1(j) += R::dnorm(Z(i, j), muijh1, sdy(j), 1);
      }
      double mlnorm = std::max(lnorm0(j), lnorm1(j));
      lnorm0(j) -= mlnorm;
      lnorm1(j) -= mlnorm;
    }
    
    p_phi0 = exp(lnorm0 + log(1 - logit.col(h) * p_constant));
    p_phi1 = exp(lnorm1 + log(logit.col(h) * p_constant));
    p_phi_sum = p_phi0 + p_phi1;
    // Update only for j in L(l_{r,-h})
    for (j = 0; j < p; j++) {
      if (in_L_h(j) == 1) {
        Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
// Update remaining Phi for a single column h (j NOT in L(l_{r,-h}))  
// positions that CANNOT be pivots. These are positions that ARE pivots in other columns
void update_phi_remaining_single(int h, arma::mat& Phi, const arma::vec& not_in_L_h,
                                  int l_h, const arma::vec& rho, 
                                  const arma::mat& logit, double p_constant, 
                                  const arma::mat& eta, const arma::mat& lambdastar, 
                                  const arma::mat& Z, const arma::vec& ps, 
                                  const arma::mat& Delta, double cMH, bool order_dependent) {
  
  int k = Phi.n_cols;
  int n = eta.n_rows;
  int p = Phi.n_rows;
  int i, j, l;
  arma::vec p_phi0(p), p_phi1(p), p_phi_sum(p), lnorm0(p), lnorm1(p);
  arma::vec sdy = sqrt(1 / ps);
  
  // Compute varpi for all j
  arma::vec varpi(p);
  for (j = 0; j < p; j++) {
    if (not_in_L_h(j) == 1) {
      if(order_dependent == TRUE) {
            varpi(j) = 1.0;
          } else {
            varpi(j) = calculate_varpi(j, h, Delta, logit, p_constant);
          }
    }
  }
  if (rho(h) == 0) {
    // Inactive factor
    for (j = 0; j < p; j++) {
      if (not_in_L_h(j) == 1) {
        double prior_prob = std::min(p_constant * logit(j, h) / varpi(j), 1.0 - cMH);
        //double prior_prob = p_constant*std::min( logit(j, h) / varpi(j), 1.0);
        p_phi0(j) = 1.0 - prior_prob;
        p_phi1(j) = prior_prob;
      } else {
        p_phi0(j) = 1.0;
        p_phi1(j) = 0.0;
      }
    }
    p_phi_sum = p_phi0 + p_phi1;
    // Update only for j NOT in L(l_{r,-h})
    for (j = 0; j < p; j++) {
      if (not_in_L_h(j) == 1) {
        Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
      }
    }
  } else {
    // Active factor
    lnorm0 = arma::zeros(p);
    lnorm1 = arma::zeros(p);
    
    for (j = 0; j < p; j++) {
      if (not_in_L_h(j) == 0) continue; // Skip if j in L(l_{r,-h})
      if (l_h < 0 || l_h > j) continue; // Will use prior only (likelihoods stay at 0)
      // Case: j > l_h - use full conditional
      for (i = 0; i < n; i++) {
        double muijh0 = 0.0;
        for (l = 0; l < k; l++) {
          if (l != h) {
            muijh0 += rho(l) * Delta(j, l) * lambdastar(j, l) * eta(i, l);
          }
        }
        lnorm0(j) += R::dnorm(Z(i, j), muijh0, sdy(j), 1);
        double muijh1 = muijh0 + rho(h) * lambdastar(j, h) * eta(i, h);
        lnorm1(j) += R::dnorm(Z(i, j), muijh1, sdy(j), 1);
      }
      double mlnorm = std::max(lnorm0(j), lnorm1(j));
      lnorm0(j) -= mlnorm;
      lnorm1(j) -= mlnorm;
    }
    
    for (j = 0; j < p; j++) {
      if (not_in_L_h(j) == 1) {
        double prior_prob = std::min(p_constant * logit(j, h) / varpi(j), 1.0-cMH);
        //double prior_prob = p_constant *std::min( logit(j, h) / varpi(j), 1.0);
        if (l_h < 0 || l_h > j) {
          p_phi0(j) = 1.0 - prior_prob;
          p_phi1(j) = prior_prob;
        } else {
          p_phi0(j) = exp(lnorm0(j) + log(1.0 - prior_prob));
          p_phi1(j) = exp(lnorm1(j) + log(prior_prob));
        }
      } else {
        p_phi0(j) = 1.0;
        p_phi1(j) = 0.0;
      }
    }
    p_phi_sum = p_phi0 + p_phi1;
    // Update only for j NOT in L(l_{r,-h})
    for (j = 0; j < p; j++) {
      if (not_in_L_h(j) == 1) {
        Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
      }
    }
  }
}

/* -------------------------------------------------------------------------- */
// Calculate varpi_jh for a single column (vector version)
// varpi_jh = pr(l_h < j | l_{r,-h})
// = 1 - prod_{m in L(l_{r,-h}), m<j} {1 - logit_h(m) * c_p}
double calculate_varpi_single(int j, int h, const arma::mat& Delta, 
                               const arma::vec& logit_h, double p_constant) {
  double prod = 1.0;
  // Loop over all m < j that are in L(l_{r,-h})
  for (int m = 0; m < j; m++) {
    // Check if m is in L(l_{r,-h}): no other column has pivot at m
    bool in_L = true;
    for (int hh = 0; hh < Delta.n_cols; hh++) {
      //if (hh != h && Delta(m, hh) == 1) {
      if (hh != h && is_pivot_row(m, hh, Delta)){
        in_L = false;
        break;
      }
    }
    if (in_L) {
      prod *= (1.0 - logit_h(m) * p_constant);
    }
  }
  return 1.0 - prod;
}

/* -------------------------------------------------------------------------- */
// compute expected values of Pólya-Gamma
// Formula: E[PG(1, c)] = tanh(|c|/2) / (2|c|)
// For c → 0, use the analytical limit: 0.25
arma::vec expected_polya_gamma(const arma::vec& pred) {
  int n = pred.n_elem;
  arma::vec exp_d(n);
  
  for (int i = 0; i < n; i++) {
    double abs_pred = std::abs(pred(i));

    if (abs_pred < 1e-8) {
      exp_d(i) = 0.25;  // lim_{c→0} tanh(c/2)/(2c) = 1/4
    } else {
      exp_d(i) = std::tanh(abs_pred / 2.0) / (2.0 * abs_pred);
    }
  }
  
  return exp_d;
}


/* -------------------------------------------------------------------------- */
// Metropolis-Hastings sampler for gamma_h (single column of Gamma matrix)
// Uses expected Polya-Gamma values for the proposal distribution
bool sample_gamma_MH(int h, arma::mat& Gamma, const arma::mat& Phi_L, 
                     const arma::mat& Phi, const arma::mat& in_L, 
                     const arma::mat& Delta, const arma::mat& wB, 
                     const arma::mat& Bh_1, double p_constant,
                     double scale_factor_MH, int p, double cMH) {
  
  arma::vec gamma_h_current = Gamma.col(h);
  
  // Compute expected Polya-Gamma values at current state
  arma::vec pred_current = wB * gamma_h_current;
  arma::vec exp_d_current = expected_polya_gamma(pred_current);
  
  // Proposal using expected Polya-Gamma values with tunable variance
  arma::vec kappa_h = Phi_L.col(h) - 0.5;
  arma::mat D_exp_current = arma::diagmat(exp_d_current);
  arma::mat V_prop_base = arma::inv(wB.t() * D_exp_current * wB + Bh_1);
  arma::mat V_prop = scale_factor_MH * V_prop_base;  // Scale variance for tuning
  arma::vec m_prop = V_prop_base * (wB.t() * kappa_h);
  arma::vec gamma_h_prop = mvrnormArma(1, m_prop, V_prop).t();
  
  // Log-likelihood ratio
  arma::vec pred_prop = wB * gamma_h_prop;
  arma::vec logit_current = arma::exp(pred_current) / (1 + arma::exp(pred_current));
  arma::vec logit_prop = arma::exp(pred_prop) / (1 + arma::exp(pred_prop));
  double log_lik_current = 0.0;
  double log_lik_prop = 0.0;
  
  for (int i = 0; i < p; i++) {
    // Use precomputed membership
    bool in_L_ih = (in_L(i, h) == 1);
    
    double prob_current, prob_prop;
    if (in_L_ih) {
      // For j in L(l_{r,-h}): pr(φ_jh=1) = c_p * logit^{-1}(w_j^T γ_h)
      prob_current = p_constant * logit_current(i);
      prob_prop = p_constant * logit_prop(i);
    } else {
       // For j NOT in L(l_{r,-h}): pr(φ_jh=1) = min{c_p * logit / varpi_jh, 1}
      double varpi_i_current = calculate_varpi_single(i, h, Delta, logit_current, p_constant);
      double varpi_i_prop = calculate_varpi_single(i, h, Delta, logit_prop, p_constant);
      prob_current = std::min(p_constant * logit_current(i) / varpi_i_current, 1.0-cMH);
      prob_prop = std::min(p_constant * logit_prop(i) / varpi_i_prop, 1.0-cMH);
      //prob_current = p_constant * std::min(logit_current(i) / varpi_i_current, 1.0);
      //prob_prop = p_constant * std::min(logit_prop(i) / varpi_i_prop, 1.0);
      // OLD no correction
      //prob_current = std::min(p_constant * logit_current(i), 1.0);
      //prob_prop = std::min(p_constant * logit_prop(i), 1.0);
    }
    
    if (Phi(i, h) == 1) {
    //if (Phi_L(i, h) == 1) {
      log_lik_current += std::log(prob_current);
      log_lik_prop += std::log(prob_prop);
    } else {
      log_lik_current += std::log(1 - prob_current);
      log_lik_prop += std::log(1 - prob_prop);
    }
  }
  
  // Log-prior ratio (Gaussian prior with precision Bh_1)
  double log_prior_current = -0.5 * arma::as_scalar(gamma_h_current.t() * Bh_1 * gamma_h_current);
  double log_prior_prop = -0.5 * arma::as_scalar(gamma_h_prop.t() * Bh_1 * gamma_h_prop);
  
  // Log-transition ratio (using expected Polya-Gamma values with scaling)
  arma::vec exp_d_prop = expected_polya_gamma(pred_prop);
  arma::mat D_exp_prop = arma::diagmat(exp_d_prop);
  
  arma::mat V_trans_to_current_base = arma::inv(wB.t() * D_exp_prop * wB + Bh_1);
  arma::mat V_trans_to_current = scale_factor_MH * V_trans_to_current_base;
  arma::vec m_trans_to_current = V_trans_to_current_base * (wB.t() * kappa_h);
  
  arma::mat V_trans_to_prop_base = arma::inv(wB.t() * D_exp_current * wB + Bh_1);
  arma::mat V_trans_to_prop = scale_factor_MH * V_trans_to_prop_base;
  arma::vec m_trans_to_prop = V_trans_to_prop_base * (wB.t() * kappa_h);

  // Log transition densities
  double log_trans_to_prop = dmvnorm_log(gamma_h_prop, m_trans_to_prop, V_trans_to_prop);
  //double log_trans_to_prop = dmvnorm_log(gamma_h_prop, m_prop, V_prop);
  double log_trans_to_current = dmvnorm_log(gamma_h_current, m_trans_to_current, V_trans_to_current);
  
  // Log acceptance ratio
  double log_alpha = (log_lik_prop - log_lik_current) + 
                     (log_prior_prop - log_prior_current) +
                     (log_trans_to_current - log_trans_to_prop);
  double alpha = std::min(1.0, std::exp(log_alpha));
  
  if (R::runif(0, 1) < alpha) {
    Gamma.col(h) = gamma_h_prop;
    return true;
  }
  return false;
}



// Optimized Adaptive Gibbs Sampler - computes pivots once per iteration
// [[Rcpp::export]]
Rcpp::List Rcpp_gibbs(double alpha, double a_sigma, double b_sigma, double a_theta, 
  double b_theta, double sd_gammaB, double p_constant,
  arma::mat y, 
  arma::mat wB,
  int burn, int nrun, int thin, int start_adapt, int kmax,
  arma::mat eta, arma::mat Gamma, arma::mat Lambda, arma::mat Lambda_star, 
  arma::vec d, int kstar,
  arma::mat logit, arma::vec rho, arma::mat Phi, arma::mat Plam, 
  arma::mat pred, arma::vec ps,
  arma::vec v, arma::vec w,
  Rcpp::List out, bool verbose,
  arma::vec uu, arma::vec prob, int sp,
  arma::vec pivots, arma::mat Delta, double scale_factor_MH, double cMH,
  bool order_dependent) {
  
  // ---------------------------------------------------------------------------
  // output
  Rcpp::List GAMMA(sp);
  Rcpp::List ETA(sp);
  Rcpp::List LAMBDA(sp);
  Rcpp::List SIG(sp);
  Rcpp::List DELTA(sp);
  Rcpp::List RHO(sp);
  arma::vec K(sp);
  arma::vec ACC_RATE(sp);
  
  // ---------------------------------------------------------------------------
  // matrix dimensions
  int k = Lambda_star.n_cols;
  int n = y.n_rows;
  int p = y.n_cols;
  int qB = wB.n_cols;
  
  // ---------------------------------------------------------------------------
  int it, i, j, h;
  int ind = 0;
  
  for (it = 0; it < nrun; it++) {
    if(verbose && (it + 1) % 50 == 0) {
      Rcout << it + 1 << " : " << kstar << " active factors\n";
    }

    // -------------------------------------------------------------------------
    // Update eta
    eta = update_eta(Lambda, ps, y);
    
    // -------------------------------------------------------------------------
    // Update Sigma
    arma::mat Z_res = y - eta * trans(Lambda);
    for (j = 0; j < p; j++) {
      ps(j) = R::rgamma(a_sigma + 0.5 * n, 1 / (b_sigma + 0.5 * arma::accu(arma::pow(Z_res.col(j), 2))));
    }

    // -------------------------------------------------------------------------
    // Update GammaB with Metropolis-Hastings
    if( order_dependent == FALSE){
        pred = wB * Gamma;
        logit = arma::exp(pred) / (1 + arma::exp(pred));
        
        // Precompute L(l_{r,-h}) membership once (used in MH step)
        arma::mat in_L = compute_L_membership(Delta, p, k);
        
        // Update phi_L for MH proposal
        arma::mat Phi_L = arma::ones(p, k);
        arma::uvec Phi0 = arma::find(Phi == 0);
        arma::vec logit_phi0 = logit.elem(Phi0);
        arma::uvec which_zero = arma::randu(logit_phi0.n_elem) < (1 - logit_phi0) / (1 - logit_phi0 * p_constant);
        Phi_L.elem(Phi0.elem(arma::find(which_zero))) -= 1;
        
        // Metropolis-Hastings update for GammaB
        arma::mat Bh_1 = arma::diagmat(arma::ones(qB) / pow(sd_gammaB, 2));
        int accepted_this_iter = 0;
        for (h = 0; h < k; h++) {
          bool accepted = sample_gamma_MH(h, Gamma, Phi_L, Phi, in_L, Delta, wB, Bh_1, 
                                    p_constant, scale_factor_MH, p, cMH);
          if (accepted) accepted_this_iter++;
        }
        double acceptance_rate = static_cast<double>(accepted_this_iter) / k;
    } else {
        double acceptance_rate = 1;
    }
    //double acceptance_rate = 1;

    // -------------------------------------------------------------------------
    // Update Phi 
    if( order_dependent == FALSE){
      pred = wB * Gamma;
      logit = arma::exp(pred) / (1 + arma::exp(pred));
    }

    // loop over columns in random order
    arma::uvec h_order = arma::randperm(k);
    for (int idx = 0; idx < k; idx++) {
      int h = h_order(idx);
      // 1. Compute L membership for column h
      arma::vec in_L_h = compute_L_membership_single(h, Delta, p);

      // 2. Fill prior probability of phi 
      
       if(order_dependent == TRUE) {
    	    int m = p - k + 1;
          int l = 0; // counter for active entries seen so far
        
          for (int j = 0; j < p; j++) {
            if (in_L_h[j] == 0) {
                logit(j, h) = 1;
            } else {
                logit(j, h) = 1.0 / (p_constant*(m - l));
                l++;
            }
        }
        }

      // 3. Update potential pivots
      update_phi_potential_pivots_single(h, Phi, in_L_h, rho, logit, p_constant, 
                                 eta, Lambda_star, y, ps);
      // 4. Identify pivot for column h
      int l_h = identify_pivot(h, Phi, Delta, p);
      pivots(h) = l_h;  
      // 5. Update Delta for column h
      update_delta_column(h, l_h, Phi, Delta);
      // 6. Recompute membership after Delta update
      arma::vec in_L_h_updated = compute_L_membership_single(h, Delta, p);
      arma::vec not_in_L_h = arma::ones(p) - in_L_h_updated;
      // 7. Update remaining positions
      update_phi_remaining_single(h, Phi, not_in_L_h, l_h, rho, logit, p_constant, 
                         eta, Lambda_star, y, ps, Delta, cMH, order_dependent);
   }


    // -------------------------------------------------------------------------
    // Update Lambda_star and Lambda
    arma::mat etarho = trans(eta.each_row() % trans(rho));
    for (j = 0; j < p; j++) {
      Lambda_star.row(j) = update_Lambda_star(j, etarho, Delta, Plam, ps, y);
    }
    Lambda = (Lambda_star.each_row() % trans(sqrt(rho))) % Delta;
    
    // -------------------------------------------------------------------------
    // Update d 
    for (h = 0; h < k; h++) {
      d(h) = update_d(h, Delta, rho, eta, Lambda_star, y, ps, w);
    }
    rho = arma::ones(k);
    rho.elem(find(d <= arma::linspace<arma::vec>(0, k - 1, k))) -= 1;
    
    // Update Plam
    arma::vec Plam_diag(k);
    for (h = 0; h < k; h++) {
      Plam_diag(h) = R::rgamma(a_theta + 0.5 * p, 1 / (b_theta + 0.5 * arma::accu(arma::pow(Lambda_star.col(h), 2))));
    }
    Plam = arma::diagmat(Plam_diag);
    
    // Update stick-breaking weights
    for (h = 0; h < k - 1; h++) {
      v(h) = R::rbeta(1 + arma::accu(d == h), alpha + arma::accu(d > h));
    }
    v(k - 1) = 1;
    w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
    
    // -------------------------------------------------------------------------
    // Save sampled values (after burn-in period)
    if ((it + 1) % thin == 0 && (it + 1) > burn) {
      if(out.containsElementNamed("gamma")) { GAMMA[ind] = Gamma; }
      if(out.containsElementNamed("eta")) { ETA[ind] = eta; }
      if(out.containsElementNamed("lambda")) { LAMBDA[ind] = Lambda; }
      if(out.containsElementNamed("preccol")) { SIG[ind] = ps; }
      if(out.containsElementNamed("numFactors")) { K[ind] = kstar; }
      if(out.containsElementNamed("activeFactors")) { RHO[ind] = rho; }
      DELTA[ind] = Delta;
      ACC_RATE[ind] = acceptance_rate;
      ind += 1;
    }
    
    // -------------------------------------------------------------------------
    // Adaptation (adjust truncation level)
    if (uu(it) < prob(it) && (it + 1) > start_adapt) {
      arma::uvec active = find(d > arma::linspace<arma::vec>(0, k - 1, k));
      int kstar_new = active.n_elem;
      kstar = kstar_new;
      
      if (kstar < k - 1) {
        // Decrease truncation to kstar
        k = kstar + 1;
        eta = arma::join_rows(eta.cols(active), rnorm_vec(n, 0, 1));
        double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
        Plam_diag = join_elem(Plam_diag.elem(active), vartheta_k);
        Plam = arma::diagmat(Plam_diag);
        Lambda_star = arma::join_rows(Lambda_star.cols(active), rnorm_vec(p, 0, sqrt(vartheta_k)));
        Phi = arma::join_rows(Phi.cols(active), rbinom_vec(p, 1, p_constant));
        // Recompute Delta and pivots from Phi deterministically
        pivots.resize(k); // Resize pivots vector to new k
        Delta = arma::zeros(p, k);
          for (int h_adapt = 0; h_adapt < k; h_adapt++) {
            int l_h_adapt = identify_pivot(h_adapt, Phi, Delta, p);
            update_delta_column(h_adapt, l_h_adapt, Phi, Delta);
            pivots(h_adapt) = l_h_adapt;
            }
        rho = join_elem(rho.elem(active), 1);
        Lambda = arma::join_rows(Lambda.cols(active), Lambda_star.col(k - 1) % Delta.col(k - 1));
        Gamma = arma::join_rows(Gamma.cols(active), rnorm_vec(qB, 0, sqrt(sd_gammaB)));
        w = join_elem(w.elem(active), 1 - sum(w.elem(active)));
        v = join_elem(v.elem(active), 1);
        d = join_elem(d.elem(active), k - 1);
      } else if (k < kmax) {
        // Increase truncation by 1
        k += 1;
        eta = arma::join_rows(eta, rnorm_vec(n, 0, 1));
        double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
        Plam_diag = join_elem(Plam_diag, vartheta_k);
        Plam = arma::diagmat(Plam_diag);
        Lambda_star = arma::join_rows(Lambda_star, rnorm_vec(p, 0, sqrt(vartheta_k)));
        Phi = arma::join_rows(Phi, rbinom_vec(p, 1, p_constant));
        // Delta = arma::join_rows(Delta, rbinom_vec(p, 1, p_constant));
        //pivots = join_elem(pivots, -1);  // Update pivots vector
        // Recompute Delta and pivots from Phi deterministically
        pivots.resize(k); // Resize pivots vector to new k
        Delta = arma::zeros(p, k);
          for (int h_adapt = 0; h_adapt < k; h_adapt++) {
            int l_h_adapt = identify_pivot(h_adapt, Phi, Delta, p);
            update_delta_column(h_adapt, l_h_adapt, Phi, Delta);
            pivots(h_adapt) = l_h_adapt;
          }
        rho = join_elem(rho, 1);
        Lambda = arma::join_rows(Lambda, Lambda_star.col(k - 1) % Delta.col(k - 1));
        Gamma = arma::join_rows(Gamma, rnorm_vec(qB, 0, sqrt(sd_gammaB)));
        v(k - 2) = R::rbeta(1, alpha);
        v = join_elem(v, 1);
        w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
        d = join_elem(d, k - 1);
      }
    }
  }
 
  // ---------------------------------------------------------------------------
  // Return results
  if(out.containsElementNamed("gamma")) { out["gamma"] = GAMMA; }
  if(out.containsElementNamed("eta")) { out["eta"] = ETA; }
  if(out.containsElementNamed("lambda")) { out["lambda"] = LAMBDA; }
  if(out.containsElementNamed("preccol")) { out["preccol"] = SIG; }
  if(out.containsElementNamed("numFactors")) { out["numFactors"] = K; }
  if(out.containsElementNamed("activeFactors")) { out["activeFactors"] = RHO; }
  out["delta"] = DELTA;
  out["accRate"] = ACC_RATE;
  
  return out;

}
