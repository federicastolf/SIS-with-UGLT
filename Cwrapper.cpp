#include <RcppArmadillo.h>
#include "rcpp_pgdraw.cpp"       /* samplepg ('pgdraw' R package) */
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

/* -------------------------------------------------------------------------- */
// Functions in this file:
//  - truncnorm_lg
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

// Update the hth column of GammaB in the Adaptive Gibbs Sampler
//
// @param h An integer number.
// @param wB A pxqB matrix.
// @param Dt A pxk matrix.
// @param Bh_1 A qBxqB matrix.
// @param Phi_L A pxk matrix.
//
// @return A qBx1 matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_gamma(const int& h, const arma::mat& wB, const arma::mat& Dt,
                        const arma::mat& Bh_1, const arma::mat& Phi_L) {
  int qB = wB.n_cols;
  arma::mat Qbeta = trans(wB.each_col() % Dt.col(h)) * wB + Bh_1;
  arma::mat Lbeta = trans(arma::chol(Qbeta));
  arma::vec bbeta = trans(wB) * (Phi_L.col(h) - 0.5);
  // mean
  arma::mat vbeta = solve(trimatl(Lbeta), bbeta);
  arma::mat mbeta = solve(trimatu(trans(Lbeta)), vbeta);
  // var
  arma::vec zbeta = rnorm_vec(qB, 0, 1);
  arma::mat ybeta = solve(trimatu(trans(Lbeta)), zbeta);
  arma::mat GammaBh = ybeta + mbeta;
  return GammaBh;
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

// Update Phi in the Adaptive Gibbs Sampler
//
// @param rho A k-dimensional vector.
// @param logit A pxk matrix.
// @param p_constant A number in (0,1).
// @param p An integer.
// @param n An integer.
// @param eta A nxk matrix.
// @param lambdastar A pxk matrix.
// @param Phi A pxk matrix.
// @param Z A nxp matrix.
// @param ps A p-dimensional vector.
// @param k An integer.
//
// @return A pxk matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_Phi(arma::mat& Phi, const arma::vec& rho, const arma::mat& logit,
                     const double& p_constant, const arma::mat& eta,
                     const arma::mat& lambdastar, const arma::mat& Z, const arma::vec& ps) {
  // define variables
  int k = Phi.n_cols;
  int n = eta.n_rows;
  int p = Phi.n_rows;
  int i, j, l, h;
  arma::uword f;
  arma::vec p_phi0(p), p_phi1(p), p_phi_sum(p), lnorm0(p), lnorm1(p);
  arma::uvec wr0 = find(rho == 0);
  arma::uvec wr1 = find(rho == 1);
  arma::vec sdy = sqrt(1 / ps);
  // update for inactive factors
  for (f = 0; f < wr0.n_elem; f++) {
    h = wr0(f);
    p_phi0 = 1 - logit.col(h) * p_constant;
    p_phi1 = logit.col(h) * p_constant;
    p_phi_sum = p_phi0 + p_phi1;
    for (j = 0; j < p; j++) {
      Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
    }
  }
  // update for active factors
  for (f = 0; f < wr1.n_elem; f++) {
    h = wr1(f);
    lnorm0 = arma::zeros(p);
    lnorm1 = arma::zeros(p);
    for (j = 0; j < p; j++) {
      for (i = 0; i < n; i++) {
        // initialize mean
        double muijh = -1 * rho(h) * Phi(j, h) * lambdastar(j, h) * eta(i, h);
        for (l = 0; l < k; l++) {
            muijh += rho(l) * Phi(j, l) * lambdastar(j, l) * eta(i, l);
        }
        lnorm0(j) += R::dnorm(Z(i, j), muijh, sdy(j), 1);
        // update mean
        muijh += rho(h) * lambdastar(j, h) * eta(i, h);
        lnorm1(j) += R::dnorm(Z(i, j), muijh, sdy(j), 1);
      }
      // adjust the scale
      double mlnorm = std::max(lnorm0(j), lnorm1(j));
      lnorm0(j) -= mlnorm;
      lnorm1(j) -= mlnorm;
    }
    p_phi0 = exp(lnorm0 + log(1 - logit.col(h) * p_constant));
    p_phi1 = exp(lnorm1 + log(logit.col(h) * p_constant));
    p_phi_sum = p_phi0 + p_phi1;
    for (j = 0; j < p; j++) {
      Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
    }
  }
  return Phi;
}

/* -------------------------------------------------------------------------- */



// Adaptive Gibbs Sampler (AGS)
// Implementation in C++ of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor model with Structured Increasing Shrinkage (SIS) prior.
//
// [[Rcpp::export]]
Rcpp::List Rcpp_cosin(double alpha, double a_sigma, double b_sigma, double a_theta, double b_theta,
  double sd_gammaB, double p_constant,
  arma::mat y, 
  arma::mat wB,
  int burn, int nrun, int thin, int start_adapt, int kmax,
  arma::mat eta, arma::mat Gamma,  arma::mat Lambda, arma::mat Lambda_star, arma::vec d, int kstar,
  arma::mat logit, arma::vec rho, arma::mat Phi, arma::mat Plam, arma::mat pred, arma::vec ps,
  arma::vec v,   arma::vec w,
  Rcpp::List out, bool verbose,
  arma::vec uu, arma::vec prob, int sp,
  arma::vec lpiv, arma::mat Delta) {
  // ---------------------------------------------------------------------------
  // output
  Rcpp::List GAMMA(sp);
  Rcpp::List ETA(sp);
  Rcpp::List LAMBDA(sp);
  Rcpp::List SIG(sp);
  arma::vec K(sp);
  // ---------------------------------------------------------------------------
  // matrix dimensions
  int k = Lambda_star.n_cols;
  int n = y.n_rows;
  int p = y.n_cols;
  int qB = wB.n_cols; //meta-covariate dim
  // ---------------------------------------------------------------------------
  int it, i, j, h;
  int ind = 0;
  for (it = 0; it < nrun; it++) {
    if(verbose && (it + 1) % 50 == 0) {
      Rcout << it + 1 << " : " << k << " active factors\n";
    }

    // update eta
    eta = update_eta(Lambda, ps, y);
    // -------------------------------------------------------------------------
    // update Sigma
    arma::mat Z_res = y - eta * trans(Lambda);
    for (j = 0; j < p; j++) {
      ps(j) = R::rgamma(a_sigma + 0.5 * n, 1 / (b_sigma + 0.5 * arma::accu(arma::pow(Z_res.col(j), 2))));
    }
    // -------------------------------------------------------------------------
    // DA MODIFICARE
    // 6 - update GammaB
    pred = wB * Gamma;
    logit = arma::exp(pred) / (1 + arma::exp(pred));
    // 6.1 Update phi_L
    arma::mat Phi_L = arma::ones(p, k);
    arma::uvec Phi0 = arma::find(Phi == 0);
    arma::vec logit_phi0 = logit.elem(Phi0);
    arma::uvec which_zero = arma::randu(logit_phi0.n_elem) < (1 - logit_phi0) / (1 - logit_phi0 * p_constant);
    Phi_L.elem(Phi0.elem(arma::find(which_zero))) -= 1;
    // 6.2 Polya gamma
    arma::mat Dt(p, k);
    for (j = 0; j < k; j++) {
      for (i = 0; i < p; i++) {
        Dt(i, j) = samplepg(pred(i, j));
      }
    }
    // 6.3 Update gammaB
    arma::mat Bh_1 = arma::diagmat(arma::ones(qB) / pow(sd_gammaB, 2));
    for (h = 0; h < k; h++) {
      Gamma.col(h) = update_gamma(h, wB, Dt, Bh_1, Phi_L);
    }
    // -------------------------------------------------------------------------
    // 7 - update Lambda_star and Lambda -> SAME
    arma::mat etarho = trans(eta.each_row() % trans(rho));
    for (j = 0; j < p; j++) {
      Lambda_star.row(j) = update_Lambda_star(j, etarho, Phi, Plam, ps, y);
    }
    Lambda = (Lambda_star.each_row() % trans(sqrt(rho))) % Phi;
    // -------------------------------------------------------------------------
    // 8.1 - update d-> SAME
    for (h = 0; h < k; h++) {
      d(h) = update_d(h, Phi, rho, eta, Lambda_star, y, ps, w);
    }
    rho = arma::ones(k);
    rho.elem(find(d <= arma::linspace<arma::vec>(0, k - 1, k))) -= 1;
    // 8.2
    arma::vec Plam_diag(k);
    for (h = 0; h < k; h++) {
      Plam_diag(h) = R::rgamma(a_theta + 0.5 * p, 1 / (b_theta + 0.5 * arma::accu(arma::pow(Lambda_star.col(h), 2))));
    }
    Plam = arma::diagmat(Plam_diag);
    // 8.3
    for (h = 0; h < k - 1; h++) {
      v(h) = R::rbeta(1 + arma::accu(d == h), alpha + arma::accu(d > h));
    }
    v(k - 1) = 1;
    w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
    // -------------------------------------------------------------------------
    // 9 - update Phi --> da anticipare dopo step 6
    pred = wB * Gamma;
    logit = arma::exp(pred) / (1 + arma::exp(pred));
    Phi = update_Phi(Phi, rho, logit, p_constant, eta, Lambda_star, y, ps);
    // -------------------------------------------------------------------------
    // save sampled values (after burn-in period)
    if ((it + 1) % thin == 0 && (it + 1) > burn) {
      if(out.containsElementNamed("gamma")) { GAMMA[ind] = Gamma; }
      if(out.containsElementNamed("eta")) { ETA[ind] = eta; }
      if(out.containsElementNamed("lambda")) { LAMBDA[ind] = Lambda; }
      if(out.containsElementNamed("sigmacol")) { SIG[ind] = ps; }
      if(out.containsElementNamed("numFactors")) { K[ind] = kstar; }
      ind += 1;
    }
    // -------------------------------------------------------------------------
    // Adaptation
    if (uu(it) < prob(it) && (it + 1) > start_adapt) {
      arma::uvec active = find(d > arma::linspace<arma::vec>(0, k - 1, k));
      int kstar_new = active.n_elem;
      kstar = kstar_new;
      if (kstar < k - 1) {
        // set truncation to kstar and subset all variables, keeping only active columns
        k = kstar + 1;
        eta = arma::join_rows(eta.cols(active), rnorm_vec(n, 0, 1));
        double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
        Plam_diag = join_elem(Plam_diag.elem(active), vartheta_k);
        Plam = arma::diagmat(Plam_diag);
        Lambda_star = arma::join_rows(Lambda_star.cols(active), rnorm_vec(p, 0, sqrt(vartheta_k)));
        Phi = arma::join_rows(Phi.cols(active), rbinom_vec(p, 1, p_constant));
        rho = join_elem(rho.elem(active), 1);
        Lambda = arma::join_rows(Lambda.cols(active), Lambda_star.col(k - 1) % Phi.col(k - 1));
        Gamma = arma::join_rows(Gamma.cols(active), rnorm_vec(qB, 0, sqrt(sd_gammaB)));
        w = join_elem(w.elem(active), 1 - sum(w.elem(active)));
        v = join_elem(v.elem(active), 1);
        d = join_elem(d.elem(active), k - 1);
      } else if (k < kmax) {
      // increase truncation by 1 and extend all variables, sampling from the prior/model
      k += 1;
      eta = arma::join_rows(eta, rnorm_vec(n, 0, 1));
      double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
      Plam_diag = join_elem(Plam_diag, vartheta_k);
      Plam = arma::diagmat(Plam_diag);
      Lambda_star = arma::join_rows(Lambda_star, rnorm_vec(p, 0, sqrt(vartheta_k)));
      Phi = arma::join_rows(Phi, rbinom_vec(p, 1, p_constant));
      rho = join_elem(rho, 1);
      Lambda = arma::join_rows(Lambda, Lambda_star.col(k - 1) % Phi.col(k - 1));
      Gamma = arma::join_rows(Gamma, rnorm_vec(qB, 0, sqrt(sd_gammaB)));
      v(k - 2) = R::rbeta(1, alpha);
      v = join_elem(v, 1);
      w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
      d = join_elem(d, k - 1);
      }
    }
    // -------------------------------------------------------------------------
  }
 
  // ---------------------------------------------------------------------------
  if(out.containsElementNamed("gamma")) { out["gamma"] = GAMMA; }
  if(out.containsElementNamed("eta")) { out["eta"] = ETA; }
  if(out.containsElementNamed("lambda")) { out["lambda"] = LAMBDA; }
  if(out.containsElementNamed("sigmacol")) { out["sigmacol"] = SIG; }
  if(out.containsElementNamed("numFactors")) { out["numFactors"] = K; }
  // ---------------------------------------------------------------------------
  return out;
  // ---------------------------------------------------------------------------
}
