// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

//' Nystrom low-rank factorisation
//'
//' Given the cross-block C (n x m) and inner block W (m x m) of a kernel
//' matrix evaluated at m landmark points, returns U such that the rank-r
//' Nystrom approximation is \eqn{\tilde K = U U^\top}, with
//' \eqn{r \le m} columns retained after eigenvalue thresholding.
//'
//' Uses the symmetric eigendecomposition \eqn{W = V \Lambda V^\top} and
//' sets \eqn{U = C V \Lambda_+^{-1/2}} where \eqn{\Lambda_+} keeps only
//' eigenvalues above \eqn{\mathrm{tol} \cdot \max(\Lambda)}. This handles
//' rank-deficient W (e.g. when landmarks are near-duplicates) without
//' the numerical instability of a direct pseudo-inverse.
//'
//' @param C Cross-block kernel matrix (n x m).
//' @param W Landmark-landmark kernel matrix (m x m), assumed symmetric PSD.
//' @param tol Relative eigenvalue threshold. Eigenvalues
//'   \eqn{\lambda_i \le \mathrm{tol} \cdot \max(\Lambda)} are dropped.
//'   Default 1e-10.
//' @return U: an n x r matrix with \eqn{r \le m} columns. Approximation
//'   \eqn{\tilde K = U U^\top} is symmetric PSD by construction.
//' @keywords internal
// [[Rcpp::export]]
arma::mat nystrom_factor_cpp(const arma::mat& C,
                             const arma::mat& W,
                             double tol = 1e-10) {
  arma::mat W_sym = 0.5 * (W + W.t());

  arma::vec eigval;
  arma::mat eigvec;
  bool ok = arma::eig_sym(eigval, eigvec, W_sym);
  if (!ok) {
    Rcpp::stop("Nystrom: eigendecomposition of W failed.");
  }

  double lam_max = eigval.max();
  if (lam_max <= 0.0) {
    return arma::mat(C.n_rows, 0);
  }

  arma::uvec keep = arma::find(eigval > tol * lam_max);
  if (keep.is_empty()) {
    return arma::mat(C.n_rows, 0);
  }

  arma::vec lam_keep = eigval.elem(keep);
  arma::mat V_keep = eigvec.cols(keep);
  arma::vec inv_sqrt = 1.0 / arma::sqrt(lam_keep);

  arma::mat U = C * V_keep;
  U.each_row() %= inv_sqrt.t();
  return U;
}
