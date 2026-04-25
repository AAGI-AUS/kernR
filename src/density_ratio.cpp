// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' RuLSIF kernel density ratio estimation (core solver)
//'
//' Solves the RuLSIF optimisation \eqn{\theta = (H + \lambda I)^{-1} h}
//' with non-negativity constraint (clamp negatives to 0).
//'
//' @param H Gram matrix (n_basis x n_basis).
//' @param h Mean kernel vector (n_basis).
//' @param lambda Regularisation parameter.
//' @return Coefficient vector theta.
//' @keywords internal
// [[Rcpp::export]]
arma::vec rulsif_solve_cpp(const arma::mat& H,
                            const arma::vec& h,
                            double lambda) {
  int p = H.n_rows;
  arma::mat reg = H + lambda * arma::eye(p, p);
  arma::vec theta = arma::solve(reg, h);
  // Non-negativity
  theta.clamp(0.0, arma::datum::inf);
  return theta;
}
