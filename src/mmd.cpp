// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Compute the unbiased MMD^2 statistic
//'
//' Unbiased estimator of \eqn{MMD^2}:
//' \deqn{\frac{1}{n(n-1)} \sum_{i \neq j} K_{xx}(i,j) + \frac{1}{m(m-1)} \sum_{i \neq j} K_{yy}(i,j) - \frac{2}{nm} \sum_{i,j} K_{xy}(i,j)}
//'
//' @param Kxx n x n kernel matrix for sample X.
//' @param Kyy m x m kernel matrix for sample Y.
//' @param Kxy n x m kernel matrix between X and Y.
//' @return Scalar unbiased MMD^2 value.
//' @keywords internal
// [[Rcpp::export]]
double mmd2_unbiased_cpp(const arma::mat& Kxx,
                          const arma::mat& Kyy,
                          const arma::mat& Kxy) {
  int n = Kxx.n_rows;
  int m = Kyy.n_rows;

  double sum_xx = arma::accu(Kxx) - arma::trace(Kxx);
  double sum_yy = arma::accu(Kyy) - arma::trace(Kyy);
  double sum_xy = arma::accu(Kxy);

  double mmd2 = sum_xx / ((double)n * (n - 1)) +
                sum_yy / ((double)m * (m - 1)) -
                2.0 * sum_xy / ((double)n * m);

  return mmd2;
}

//' Compute the biased MMD^2 statistic
//'
//' @param Kxx n x n kernel matrix for sample X.
//' @param Kyy m x m kernel matrix for sample Y.
//' @param Kxy n x m kernel matrix between X and Y.
//' @return Scalar biased MMD^2 value.
//' @keywords internal
// [[Rcpp::export]]
double mmd2_biased_cpp(const arma::mat& Kxx,
                        const arma::mat& Kyy,
                        const arma::mat& Kxy) {
  int n = Kxx.n_rows;
  int m = Kyy.n_rows;

  double mmd2 = arma::accu(Kxx) / ((double)n * n) +
                arma::accu(Kyy) / ((double)m * m) -
                2.0 * arma::accu(Kxy) / ((double)n * m);

  return mmd2;
}
