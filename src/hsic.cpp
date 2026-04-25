// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

//' Compute the HSIC statistic (unweighted)
//'
//' Biased HSIC estimator: (1/n^2) * trace(K_x_c %*% K_y_c)
//' where K_c = H K H with H = I - 1/n.
//'
//' @param Kx n x n kernel matrix for X.
//' @param Ky n x n kernel matrix for Y.
//' @return Scalar HSIC value.
//' @keywords internal
// [[Rcpp::export]]
double hsic_stat_cpp(const arma::mat& Kx, const arma::mat& Ky) {
  int n = Kx.n_rows;
  double n2 = (double)(n * n);

  // Centre kernel matrices: Kc = H K H where H = I - 1/n
  arma::vec ones_n = arma::ones<arma::vec>(n);
  arma::vec Kx_row_mean = Kx * ones_n / n;
  arma::vec Ky_row_mean = Ky * ones_n / n;
  double Kx_mean = arma::accu(Kx) / n2;
  double Ky_mean = arma::accu(Ky) / n2;

  // HSIC = (1/n^2) * sum_{i,j} Kxc_{ij} * Kyc_{ij}
  // Kxc_{ij} = Kx_{ij} - Kx_row_mean_i - Kx_row_mean_j + Kx_mean
  double hsic = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double kxc = Kx(i, j) - Kx_row_mean(i) - Kx_row_mean(j) + Kx_mean;
      double kyc = Ky(i, j) - Ky_row_mean(i) - Ky_row_mean(j) + Ky_mean;
      hsic += kxc * kyc;
    }
  }

  return hsic / n2;
}

//' Compute the weighted HSIC statistic (for bd-HSIC)
//'
//' Weighted version: \eqn{\sum_{i,j} w_i \, w_j \, K_{xc,ij} \, K_{yc,ij}}
//'
//' @param Kx n x n kernel matrix for X.
//' @param Ky n x n kernel matrix for Y.
//' @param w Weight vector of length n.
//' @return Scalar weighted HSIC value.
//' @keywords internal
// [[Rcpp::export]]
double weighted_hsic_stat_cpp(const arma::mat& Kx,
                               const arma::mat& Ky,
                               const arma::vec& w) {
  int n = Kx.n_rows;

  // Weighted centring
  double W = arma::accu(w);
  arma::vec wKx_row = (Kx * w) / W;
  arma::vec wKy_row = (Ky * w) / W;
  double wKx_mean = arma::dot(w, Kx * w) / (W * W);
  double wKy_mean = arma::dot(w, Ky * w) / (W * W);

  double hsic = 0.0;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      double kxc = Kx(i, j) - wKx_row(i) - wKx_row(j) + wKx_mean;
      double kyc = Ky(i, j) - wKy_row(i) - wKy_row(j) + wKy_mean;
      hsic += w(i) * w(j) * kxc * kyc;
    }
  }

  return hsic / (W * W);
}
