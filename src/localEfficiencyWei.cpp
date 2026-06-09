//' @export

#include <RcppArmadillo.h>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

arma::mat floydWarshallRcpp(const arma::mat& graph);

// [[Rcpp::export]]
arma::vec localEfficiencyWei(const arma::mat& W_original) {
  arma::mat W = W_original;
  W.diag().zeros();
  int n = W.n_rows;

  arma::vec E(n, arma::fill::zeros);
  arma::mat A = arma::conv_to<arma::mat>::from(W > 0);

  // Precompute cube-root weights: W^(1/3)
  arma::mat cubeW = arma::pow(W, 1.0 / 3.0);

  // Precompute cube-root lengths: (1/w)^(1/3) for non-zero w, else 0
  arma::mat cubeL(n, n, arma::fill::zeros);
  for (int i = 0; i < n; i++)
    for (int j = 0; j < n; j++)
      if (W(i, j) > 0) cubeL(i, j) = std::pow(1.0 / W(i, j), 1.0 / 3.0);

  for (int node = 0; node < n; node++) {
    arma::uvec nbrs = arma::find(W.row(node) > 0);
    if (nbrs.n_elem <= 1) continue;

    // Symmetrised cube-root weight vector (BCT formula: sw = W(i,V).^(1/3) + W(V,i).^(1/3))
    arma::vec sw = cubeW.row(node).cols(nbrs).t() + cubeW.col(node).rows(nbrs);

    // Subgraph cube-root lengths and shortest paths
    arma::mat D = floydWarshallRcpp(cubeL(nbrs, nbrs));
    arma::mat e = 1.0 / D;
    e.diag().zeros();
    arma::mat se = e + e.t();

    // Numerator: sum((sw * sw^T) .* se) / 2
    double numer = arma::accu((sw * sw.t()) % se) / 2.0;

    if (numer > 0.0) {
      arma::vec sa = A.row(node).cols(nbrs).t() + A.col(node).rows(nbrs);
      double denom = std::pow(arma::accu(sa), 2.0) - arma::dot(sa, sa);
      E(node) = numer / denom;
    }
  }

  return E;
}
