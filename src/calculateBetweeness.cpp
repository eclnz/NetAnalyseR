#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector calculateBetweennessCentrality(const NumericMatrix& D, const NumericMatrix& NP) {
    int n = D.nrow();
    NumericVector BC(n, 0.0);

    // Calculate betweenness centrality
    for (int s = 0; s < n; s++) {
        for (int t = 0; t < n; t++) {
            if (s != t) {
                for (int v = 0; v < n; v++) {
                    if (s != v && t != v) {
                        if (D(s, v) + D(v, t) == D(s, t)) {  // Check if v is on the shortest path from s to t
                            BC[v] += (NP(s, v) * NP(v, t)) / NP(s, t);
                        }
                    }
                }
            }
        }
    }

    return BC;
}
