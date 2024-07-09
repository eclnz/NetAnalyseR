#include <Rcpp.h>
#include <algorithm> // For std::swap
#include <random>    // For std::mt19937, std::random_device

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;


//' Rewire a Network Matrix
//'
//' This function takes a network represented as an adjacency matrix and
//' performs a specified number of edge rewiring iterations to randomize the network.
//' It checks for network density and skips rewiring if the network is fully connected.
//' @name rewireNetworkCpp
//' @param R A numeric matrix representing the network's adjacency matrix.
//' @param initialIter int number of iterations for the rewiring process.
//' @param validate bool specifying whether to check the validity of the input matrix
//'
//' @return A rewired numeric matrix representing the network's adjacency matrix.
//' @export
// [[Rcpp::export]]
NumericMatrix rewireNetworkCpp(NumericMatrix R, int initialIter = 2, bool validate = true) {
    int n = R.nrow(); // Size of the matrix
    std::vector<int> i, j;
    int K = 0; // Count of edges
    
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    // Count edges in the lower triangular part of 'R'
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < row; ++col) {
            if (R(row, col) > 0) {
                i.push_back(row);
                j.push_back(col);
                ++K;
            }
        }
    }
    
    // Edge cases handling
    if (K == 0) {
        // Network has no connections. Rewiring is not applicable.
        return R; // Return the original matrix unmodified
    }
    
    if (K == 1) {
        // Network has only one connection. Rewiring would not alter the structure.
        return R; // Return the original matrix unmodified
    }
    
    // Calculate network density
    if(validate){
        double density = 2 * K / static_cast<double>(n * (n - 1));
        if (density == 1) {
            // Network density is 1. Rewiring is not possible.
            return R; // Return the original matrix unmodified
        }
    }
    
    int ITER = K * initialIter;
    int maxAttempts = std::round(static_cast<double>(n) * K / (n * (n - 1)));
    
    for (int iter = 0; iter < ITER; ++iter) {
        bool rewired = false;
        for (int att = 0; att <= maxAttempts && !rewired; ++att) {
            int e1, e2, a, b, c, d;
            do {
                e1 = g() % K;
                e2 = g() % K;
            } while (e1 == e2);
            
            a = i[e1]; b = j[e1];
            c = i[e2]; d = j[e2];
            
            // Ensure all vertices are different and no self-loops
            if (a != c && a != d && b != c && b != d) {
                // Flip edge c-d with 50% probability
                if (dis(g) > 0.5) {
                    std::swap(i[e2], j[e2]);
                    c = i[e2]; d = j[e2];
                }
                
                // Check rewiring condition
                if (R(a, d) == 0 && R(c, b) == 0) {
                    R(a, d) = R(a, b); R(a, b) = 0;
                    R(d, a) = R(b, a); R(b, a) = 0;
                    R(c, b) = R(c, d); R(c, d) = 0;
                    R(b, c) = R(d, c); R(d, c) = 0;
                    
                    // Update edge indices to reflect rewiring
                    j[e1] = d;
                    j[e2] = b;
                    rewired = true; // Mark successful rewiring
                }
            }
        }
        // If no rewiring occurred in this iteration, break out of the loop
    }
    
    return R; // Return the rewired matrix
}

//' Generate a Series of Rewired Matrices
//'
//' Rewires a matrix n times, performing 10 iterations each time.
//' Saves all random matrices generated as a list.
//' The output is in the format of a list. If manually using the output in normalised mesures
//' the output must be converted to an array with `abind::abind(matrix, along = 3)`
//' @name generateRewiredMatrices
//' @param initialMatrix A numeric matrix representing the initial network's adjacency matrix.
//' @param n The number of rewired matrices to generate.
//'
//' @return A list of 'n' rewired matrices.
//' @export
// [[Rcpp::export]]
std::vector<NumericMatrix> generateRewiredMatrices(NumericMatrix initialMatrix, int n = 100) {
    // Preallocate the vector to hold the matrices
    std::vector<NumericMatrix> matrices;
    matrices.reserve(n);
    
    // Clone the initial matrix to avoid modifying the original
    NumericMatrix matrix = clone(initialMatrix);

    // Generate the remaining n-1 networks
    for (int i = 0; i < n; ++i) {
        // Perform 10 rewiring iterations for each network
        matrix = rewireNetworkCpp(matrix, 2, false); // Perform 2 rewirings in each iteration
        matrices.push_back(matrix); // Save the matrix after every 10 rewirings
    }
    
    return matrices;
}
