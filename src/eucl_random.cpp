
#include <Rcpp.h>
#include <algorithm> // For std::swap
#include <random>    // For std::mt19937, std::random_device
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List rewireTwoNetworksCpp(NumericMatrix originalMatrix, NumericMatrix secondaryMatrix, int initialIter) {
    int n = originalMatrix.nrow(); // Assuming both matrices are square and of the same dimension
    std::vector<int> i, j;
    int K = 0; // Count of edges in the original matrix
    
    // Clone the secondaryMatrix to avoid altering the original
    NumericMatrix secondaryMatrixClone = clone(secondaryMatrix);
    
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_real_distribution<> dis(0, 1);
    
    // Identifying edges in the lower triangular part of the original matrix
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < row; ++col) {
            if (originalMatrix(row, col) > 0) {
                i.push_back(row);
                j.push_back(col);
                ++K;
            }
        }
    }
    
    int ITER = K * initialIter;
    
    for (int iter = 0; iter < ITER; ++iter) {
        int e1, e2, a, b, c, d;
        do {
            e1 = g() % K;
            e2 = g() % K;
        } while (e1 == e2);
        
        a = i[e1]; b = j[e1];
        c = i[e2]; d = j[e2];
        
        // Ensuring all vertices are different
        if (a != c && a != d && b != c && b != d) {
            // Flip edge c-d with 50% probability
            if (dis(g) > 0.5) {
                std::swap(i[e2], j[e2]);
                c = i[e2]; d = j[e2];
            }
            
            // Rewiring condition based on the original matrix
            if (originalMatrix(a, d) == 0 && originalMatrix(c, b) == 0) {
                // Apply the rewiring to both matrices
                originalMatrix(a, d) = originalMatrix(a, b); originalMatrix(a, b) = 0;
                originalMatrix(d, a) = originalMatrix(b, a); originalMatrix(b, a) = 0;
                originalMatrix(c, b) = originalMatrix(c, d); originalMatrix(c, d) = 0;
                originalMatrix(b, c) = originalMatrix(d, c); originalMatrix(d, c) = 0;
                
                secondaryMatrixClone(a, d) = secondaryMatrixClone(a, b); secondaryMatrixClone(a, b) = 0;
                secondaryMatrixClone(d, a) = secondaryMatrixClone(b, a); secondaryMatrixClone(b, a) = 0;
                secondaryMatrixClone(c, b) = secondaryMatrixClone(c, d); secondaryMatrixClone(c, d) = 0;
                secondaryMatrixClone(b, c) = secondaryMatrixClone(d, c); secondaryMatrixClone(d, c) = 0;
            }
        }
    }
    
    // Return both matrices as a List
    return List::create(Named("originalMatrix") = originalMatrix, Named("secondaryMatrixClone") = secondaryMatrixClone);
}

