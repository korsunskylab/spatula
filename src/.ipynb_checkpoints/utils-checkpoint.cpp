#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <vector>
using namespace std; 
using namespace Rcpp; 

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
IntegerMatrix triplets_to_pairs(const IntegerMatrix& triplets) {
    unsigned N = triplets.nrow(); 
    IntegerMatrix res(6*N, 2);
    for (int i = 0; i < N; i++) {
        res(6*i, 0) = triplets(i, 0);
        res(6*i, 1) = triplets(i, 1);
        
        res(6*i+1, 0) = triplets(i, 0);
        res(6*i+1, 1) = triplets(i, 2);
        
        res(6*i+2, 0) = triplets(i, 1);
        res(6*i+2, 1) = triplets(i, 2);   
        
        res(6*i+3, 0) = triplets(i, 1);
        res(6*i+3, 1) = triplets(i, 0);
        
        res(6*i+4, 0) = triplets(i, 2);
        res(6*i+4, 1) = triplets(i, 0);
        
        res(6*i+5, 0) = triplets(i, 2);
        res(6*i+5, 1) = triplets(i, 1);   
    }
    return res;
}

// [[Rcpp::export]]
IntegerMatrix dgCMat_to_pairs_cpp(const IntegerVector& ivec, const IntegerVector& pvec, const IntegerVector&xvec) {
    int i = 0;
    IntegerMatrix pairs(ivec.size(), 3); 
    for (int p = 0; p < pvec.size(); p++) {
//         cout << p_cusum << "," << i << endl;
        while(i < pvec[p]) { 
            // 1-indexing for R 
            pairs(i, 0) = p + 1; 
            pairs(i, 1) = ivec(i) + 1;
            pairs(i, 2) = xvec(i);
            i++; 
        }
    }
    return pairs;     
}




// [[Rcpp::export]]
std::vector<arma::umat> get_idx_cpp(const arma::umat& X, bool only_boundary=true) {
    size_t ncells = X.max(); 
    std::vector<arma::umat> res(ncells);     
    arma::umat X_bound = arma::zeros<arma::umat>(arma::size(X));  
    
    // Set all internal pixels to 0 
    if (only_boundary) {
        for (int x = 0; x < X.n_rows; x++) {
            for (int y = 0; y < X.n_cols; y++) {
                // Skip background pixel
                if (X(x, y) == 0) {
                    // do nothing 
                } else if (x == 0 || y == 0 || x == X.n_rows-1 || y == X.n_cols-1) {
                    // Keep image boundary pixels as is 
                    X_bound(x, y) = X(x, y); 
                } else if (
                    // Zero out internal pixels
                    X(x-1, y-1) == X(x, y) &&
                    X(x-1, y) == X(x, y) &&
                    X(x-1, y+1) == X(x, y) &&
                    X(x, y-1) == X(x, y) &&
                    X(x, y+1) == X(x, y) &&
                    X(x+1, y-1) == X(x, y) &&
                    X(x+1, y) == X(x, y) &&
                    X(x+1, y+1) == X(x, y) 
                ) {
                    // do nothing 
                } else {
                    X_bound(x, y) = X(x, y); 
                }
            }
        }
        for (int i = 1; i <= ncells; i++) {
            res[i-1] = ind2sub(size(X), find(X_bound == i)) + 1; 
        }        
    } else {
        for (int i = 1; i <= ncells; i++) {
            res[i-1] = ind2sub(size(X), find(X == i)) + 1; 
        }
        
    }
    
    return res;
}


