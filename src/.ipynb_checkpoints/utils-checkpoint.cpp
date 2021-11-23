#include <Rcpp.h>
using namespace std; 
using namespace Rcpp; 

// [[Rcpp::export]]
IntegerMatrix triplets_to_pairs(const IntegerMatrix& triplets) {
    unsigned N = triplets.nrow(); 
    IntegerMatrix res(6*N, 2);
    for (int i = 0; i < N; i++) {
        cout << "i" << endl;
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

