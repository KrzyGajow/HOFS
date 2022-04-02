// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat mult_(const arma::mat& X, const arma::mat& Y) {

  arma::mat result( X * Y );

  return result;

}

// [[Rcpp::export]]
List array_mult( const arma::mat& X, const arma::mat& Y ) {

  int n = X.n_cols;
  List out(n);

  for( int i = 0; i < n; i++ ) {

    arma::mat result( X.col(i) * Y.row(i) );
    out[ i ] = result;

  }

  return out;

}

// [[Rcpp::export]]
NumericMatrix Reduce_cpp(List x) {

  int n = x.size();
  NumericMatrix result = as<NumericMatrix>(x[0]);
  
  for ( int i = 1; i < n; ++i ) {
  
    result += as<NumericMatrix>(x[i]);
    
  }
  
  return result;
  
}
