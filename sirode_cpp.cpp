#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
List sir_cpp(double time, NumericVector y, NumericVector params) {
  int n = y.size();
  NumericVector dydt( n );
  double beta = params[ 0 ];
  double gamma = params[ 1 ];
  double N = 0;
  for(int i = 0; i < n; ++i) {
    N += y[i];
  }
  double foi = beta * y[1]/N;  
  
  dydt[0] = - foi*y[0]; //S
  dydt[1] = foi*y[0] - gamma*y[1];
  dydt[2] = gamma*y[1];

  return List::create(dydt);
}

// [[Rcpp::export]]
List sird_cpp(double time, NumericVector y, NumericVector params) {
  int n = y.size();
  NumericVector dydt( n );
  double beta = params[ 0 ];
  double gamma = params[ 1 ];
  double N = 0;
  for(int i = 0; i < n; ++i) {
    N += y[i];
  }
  double foi = beta * y[1]/N;  
  double mu = 1/(60*365);
  dydt[0] = N*mu - foi*y[0] - mu*y[0]; //S
  dydt[1] = foi*y[0] - gamma*y[1] - mu*y[1];
  dydt[2] = gamma*y[1] - mu*y[2];
  
  return List::create(dydt);
}
