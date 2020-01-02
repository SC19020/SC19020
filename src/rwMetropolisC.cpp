#include <Rcpp.h>
using namespace Rcpp;

//' @title random walk Metropolis sampler using Rcpp
//' @description Use three inputs to generate random walk metropolis samples with the proposal distribution of normal distribution
//' @param sigma the first predictor (numeirc)
//' @param a the second predictor (numeric)
//' @param N the third predictor (numeric)
//' @return a random sample of size \code{N}
//' @examples
//' \dontrun{
//' rw <- rwMetropolisC(2,25,1000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector rwMetropolisC(double sigma, double a, int N) {
  NumericVector x(N);
  x[0] = a;
  int k = 0;
  for(int i = 1; i < N; i++){
    double u = as<double>(runif(1));
    double y = as<double>(rnorm(1, x[i-1], sigma));
    if(u <= exp(-abs(y)+abs(x[i-1]))) 
      x[i] = y;
    else {
      x[i] = x[i-1];
      k += 1;
    }
  }
  return x;
}