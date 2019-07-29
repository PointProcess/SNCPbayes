#include <Rcpp.h>
#include <Rmath.h>

using namespace Rcpp;

//' Sums Z?
//' @description Alex isn't quite clear what this does but needed to benchmark it so he started an example.
//' @param x An n x 2 matrix of locations of points
//' @param Cext An m x 2 matrix of centers
//' @param w What appears to be a standard deviation
//' @export
//' @return A vector of length n which has something to do with the m cluster centers
//' @examples
//' set.seed(1)
//' n <- 1e3
//' m <- 5
//' x <- matrix(rnorm(n * 2),n,2)
//' Cext <- matrix(rnorm(m * 2), m, 2)
//' w <- .01
//' a <- sumZ(x, Cext, w)
// [[Rcpp::export]]
NumericVector sumZ(NumericMatrix X, NumericMatrix Cext, double w)
{
  int n = X.nrow();
  int m = Cext.nrow();
  NumericVector result(n);
  int i,j;
  double temp;
  double a,b;
  double two_w2 = 2 * w * w;
  for(i = 0; i < m; ++i)
    {
      for(j = 0; j < n; ++j)
	{
	  a = X(j,0) - Cext(i,0);
	  b = X(j,1) - Cext(i,1);
	  result[j] += exp(-(a * a + b * b) / two_w2);
	}
    }
  return result;
}
