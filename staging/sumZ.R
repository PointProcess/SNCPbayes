library(inline)
library(Rcpp)

sumZ <- function(x,Cext,w)
    {
        m <- dim(Cext)[1]
        if(length(m)==0)
            {
                m <- 1
                Cext <- t(as.matrix(Cext))
            }
        return(sumZC(x,Cext,w))
    }

sumZC <- cppFunction("NumericVector sumZC(NumericMatrix x, NumericMatrix Cext,double w)
{
  int i,j;
  int n = x.nrow();
  int m = Cext.nrow();
  double a,b;
  double two_w2 = 2 * w * w;
  NumericVector result(n);

  for(i = 0; i < m; ++i)
      {
        for(j = 0; j < n; ++j)
          {
              a = x(j,0) - Cext(i,0);
              b = x(j,0) - Cext(i,0);
              result[j] += exp(- (a * a + b * b) / two_w2);
	}
    }
  return(result);

}")
