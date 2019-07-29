/*#include <Rmath.h>

void sumZ_aniso_R(double *x1, double *x2, int *n, double *Cext1, double* Cext2, 
		  int *m, double *w, double *rho, double *result)
{
  int i,j;
  double temp;
  double a,b;
  double two_w2 = 2 * *w * *w * (1 - *rho * *rho);
  for(i = 0; i < *n; i++)result[i] = 0;
  for(i = 0; i < *m; i++)
    {
      for(j = 0; j < *n; j++)
	{
	  a = x1[j] - Cext1[i];
	  b = x2[j] - Cext2[i];
	  result[j] += exp(- (a * a + b * b - 2 * *rho * a * b) / two_w2);
	}
    }
  return;
}
*/
