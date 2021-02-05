#ifndef NUMERICAL_METHODS_CPP
#define NUMERICAL_METHODS_CPP

#include "NumericalMethods.h"
#include <assert.h>

/* Iterative methods.
   TODO: Add preconditioners
*/

Optimized CG(const CSR& A,
	     const std::vector<double>& b,
	     const std::vector<double>& x0 = NULL,
	     const double eps, const int max_iter)
{
  int n = b.size();
  std::vector<double> xn;

  assert(A.get_col() == n);
  assert(A.get_row() == n);
  
  if (!x0)
    xn.assign(n, 0.0);
  else
    {
      assert(x0.size()==n);
      xn = x0;
    }

  // Iterative search for solution
  std::vector<double> Ap(n,0.0), r(n,0.0), history;
  double residual, alpha, beta;
  int Iter = 0;
  
  r = b - A*xn;
  p = r;
  residual = r*r;
  history.push_back(sqrt(residual));

  do 
    {
      Ap = A*p;
      alpha = residual/(p*Ap);
      xn = xn + p*alpha;

      if (k%50 == 0) // Get rid of rounding error
	r = b - A*xnext;

      residual = r*r;
      if (sqrt(residual) < eps)
	{
	  history.push_back(sqrt(residual));
	  break;
	}

      beta = residual/pow(history.back(), 2);
      history.push_back(sqrt(residual));
      
      p = r + p*beta;
      history.push_back(sqrt(residual));
      k++;
    }
  while(k < max_iter);

  Optimized sol {xn, residual, k, history};
  return sol;  
}

Optimized CG(const Matrix& A,
	     const std::vector<double>& b,
	     const std::vector<double>& x0 = NULL,
	     const double eps, const int max_iter)
{
  int n = b.size();
  std::vector<double> xn;

  assert(A.get_col() == n);
  assert(A.get_row() == n);
  
  if (!x0)
    xn.assign(n, 0.0);
  else
    {
      assert(x0.size()==n);
      xn = x0;
    }

  // Iterative search for solution
  std::vector<double> Ap(n,0.0), r(n,0.0), history;
  double residual, alpha, beta;
  int Iter = 0;
  
  r = b - A*xn;
  p = r;
  residual = r*r;
  history.push_back(sqrt(residual));

  do 
    {
      Ap = A*p;
      alpha = residual/(p*Ap);
      xn = xn + p*alpha;

      if (k%50 == 0) // Get rid of rounding error
	r = b - A*xnext;

      residual = r*r;
      if (sqrt(residual) < eps)
	{
	  history.push_back(sqrt(residual));
	  break;
	}

      beta = residual/pow(history.back(), 2);
      history.push_back(sqrt(residual));
      
      p = r + p*beta;
      history.push_back(sqrt(residual));
      k++;
    }
  while(k < max_iter);

  Optimized sol {xn, residual, k, history};
  return sol;  
}


Optimized GMRES(const CSR& A,
		const std::vector<double>& b,
		const std::vector<double>& x0,
		const double eps, const int max_iter);

Optimized GMRES(const Matrix& A,
		const std::vector<double>& b,
		const std::vector<double>& x0,
		const double eps, const int max_iter);

Optimized GaussSeidel(const CSR& A,
		      const std::vector<double>& b,
		      const std::vector<double>& x0,
		      const double eps, const int max_iter);

Optimized GaussSeidel(const Matrix& A,
		      const std::vector<double>& b,
		      const std::vector<double>& x0,
		      const double eps, const int max_iter);

/*
  Direct methods
*/

Optimized LU(const CSR& A,
	     const std::vector<double> b);

Optimized LU(const Matrix& A,
	     const std::vector<double> b);

Optimized Cholesky(const CSR& A,
		   const std::vector<double> b);

Optimized Cholesky(const Matrix& A,
		   const std::vector<double> b);





#endif
