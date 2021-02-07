#ifndef NUMERICAL_METHODS_CPP
#define NUMERICAL_METHODS_CPP

#include "NumericalMethods.h"
#include <assert.h>

/* Iterative methods.
   TODO: Add preconditioners
*/

Optimized CG(const CSR& A,
	     const std::vector<double>& b,
	     const std::vector<double>* x0,
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
      assert(x0->size()==n);
      xn = *x0;
    }

  // Iterative search for solution
  std::vector<double> Ap(n,0.0), r(n,0.0), history, p;
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
      r = r - alpha*(A*p);

      if (Iter%50 == 0) // Get rid of rounding errors
	r = b - A*xn;

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
      Iter++;
    }
  while(Iter < max_iter);

  Optimized sol {xn, residual, Iter, history};
  return sol;  
}

Optimized CG(const Matrix& A,
	     const std::vector<double>& b,
	     const std::vector<double>* x0,
	     const double eps, const int max_iter)

{
  int n = b.size();
  std::vector<double> xn;

  assert(A.m_cols == n);
  assert(A.m_rows == n);
  
  if (!x0)
    xn.assign(n, 0.0);
  else
    {
      assert(x0->size()==n);
      xn = *x0;
    }

  // Iterative search for solution
  std::vector<double> Ap(n,0.0), r(n,0.0), history, p;
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

      r = r - alpha*(A*p);
      if (Iter%50 == 0) // Get rid of rounding error
	r = b - A*xn;

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
      Iter++;
    }
  while(Iter < max_iter);

  Optimized sol {xn, residual, Iter, history};
  return sol;  
}

Arnoldi AnorldiProcess(const CSR& A,
		       const std::vector<std::vector<double>>& Q,
		       const int k)
{
  std::vector<double> q = A*(Q[k]);
  std::vector<double> h(k+1, 0.0); // Hessenberg matrice's new column
  
  for (int i = 0; i<k; i++)
    {
      h[i] = q*(Q[i]);
      q    = q - h[i] * (Q[i]);
    }
  h[k] = sqrt(q*q);

  Arnoldi HQ {h, q};
  return HQ;
}

Arnoldi AnorldiProcess(const Matrix& A,
		       const std::vector<std::vector<double>>& Q,
		       const int k)
{
  std::vector<double> q = A*(Q[k]);
  std::vector<double> h(k+1, 0.0); // Hessenberg matrice's new column
  
  for (int i = 0; i<k; i++)
    {
      h[i] = q*(Q[i]);
      q    = q - h[i] * (Q[i]);
    }
  h[k] = sqrt(q*q);

  Arnoldi HQ {h, q};
  return HQ;
}


Optimized GMRES(const CSR& A,
		const std::vector<double>& b,
		const std::vector<double>* x0,
		const double eps, const int max_iter,
		const int restart)
{
  // int n = b.size();
  // int m = max_iter;
  
  // assert(A.get_col() == n);
  // assert(A.get_row() == n);

  // if (!x0)
  //   xn.assign(n, 0.0);
  // else
  //   {
  //     assert(x0.size()==n);
  //     xn = x0;
  //   }

  // std::vector<double> r = b - A*xn;
  // double norm_b         = sqrt(b*b);
  // double norm_r         = sqrt(r*r);
  // double residual       = norm_r/norm_b;
  // std::vector<double>   beta(m+1, 0.0), history;
  // std::vector<std::vector<double>> H, Q;
  // int k = 0;

  // std::vector<double>
  
  // history.puch_back(residual);
  // beta[0] = norm_r;
  // Q[0] = b/norm_b;
  
  // while (k < m)
  //   {
      
  //     if (residual < eps)
  // 	{
  // 	  break;
  // 	}

  //     beta.push_back(-sn[k] * beta[k]);
  //     beta[k] *= cs[k];
  //     residual = abs(beta[k+1])/norm_b;
  //     history.push_back(residual);

  //     if (residual < eps)
  // 	break;           
  
  
  // Optimized sol {xn, residual, iter, history};
  // return sol;
}

Optimized GMRES(const Matrix& A,
		const std::vector<double>& b,
		const std::vector<double>* x0,
		const double eps, const int max_iter,
		const int restart)
{}

Optimized GaussSeidel(const CSR& A,
		      const std::vector<double>& b,
		      const std::vector<double>* x0,
		      const double eps, const int max_iter)
{}

Optimized GaussSeidel(const Matrix& A,
		      const std::vector<double>& b,
		      const std::vector<double>* x0,
		      const double eps, const int max_iter)
{}

/*
  Direct methods
*/

Optimized LU(const CSR& A,
	     const std::vector<double> b)
{}

Optimized LU(const Matrix& A,
	     const std::vector<double> b)
{}

Optimized Cholesky(const CSR& A,
		   const std::vector<double> b)
{}

Optimized Cholesky(const Matrix& A,
		   const std::vector<double> b)
{}





#endif
