#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "Sparse.h"
#include <math.h>

struct Optimized
{
  std::vector<double> m_x;
  std::vector<double> m_history;
  double              m_fun;
  int                 m_evals;
};

/* Iterative methods.
   TODO: Add preconditioners
*/

Optimized CG(const CSR& A,
	     const std::vector<double>& b,
	     const std::vector<double>& x0 = NULL,
	     const double eps = 1e-6, const int max_iter = 1000);

Optimized CG(const Matrix& A,
	     const std::vector<double>& b,
	     const std::vector<double>& x0 = NULL,
	     const double eps = 1e-6, const int max_iter = 1000);


Optimized GMRES(const CSR& A,
		const std::vector<double>& b,
		const std::vector<double>& x0,
		const double eps = 1e-6, const int max_iter = 20);

Optimized GMRES(const Matrix& A,
		const std::vector<double>& b,
		const std::vector<double>& x0 = NULL,
		const double eps = 1e-6, const int max_iter = 20);

Optimized GaussSeidel(const CSR& A,
		      const std::vector<double>& b,
		      const std::vector<double>& x0 = NULL,
		      const double eps = 1e-6, const int max_iter = 1000);

Optimized GaussSeidel(const Matrix& A,
		      const std::vector<double>& b,
		      const std::vector<double>& x0 = NULL,
		      const double eps = 1e-6, const int max_iter = 1000);

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
