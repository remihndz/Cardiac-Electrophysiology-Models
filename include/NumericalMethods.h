#ifndef NUMERICAL_METHODS_H
#define NUMERICAL_METHODS_H

#include "Sparse.h"
#include <math.h>

struct Optimized
{
  std::vector<double> m_x;
  double              m_fun;
  int                 m_evals;
  std::vector<double> m_history;
};

struct Arnoldi
{
  std::vector<double> h;
  std::vector<double> q;
};

/* Iterative methods.
   TODO: Add preconditioners
*/

Optimized CG(const CSR& A,
	     const std::vector<double>& b,
	     const std::vector<double>* x0 = nullptr,
	     const double eps = 1e-6, const int max_iter = 1000);

Optimized CG(const Matrix& A,
	     const std::vector<double>& b,
	     const std::vector<double>* x0 = nullptr,
	     const double eps = 1e-6, const int max_iter = 1000);


Arnoldi AnorldiProcess(const CSR& A,
		       const std::vector<std::vector<double>>& Q,
		       const int k);

Arnoldi AnorldiProcess(const Matrix& A,
		       const std::vector<std::vector<double>>& Q,
		       const int k);



Optimized GMRES(const CSR& A,
		const std::vector<double>& b,
		const std::vector<double>* x0 = nullptr,
		const double eps = 1e-6, const int max_iter = 200,
		const int restart = 10);

Optimized GMRES(const Matrix& A,
		const std::vector<double>& b,
		const std::vector<double>* x0 = nullptr,
		const double eps = 1e-6, const int max_iter = 200,
		const int restart = 10);

Optimized GaussSeidel(const CSR& A,
		      const std::vector<double>& b,
		      const std::vector<double>* x0 = nullptr,
		      const double eps = 1e-6, const int max_iter = 1000);

Optimized GaussSeidel(const Matrix& A,
		      const std::vector<double>* b,
		      const std::vector<double>* x0 = nullptr,
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
