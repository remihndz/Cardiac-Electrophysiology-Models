#include "Matrices.h"
#include "Sparse.h"
#include "NumericalMethods.h"

#include <fstream>
#include <cmath>

double f(double x,double y)
{
  double f = -exp(-x*x-y*y)*((4*pow(x,4) - 4*pow(x,3) - 10*pow(x,2) + 6*x + 2)*y*(y-1)
			     + (4*pow(y,4) - 4*pow(y,3) - 10*pow(y,2) + 6*y + 2)*x*(x-1));
  return f;
}

double u(double x, double y)
{
  return exp(-x*x-y*y)*(x-1)*(y-1)*x*y;
}
  

int main()
{
  /* Creates a classical Laplacian finite differences
     matrix and solve the Laplace equation 
        u_xx + u_yy = 0    0<x,y<1 
	u = 1              x = 0,1
	u_y = 0            y = 0,1
  */

  int n;
  std::cout << "Enter the value of n: ";
  std::cin >> n;
  
  int nn = n+1, N = nn*nn;
  double h = 1./n, hh = h*h;
  std::vector<double> vals {-1.0/hh,-1.0/hh,4.0/hh,-1.0/hh,-1.0/hh};
  std::vector<int>   diags {-nn, -1, 0, 1, nn};

  Matrix A(vals, diags, N, N);	// Discrete Laplace

  // Add boundary condition for x=0 and x=1
  for (int i = 1; i<n; i++)
    {
      // x = 0
      int RowIndex = i*(n+1);
      *A(RowIndex, RowIndex) = 1.0;
      *A(RowIndex, RowIndex-1) = 0.0;
      *A(RowIndex, RowIndex+1) = 0.0;
      *A(RowIndex, RowIndex-n-1) = 0.0;
      *A(RowIndex, RowIndex+n+1) = 0.0;
      
      // x = 1
      RowIndex += n;
      *A(RowIndex, RowIndex) = 1.0;
      *A(RowIndex, RowIndex-1) = 0.0;
      *A(RowIndex, RowIndex+1) = 0.0;
      *A(RowIndex, RowIndex-n-1) = 0.0;
      *A(RowIndex, RowIndex+n+1) = 0.0;
    }

  for (int i = 1; i<n; i++)
    {
      // y = 0
      *A(i, i) = 1.0;
      *A(i, i-1) = 0.0;
      *A(i, i+1) = 0.0;
      *A(i, i+n+1) = 0.0;

      // y = 1
      *A(N-1-i, N-1-i) = 1.0;
      *A(N-1-i, N-1-i-1) = 0.0;
      *A(N-1-i, N-1-i+1) = 0.0;
      *A(N-1-i, N-1-i-n-1) = 0.0;
    }

  // Corner points
  {
  *A(0,0)     = 1.0;
  *A(0,1)     = 0.0;
  *A(0,n+1)   = 0.0;

  *A(n,n)     = 1.0;
  *A(n,n+1)   = 0.0;
  *A(n,n-1)   = 0.0;
  *A(n,n+n+1) = 0.0;

  *A(N-1, N-1)   = 1.0;
  *A(N-1, N-2)   = 0.0;
  *A(N-1, N-n-2) = 0.0;

  *A(N-1-n, N-1-n)     = 1.0;
  *A(N-1-n, N-1-n+1)   = 0.0;
  *A(N-1-n, N-1-n-1)   = 0.0;
  *A(N-1-n, N-1-n-n-1) = 0.0;
  }


  std::vector<double> F;
  
  for (int i = 0; i<n+1; i++)
    {
      for (int j = 0; j<n+1; j++)
  	{
	  if (i == 0 || i == n || j == 0 || j == n)
	    F.push_back(0.0);
	  else
	    F.push_back(f(i*h, j*h));
  	}
    }

  
  CSR M(A);

  std::vector<double>* x0 = nullptr;
  Optimized sol;

  sol = CG(M,F);

  std::cout << "Conjugate gradient reached residual of " << sol.m_fun << " after " << sol.m_evals << std::endl;

  std::ofstream fout("./Results/TestLaplaceEquation/u.dat");

  int k = 0;
  for (int i = 0; i < n+1; i++)
    {
      for (int j = 0; j < n+1; j++)
	{
	  {
	    fout << i*h << " " << j*h << " " << sol.m_x[k] << '\n';
	    k++;
	  }
	}
    } 

  fout.close();
  fout.open("./Results/TestLaplaceEquation/history.dat");

  for (int i = 0; i<sol.m_history.size(); i++)
    fout << sol.m_history[i] << '\n';

  fout.close();

  fout.open("./Results/TestLaplaceEquation/F.dat");

  k = 0;
  for (int i = 0; i < n+1; i++)
    {
      for (int j = 0; j < n+1; j++)
	{
	  {
	    fout << i*h << " " << j*h << " " << F[k] << '\n';
	    k++;
	  }
	}
    } 

  fout.close();
  fout.open("./Results/TestLaplaceEquation/u_ex.dat");

  k = 0;
  for (int i = 0; i < n+1; i++)
    {
      for (int j = 0; j < n+1; j++)
	{
	  double x = i*h, y = j*h;
	  fout << x << " " << y << " " << pow(sol.m_x[k]-u(x,y),2) << '\n';
	  k++;
	}
	  
    }   
  
  return 0;
}
