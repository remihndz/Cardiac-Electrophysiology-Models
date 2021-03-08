/*
  This file contains the routine to assemble the 
  finite difference matrices for the bidomain 
  in two spatial dimensions. 
  Neumann homogeneous boundary conditions are 
  always assumed and implemented with 
  ghost points method.
*/


#ifndef FDM_CPP
#define FDM_CPP

#include "fdm.h"
#include <iostream>

void add_coefficient(int id, const int i, const int j,
		     double w, std::vector<T>& coeffs,
		     const int nx, const int ny)
{
  int id1 = i + j*nx; // Column in the matrix
  double v = 0.0;
  int k = id%nx, l = id/ny;
  
  v = w;
  
  if (id1 != id)
    {
      if      (i==1 && k == 0)   v = 2*w;
      else if (i==nx-2 && k==nx-1) v = 2*w;
      else if (j==1 && l == 0)   v = 2*w;
      else if (j==ny-2 && l==ny-1) v = 2*w;
    }

  // Divide to keep symetry of the matrix? 
  if (k == 0 || k == nx-1)
    v /= 2.0;
  if (l == 0 || l == ny-1)
    v /= 2.0; 

  if (id1 != id)
    {
      if      (i==-1 && k == 0) v = 0.0;
      else if (i==nx && k==nx-1)  v = 0.0;
      else if (j==-1 && l == 0) v = 0.0;
      else if (j==ny && l==ny-1)  v = 0.0;
    } 
  
  if (id < nx*ny && id1 < nx*ny && id > -1 && id1 > -1 && v!=0.0)
    coeffs.push_back(T(id,id1,v));
}

/* coeffs is a vector of Eigen::Triplet to be filled with
       matrix nonzero coefficients
   sx/sy are the x/y direction diffusivity coefficients
   nx,ny,hx,hy are the domain's discretization parameters */
void build2DLaplacian(std::vector<T>& coeffs,
		      double sx, double sy,
		      int nx, int ny, double hx, double hy)
{
  for (int j = 0; j < nx; j++)
    {
      for (int i = 0; i < nx; i++)
	{
	  int id = i + j*nx; // Line in the matrix
	  add_coefficient(id, i-1, j, sx/hx/hx, coeffs, nx, ny);
	  add_coefficient(id, i+1, j, sx/hx/hx, coeffs, nx, ny);
	  add_coefficient(id, i, j, -2*sx/hx/hx -2*sy/hy/hy, coeffs, nx, ny);
	  add_coefficient(id, i, j-1, sy/hy/hy, coeffs, nx, ny);
	  add_coefficient(id, i, j+1, sy/hy/hy, coeffs, nx, ny);	  
	}
    }
}

/*
  Assemble the matrix of the whole problem
  TODO: add a function overload in case
        the diffusion tensor is not diagonal
	(anisotropic diffusion) or if 
	the diffusion is function of space.
	In either cases, need to implement 
	the finite element matrix for the
	second order term of the equation
	(namely ∇·(σ·∇·))
*/
void build_bidomain2D(SpMat&  A, SpMat& B, SpMat& C,
		      const double si_x, const double se_x,
		      const double si_y, const double se_y,
		      const int nx, const int ny,
		      const double hx, const double hy)
{
  /*
    Build the semi-discrete system:
      Aₘ(∂ₜVₘ + I(t,Vₘ)) = AVₘ + Bφₑ
      0                  = BᵀVₘ + Cu
  */
  std::vector<T> coeffs;
  
  build2DLaplacian(coeffs, si_x,si_y,nx,ny,hx,hy);  

  A.setFromTriplets(coeffs.begin(), coeffs.end());
  B.setFromTriplets(coeffs.begin(), coeffs.end());

  coeffs.empty();
  build2DLaplacian(coeffs, si_x+se_x, si_y+se_y, nx, ny, hx,hy);
  C.setFromTriplets(coeffs.begin(), coeffs.end());  
}


#endif
