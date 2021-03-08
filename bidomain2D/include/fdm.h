/* 
   This file contains the routine to assemble the 
   finite difference matrices for the bidomain 
   in two spatial dimensions.  
*/

#ifndef FDM_H
#define FDM_H

#include "cellmodel.h"
#include <eigen3/Eigen/Sparse>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;



void add_coefficient(int id, const int i, const int j,
		     double w, std::vector<T>& coeffs,
		     const int nx, const int ny);
		     

/* coeffs is a vector of Eigen::Triplet to be filled with
       matrix nonzero coefficients
   sx/sy are the x/y direction diffusivity coefficients
   nx,ny,hx,hy are the domain's discretization parameters */
void build2DLaplacian(std::vector<T>& coeffs,
		      double sx, double sy,
		      int nx, int ny, double hx, double hy);

/*
  Assemble the matrices of the whole problem
  TODO: add a function overload in case
        the diffusion tensor is not diagonal
	(anisotropic diffusion) or if 
	the diffusion is function of space (non-homogeneous)
*/
void build_bidomain2D(SpMat& A, SpMat& B, SpMat& C,
		      const double si_x, const double se_x,
		      const double si_y, const double se_y,
		      const int nx, const int ny,
		      const double hx, const double hy);

#endif 
