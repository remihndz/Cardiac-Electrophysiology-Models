/*
  Finite difference solution of the bidomain equations
  on a square domain. 
*/ 

#include "fdm.h"


int main()
{
  int nx=10, ny = 10;
  double hx = 1.0/nx, hy = 1.0/ny;
  double sx = 2.5e-4, sy = 2.5e-3;

  SpMat M;
  build_bidomain2D(M, sx,sy,sx,sy, nx,ny, hx,hy);

}
