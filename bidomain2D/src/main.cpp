/*
  Finite difference solution of the bidomain equations
  on a square domain. 
*/ 

#include "fdm.h"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "savesolution.h"
#include <eigen3/Eigen/SparseCholesky>

int main()
{
  const int nx=50, ny = 50;
  double lx = 1, ly = 1;
  double hx = lx/(nx-1), hy = ly/(ny-1);
  double sx = 2.5e-3, sy = 1e-4;

  std::string rootResults = "./Results/bidomain2D/";
  std::string s = "rm -r " + rootResults;
  const char *command = s.c_str();
  system(command);
  s = "mkdir -p " + rootResults;
  const char *command2 = s.c_str();
  system(command2);



  // Spatial discretisation
  std::cout << "Building the finite differences matrices... " << std::endl;
  SpMat A(nx*ny, nx*ny), B(nx*ny, nx*ny), C(nx*ny, nx*ny);
  build_bidomain2D(A,B,C, sx,sy,sx,sy, nx,ny, hx,hy);

  Eigen::SimplicialLDLT<SpMat> Solver_ue;
  Solver_ue.compute(C); // Factorizes C
  if (Solver_ue.info()==Eigen::Success)
    std::cout << "Decomposition succeeded." << std::endl;
  else
    std::cout << "Decomposition failed." << std::endl;
  std::cout << "Spatial discretisation done!" << std::endl;
  
  // Time parameters
  double t = 0, dt = 0.5, tmax = 3000;
  int time_iter = 0, save_iter = 0;

  // Initialization
  Vec vm(nx*ny), vmnext(nx*ny), ue(nx*ny), aux(2*nx*ny);

  const double Am = 200, Cm = 1e-3;
  Vec w(nx*ny), wnext(nx*ny);
  MS ionic(w, vm);

  // const double Am = 200, Cm = 1;
  // MatW  w(nx*ny,7), wnext(nx*ny,7);
  // BR ionic(w, vm);
  
  Vec F = Eigen::ArrayXd::Zero(nx*ny) ;
  Mesh mesh(nx*ny, 3);
  
  for (int i = 0; i < nx*ny; i++)
    {
      double x = (i%nx)*hx, y = (i/nx)*hy;
      mesh.row(i) << x,y,0.0;
      if (sqrt(x*x+y*y) < 5*hx)
	vm[i] = 20.0;
      ue[i] = 0.0;
    }

  {
    std::string savefile = rootResults+"vm_"+std::to_string(save_iter)+".txt";
    savetxt(savefile, vm, mesh);
    std::string savefileue = rootResults+"ue_"+std::to_string(save_iter)+".txt";
    savetxt(savefileue, ue, mesh);
    save_iter++;
  }


  // Time loop 
  while (t < tmax)
    {

      if (mod(t,10)<dt/2. || mod(t,10)>10-dt/2.)
	{
	  std::cout << "======== Time = " << t << "ms ==========" << std::endl;
	}

      wnext = w + dt*ionic.g(w, vm);
      w = wnext;
      F = dt*ionic.Iion(w, vm)/Cm;
      
      vmnext = vm + (dt/Am)*(A*vm + B*ue) + dt*F;
      ue = Solver_ue.solve(-1.0*(B.transpose()*vm));
      vm = vmnext;
      
      if (time_iter%10==0)
	{
	  std::string savefilevm = rootResults+"vm_"+std::to_string(save_iter)+".txt";
	  savetxt(savefilevm, vm, mesh);
	  std::string savefileue = rootResults+"ue_"+std::to_string(save_iter)+".txt";
	  savetxt(savefileue, ue, mesh);
	  save_iter++;
	}

      if (mod(t,700)<dt)
	{
	  for (int i = 0; i < nx*ny; i++)
	    {
	      double x = (i%nx)*hx, y = (i/nx)*hy;
	      if (sqrt(x*x+y*y) < 5*hx)
		vm[i] = 20.0;
	    }
	}

      // if (t > 500 and t < 600+3*dt)
      // 	{
      // 	  for (int i = 0; i < nx*ny; i++)
      // 	    {
      // 	      double x = (i%nx)*hx, y = (i/nx)*hy;
      // 	      if (sqrt(x*x+pow(y-0.4,2)) < 5*hx)
      // 		vm[i] = 20.0;
      // 	    }
      // 	}
	
      
      time_iter++;
      t+=dt;
    }  
}
