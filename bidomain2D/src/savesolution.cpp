/* Routines to save solutions in vector form 
   into either vtk or txt format */
#ifndef SAVESOLUTION_CPP
#define SAVESOLUTION_CPP

#include "savesolution.h"

/* 
   sol is a vector of the solution to save
   X is an array of triplet (x_k, y_k, z_k)
   giving the position of node k
   sol[k] = sol(X(k,0),X(k,1),X(k,2))
*/

int savevtk(std::string filename, Eigen::VectorXd& sol,
	     Mesh& X)
{
  std::ofstream fout(filename);

  if (!fout.is_open())
    {
      std::cout << "The output file could not be opened. The file location was: "
		<< filename << std::endl;
      return 1;
    }

  int n = sol.size();
  
  fout << "#vtk DataFile Version 2.0" << std::endl;
  fout << "vtk format created via C++" << std::endl;
  fout << "ASCII" << std::endl;
  fout << "DATASET UNSTRUCTURED_GRID " << std::endl;

  fout << "POINTS " << n << " double" << std::endl;

  for (int k = 0; k < n; k++)
    fout << X(k,0) << " " << X(k,1) << " " << X(k,2) << std::endl;

  fout << "POINT_DATA " << n << std::endl;
  fout << "FIELD" << " fieldata " << " 1 " << std::endl;
  fout << "u " << " 1 " << n << " double " << std::endl;
  for (int k = 0; k < n; k++)
    fout << sol[k] << std::endl;

  fout.close();
  
  return 0;
}
  
int savetxt(std::string filename, Eigen::VectorXd& sol,
	     Mesh& X)
{
  std::ofstream fout(filename);

  if (!fout.is_open())
    {
      std::cout << "The output file could not be opened. The file location was: "
		<< filename << std::endl;
      return 1;
    }
  
  int n = sol.size();

  for (int k = 0; k < n; k++)
    fout << X(k,0) << " " << X(k,1) << " " << X(k,2) << " " << sol[k] << std::endl;

  fout.close();
  return 0;
}

#endif
