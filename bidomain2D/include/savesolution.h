/* Routines to save solutions in vector form 
   into either vtk or txt format */
#ifndef SAVESOLUTION_H
#define SAVESOLUTION_H


#include <fstream>
#include <eigen3/Eigen/Dense>
#include <string>
#include <iostream>

typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Mesh;

/* 
   sol is a vector of the solution to save
   X is an array of triplet (x_k, y_k, z_k)
   giving the position of node k
   sol[k] = sol(x_k,y_k,z_k)
*/
int savevtk(std::string filename, Eigen::VectorXd& sol,
	    Mesh& X);
int savetxt(std::string filename, Eigen::VectorXd& sol,
	    Mesh& X);



#endif
