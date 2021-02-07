#ifndef MATRICES_H
#define MATRICES_H

#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <stdexcept>
#include <unordered_set>

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
double operator*(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator*(const std::vector<double>& v, const double c);
std::vector<double> operator*(const double c, const std::vector<double>& v);

class Matrix
{
 public:
  std::vector<std::vector<double>> m_mat;
  unsigned int m_rows;
  unsigned int m_cols;

  /* Creates a multidiagonal matrix of size rowsxcols
     with value vals[i] on the diagonal diags[i]
  */
  Matrix(const std::vector<double> vals,
	 const std::vector<int> diags,
	 const unsigned int rows,
	 const unsigned int cols);
  // Construct a rowsxcols matrix  with value init everywhere
  Matrix(unsigned int rows, unsigned int cols, const double& init);
  // Initiate matrix using another matrix
  Matrix(const Matrix& A);

  // Matrix/Matrix operations
  Matrix& operator=(const Matrix& rhs);
  Matrix operator+(const Matrix &rhs);
  Matrix operator-(const Matrix &rhs);

  // Matrix/Scalar operations
  Matrix operator*(const double c);
  Matrix operator/(const double c);  
		   
  // Matrix/Vector operations
  std::vector<double> operator*(const std::vector<double> &rhs);
  std::vector<double> operator*(const std::vector<double> &rhs) const;

  // Access individual elements
  double* operator()(const unsigned int row, const unsigned int col);

  void SwapRows(std::vector<unsigned int> &I,
		std::vector<unsigned int> &J);
  
};

#endif
