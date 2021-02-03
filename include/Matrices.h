#ifndef MATRICES_H
#define MATRICES_H

#include <vector>
#include <iostream>
#include <algorithm>
#include <stdexcept>


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
  
  void MoveRows(const std::vector<int> &I, const std::vector<int> &J);
  void CreateDiagonalMatrix(const std::vector<double> vals,
			    const std::vector<int> diags);
  
  Matrix(unsigned int rows, unsigned int cols, const double& init);
  Matrix(const Matrix& A);

  Matrix& operator=(const Matrix& rhs);
  
  // Matrix/Matrix operations
  Matrix operator+(const Matrix &rhs);
  Matrix operator-(const Matrix &rhs);

  // Matrix/Scalar operations
  Matrix operator*(const double c);
  Matrix operator/(const double c);  
		   
  // Matrix/Vector operations
  std::vector<double> operator*(const std::vector<double> &rhs);

  // Access individual elements
  double operator()(const unsigned int row, const unsigned int col);
};

#endif
