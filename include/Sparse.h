#ifndef SPARSE_H
#define SPARSE_H

#include "Matrices.h"
#include <fstream>

class CSR
{
 private:
  unsigned int m_row;
  unsigned int m_col;
  double m_sparsity;
  
  std::vector<unsigned int> m_rows;
  std::vector<unsigned int> m_cols;
  std::vector<double> m_vals;


 public:
  CSR(const Matrix &mat);
  CSR(const CSR &mat);
  CSR(const std::string filename);
  CSR();
  
  // Sparse matrix/scalar operations
  CSR& operator=(const CSR& rhs);
  CSR operator*(const double c);
  CSR operator/(const double c);

  // Sparse matrix/vector product
  std::vector<double> operator*(const std::vector<double> &rhs);

  int get_row() const;
  int get_col() const;
  double get_sparsity() const;

  const std::vector<unsigned int> * get_rows() const;
  const std::vector<unsigned int> * get_cols() const;
  const std::vector<double> * get_vals() const;

  void Print_Details() const;
  void SaveMatrix(const std::string filename) const;  
};

#endif
