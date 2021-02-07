#ifndef MATRICES_CPP
#define MATRICES_CPP 

#include "Matrices.h"
#include <assert.h>

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());

  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(),
		 std::back_inserter(result), std::plus<double>());
  return result;
}

double operator*(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());
  double result=0.0;
  
  result = std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
  return result;
}

std::vector<double> operator*(const std::vector<double>& v, const double c)
{
  std::vector<double> result = v;
  for (int i = 0; i<result.size(); i++)
    result[i] = result[i]*c;
  return result;
}

std::vector<double> operator*(const double c, const std::vector<double>& v)
{
  return v*c;
}

std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());

  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(),
		 std::back_inserter(result), std::minus<double>());
  return result;
}





// Matrix constructors

Matrix::Matrix(const unsigned int rows,
	       const unsigned int cols, const double& init)
{
  unsigned int i;
  
  m_rows = rows;
  m_cols = cols;
  m_mat.resize(rows);
  for (i = 0; i<rows; i++)
    m_mat[i].resize(cols, init);
}

Matrix::Matrix(const Matrix& rhs)
{
  this->m_mat = rhs.m_mat;
  this->m_rows = rhs.m_rows;
  this->m_cols = rhs.m_cols;
}


Matrix::Matrix(const std::vector<double> vals,
	       const std::vector<int> diags,
	       const unsigned int rows,
	       const unsigned int cols)
{
  if (vals.size() != diags.size())
    {
      throw std::invalid_argument("The number of value differs from the number of diagonals given!\n");
    }

  this->m_rows = rows;
  this->m_cols = cols;
  this->m_mat.resize(rows);

  std::vector<std::vector<double>>::iterator row;
  for (row = this->m_mat.begin(); row != this->m_mat.end(); row++)
    {
      row->resize(cols, 0.0);
    }

  for (int i = 0; i < rows; i++)
    {      
      for (unsigned int k = 0; k < diags.size(); k++)
	{
	  double val = vals[k];
	  int      l = diags[k];
	  
	  if (l >= 0 && i+l < cols)
	    {
	      this->m_mat[i][i+l] = val;
	    }

	  else if (l < 0 && i+l >= 0)
	    {
	      this->m_mat[i][i+l] = val;
	    }
	}
    }
}


Matrix& Matrix::operator=(const Matrix& rhs)
{  
  if (&rhs == this)
    return *this;


  std::vector<std::vector<double>>::const_iterator rowA;
  std::vector<std::vector<double>>::iterator rowB;
  
  m_rows = rhs.m_rows;
  m_cols = rhs.m_cols;

  m_mat.resize(m_rows);
  rowB = m_mat.begin();
  for (rowA = rhs.m_mat.begin(); rowA != rhs.m_mat.end(); rowA++)
    {
      *rowB = *rowA;
      rowB++;
    }
  return *this;
}
  
 
// Matrix/Matrix operations
Matrix Matrix::operator+(const Matrix& rhs)
{
  if (this->m_rows != rhs.m_rows || this->m_cols != rhs.m_cols)
    {
      std::cout << "Shape \t(" << this->m_rows << ","<< this->m_cols
		<< ") \t does not match shape\t("<<rhs.m_rows
		<< ","<<rhs.m_cols << "!\n";
      throw std::invalid_argument("Mismatch in matrices' dimensions.");
    }

  Matrix B(rhs);
  std::vector<std::vector<double>>::const_iterator rowA;
  std::vector<std::vector<double>>::iterator rowB;
    
  rowB = B.m_mat.begin();
  for (rowA = this->m_mat.begin(); rowA != this->m_mat.end(); rowA++)
    {
      *rowB = *rowB + *rowA;
      rowB++;
    }
  return B;
}

Matrix Matrix::operator-(const Matrix& rhs)
{
  if (this->m_rows != rhs.m_rows || this->m_cols != rhs.m_cols)
    {
      std::cout << "Shape \t(" << this->m_rows << ","<< this->m_cols
		<< ") \t does not match shape\t("<<rhs.m_rows
		<< ","<<rhs.m_cols << "!\n";
      throw std::invalid_argument("Mismatch in matrices' dimensions.");
    }

    Matrix B(rhs);
  std::vector<std::vector<double>>::const_iterator rowA;
  std::vector<std::vector<double>>::iterator rowB;
    
  rowB = B.m_mat.begin();
  for (rowA = this->m_mat.begin(); rowA != this->m_mat.end(); rowA++)
    {
      *rowB = *rowA - *rowB;
      rowB++;
    }
  return B;
}


// Matrix/Scalar operations
Matrix Matrix::operator*(const double c)
{
  std::vector<std::vector<double>>::iterator row;
  Matrix Ac = *this;
  
  for (row = Ac.m_mat.begin(); row != Ac.m_mat.end(); row++)
    {
      std::transform(row->begin(), row->end(), row->begin(), [c](double &x){return c*x;});
    }
  return Ac;
}

Matrix Matrix::operator/(const double c)
{
  if (c==0.0)
    throw std::invalid_argument("Division by zero!");
  double d = 1.0/c;
  return (*this)*d;  
}


// Matrix/Vector operations
std::vector<double> Matrix::operator*(const std::vector<double> &rhs)
{
  if (m_cols != rhs.size())
    {
      std::cout << "Shape\t" << m_cols << "\t does not match\t" << rhs.size() << '\n';
      throw std::invalid_argument("Mismatch in arrays shapes.\n");
    }

  std::vector<std::vector<double>>::const_iterator row;  
  std::vector<double> v(m_rows, 0.0);
  unsigned int i = 0;
  
  for (row = this->m_mat.begin(); row != this->m_mat.end(); row++)
    {
      v[i] = (*row) * rhs;
      i++;
    }

  return v;
}

std::vector<double> Matrix::operator*(const std::vector<double> &rhs) const
{
  if (m_cols != rhs.size())
    {
      std::cout << "Shape\t" << m_cols << "\t does not match\t" << rhs.size() << '\n';
      throw std::invalid_argument("Mismatch in arrays shapes.\n");
    }

  std::vector<std::vector<double>>::const_iterator row;  
  std::vector<double> v(m_rows, 0.0);
  unsigned int i = 0;
  
  for (row = this->m_mat.begin(); row != this->m_mat.end(); row++)
    {
      v[i] = (*row) * rhs;
      i++;
    }

  return v;
}

// Access individual elements
double * Matrix::operator()(const unsigned int i, const unsigned int j)
{
  if (i >= m_rows || j >= m_cols)
    {
      std::cout << "Matrix of size (" << m_rows << ","<< m_cols
		<< ") and (i,j)=("<<i<<","<<j<<'\n';
      throw std::invalid_argument("Out of bound in Matrix(i,j).\n"); 
    }
  return &m_mat[i][j];
}


void Matrix::SwapRows(std::vector<unsigned int> &I, std::vector<unsigned int> &J)
{
  if (I.size() != J.size())
    {
      throw std::invalid_argument("The number of indexes in I and J don't match!\n");
    }

  // Remove potential duplicates
  std::vector<unsigned int>::iterator itr = I.begin();
  std::unordered_set<unsigned int> s;

  for (auto curr = I.begin(); curr != I.end(); ++curr) {
    if (s.insert(*curr).second)
      *itr++ = *curr;
  }

  I.erase(itr, I.end());
  
  itr = J.begin();
  for (auto curr = J.begin(); curr != J.end(); ++curr) {
    if (s.insert(*curr).second)
      *itr++ = *curr;
  }
  
  J.erase(itr, J.end());

  for (unsigned int k = 0; k<I.size(); k++)
    {
      unsigned int i = I[k], j = J[k];
      this->m_mat[i].swap(this->m_mat[j]);
    }
}



#endif
