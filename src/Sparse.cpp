#ifndef SPARSE_CPP
#define SPARSE_CPP

#include "Sparse.h"
#include <assert.h>
#include <algorithm>

// Constructors
CSR::CSR(const Matrix &mat)
{
  m_row = mat.m_rows;
  m_col = mat.m_cols;

  m_rows.push_back(0);

  std::vector<std::vector<double>>::const_iterator row;
  std::vector<double>::const_iterator col;
  
  for (row = mat.m_mat.begin(); row != mat.m_mat.end(); row++)
    {
      int NbNonZero = 0, j = 0;
      
      for (col = row->begin(); col != row->end(); col++)
	{
	  double v = *col;
	  if (v != 0.0)
	    {
	      NbNonZero++;
	      m_vals.push_back(v);
	      m_cols.push_back(j);
	    }
	  j++;
	}
      NbNonZero = m_rows.back() + NbNonZero;
      m_rows.push_back(NbNonZero);
    }

  m_sparsity = 100.0*(1.0 - 1.0*m_rows.back()/(m_row*m_col));
  //this->Print_Details();     
}

CSR::CSR(const CSR &mat)
{
  m_row = mat.get_row();
  m_col = mat.get_col();
  m_sparsity = mat.get_sparsity();

  m_rows = *mat.get_rows();
  m_cols = *mat.get_cols();
  m_vals = *mat.get_vals();
}

CSR::CSR(const std::string filename)
{
  std::ifstream fin(filename);
  int n,m, nnz;
  
  fin >> n >> m >> nnz;

  m_row = n;
  m_col = m;
  m_rows.resize(n+1);
  m_cols.resize(nnz);
  m_vals.resize(nnz);
  
  if (n+1 <= nnz)
    {
      for (int i = 0; i < n+1; i++)
	{
	  unsigned int a,b;
	  double v;
	  fin >> v >> a >> b;
	  m_vals[i] = v;
	  m_cols[i] = a;
	  m_rows[i] = b;
	}
      for (int i = n+1; i < nnz; i++)
	{
	  unsigned int a,b;
	  double v;
	  fin >> v >> a >> b;
	  m_vals[i] = v;
	  m_cols[i] = a;
	}
    }

  else
        {
      for (int i = 0; i < nnz; i++)
	{
	  unsigned int a,b;
	  double v;
	  fin >> v >> a >> b;
	  m_vals[i] = v;
	  m_cols[i] = a;
	  m_rows[i] = b;
	}
      for (int i = nnz; i < n+1; i++)
	{
	  unsigned int a,b;
	  double v;
	  fin >> b;
	  m_rows[i] = b;
	}
    }

  fin.close();
  m_sparsity = 100.0*(1.0 - 1.0*m_rows.back()/(m_row*m_col));
  
  std::cout << "Matrix successfully read from " << filename << '\n';

  Print_Details(); // Prints a short summary of the matrix characteristics  
}

CSR::CSR()
{
  m_row = 0;
  m_col = 0;
  m_sparsity = 0;

  m_rows = {0};
  m_cols = {0};
  m_vals = {0};
  
}



// Allocator
CSR& CSR::operator=(const CSR& rhs)
{
  if (&rhs == this)
    return *this;

  m_row = rhs.get_row();
  m_col = rhs.get_col();
  m_sparsity = rhs.get_sparsity();

  m_rows = *rhs.get_rows();
  m_cols = *rhs.get_cols();
  m_vals = *rhs.get_vals();

  return *this;
}


CSR CSR::operator*(const double c)
{
  if (c==0)
    {
      std::cout << "Warning: Multiplication by 0 returns an empty CSR matrix!\n";
      CSR Res;
      return Res;
    }
  CSR Res(*this);
  Res.m_vals = (this->m_vals)*c;
  return Res;
}

CSR CSR::operator/(const double c)
{
  if (c==0)
    throw std::invalid_argument("Division by zero!");     

  double d = 1.0/c;

  return (*this)*d;
}

std::vector<double> CSR::operator*(const std::vector<double> &rhs)
{

  int n = this->m_row, m = this->m_col;
  if (m != rhs.size())
    throw std::invalid_argument("Shape mismatch between CSR and vector!");

  std::vector<double> res(n,0.0);
  std::vector<unsigned int>::const_iterator col = m_cols.begin();
  std::vector<double>::const_iterator val = m_vals.begin();

  for (int i = 0; i < n; i++)
    {
      int nnz = m_rows[i+1] - m_rows[i];

      for (int j = 0; j < nnz; j++)
	{
	  res[i] += (*val) * rhs[*col];
	  col++;
	  val++;
	}
    }
  return res;  
}

std::vector<double> CSR::operator*(const std::vector<double> &rhs) const
{

  int n = this->m_row, m = this->m_col;
  if (m != rhs.size())
    throw std::invalid_argument("Shape mismatch between CSR and vector!");

  std::vector<double> res(n,0.0);
  std::vector<unsigned int>::const_iterator col = m_cols.begin();
  std::vector<double>::const_iterator val = m_vals.begin();

  for (int i = 0; i < n; i++)
    {
      int nnz = m_rows[i+1] - m_rows[i];

      for (int j = 0; j < nnz; j++)
	{
	  res[i] += (*val) * rhs[*col];
	  col++;
	  val++;
	}
    }
  return res;  
}

	 
int CSR::get_row() const
{
  return this->m_row;
}

int CSR::get_col() const
{
  return this->m_col;
}

double CSR::get_sparsity() const
{
  return this->m_sparsity;
}

const std::vector<unsigned int>* CSR::get_rows() const
{
  return &(m_rows);
}

const std::vector<unsigned int>* CSR::get_cols() const
{
  return &(m_cols);
}

const std::vector<double>* CSR::get_vals() const
{
  return &(m_vals);
} 

void CSR::Print_Details() const
{
  std::cout << "Matrix of size\t(" << m_row << ',' << m_col << ").\n";
  std::cout << "Zero elements represent\t" << m_sparsity << "%\tof the matrix.\n";
  std::cout << "Storage type is CSR.\n";
}


void CSR::SaveMatrix(const std::string filename) const
{
  std::ofstream fout(filename);

  fout << m_row << " " << m_col << " " << m_rows.back() << "\n\n";

  if (m_rows.size() < m_cols.size())
    {
      for (int i = 0; i<m_rows.size(); i++)
	fout << m_vals[i] << " " << m_cols[i] << " " << m_rows[i] << '\n';

      for (int i = m_rows.size(); i < m_cols.size(); i++)
	fout << m_vals[i] << " " << m_cols[i] << "\n";
    }

  else
    {
      for (int i = 0; i<m_cols.size(); i++)
	fout << m_vals[i] << " " << m_cols[i] << " " << m_rows[i] << '\n';
      
      for (int i = m_cols.size(); i < m_rows.size(); i++)
	fout << m_rows[i] << "\n";
    }
  fout.close();
  std::cout << "Matrix has been saved in " << filename << '\n';
}

#endif
