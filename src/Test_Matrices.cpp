#include "Matrices.h"
#include "Sparse.h"

#include <ctime>
#include <cstdlib>

int main()
{

  {
    unsigned int n=5;
    int ret = std::system("mkdir -p ./Results/SpeedTests/");

    std::ofstream Crea;
    Crea.open("./Results/SpeedTests/Matrix_Creation.dat");
    std::ofstream Mult;
    Mult.open("./Results/SpeedTests/MatrixVect_Mult.dat");
  
    std::cout << "Testing the speed of matrix creation and matrix-vector product (see ./Results/SpeedTests/)..." << std::endl;
    for (n = 20; n<5000; n += 100)
      {
    	std::clock_t start = std::clock(); 
    	Matrix A(n,n,10.0);
    	std::clock_t end = std::clock();
    	long double time_elapsed_ms = 1000.0 * (end-start)/CLOCKS_PER_SEC;
    	Crea << n << " " << time_elapsed_ms << '\n';

    	std::vector<double> v(n, 2.0), u;
    	start = std::clock();
    	u = A*v;
    	end = std::clock();
    	time_elapsed_ms = 1000.0 * (end-start)/CLOCKS_PER_SEC;
    	Mult << n << " " << time_elapsed_ms << '\n';
      }

    Crea.close();
    Mult.close();    

    std::vector<int> diags{-1,0,1};
    std::vector<double> vals{-1.,2.,1.};
    n = 50;
    Matrix A(vals, diags, n, n);
  
    std::vector<unsigned int> I{0,2,3}, J{n-3,n-2,n-1};
    A.SwapRows(I,J);      
    std::cout << "SwapRows() works!\n";

    std::cout << "Testing if SwapRows messes with multiplication...\n" << std::endl;
    
    std::vector<double> v(n, 1.0), u, r(n, 2.0);
    
    r[I.back()] = 1.0;
    r[J[0]] = 3.0;
    u = A*v;
    
    bool Correct = true;
    int i=0;
    while (Correct && i<n)
      {
	Correct = u[i]==r[i];
	i++;
      }


    if (Correct)
      std::cout << "The multiplication stays correct after SwapRows()!\n" << std::endl;
    else
      std::cout << "Somethind went wrong at index " << i-1 << " in the result!\n" << std::endl;
  }

  std::cout << "\n\nTesting CSR storage for sparse matrices...\n" << std::endl;

  
  int n = 10;
  std::vector<int> diags{-1,0,1};
  std::vector<double> vals{-1.,2.,1.};
  Matrix A(vals, diags, n,n);
  CSR B(A);
  
  std::vector<double> v(n,1.0), u(n, 0.0);

  std::ofstream Crea;
  Crea.open("./Results/SpeedTests/CSR_Creation.dat");
  std::ofstream Mult;
  Mult.open("./Results/SpeedTests/CSRVect_Mult.dat");
  
  std::cout << "Testing the speed of matrix creation and CSR matrix-vector product (see ./Results/SpeedTests/)..." << std::endl;
  for (n = 20; n<5000; n += 100)
    {
      Matrix M(vals, diags, n,n);
      std::clock_t start = std::clock();
      CSR A(M);
      std::clock_t end = std::clock();
      long double time_elapsed_ms = 1000.0 * (end-start)/CLOCKS_PER_SEC;
      Crea << n << " " << time_elapsed_ms << '\n';
      
      std::vector<double> v(n, 2.0), u;
      start = std::clock();
      u = A*v;
      end = std::clock();
      time_elapsed_ms = 1000.0 * (end-start)/CLOCKS_PER_SEC;
      Mult << n << " " << time_elapsed_ms << '\n';
    }
  
  Crea.close();
  Mult.close();    
  
  
  return 0;
}
