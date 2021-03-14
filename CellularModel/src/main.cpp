#include "cellmodel.h"

#include <fstream>
#include <iostream>


int main()
{
  /*
    Runs a simulation using the Mitchel and Schaeffer
    two currents model
  */
  {
    system("mkdir -p ./Results/MS/");
    
    double dt = 0.1; // For i_stim to work, dt <= 0.5
    double tmax = 1000, time = 0.0;
    std::vector<double> T;
    
    double vm, w; // Solutions at time t
    std::vector<double> V, W; // Solutions at all time
    MS ms;		    // Model used for simulation
    
    // Initial state
    vm = ms.vmin;
    w  = 1./pow(ms.vmax-ms.vmin, 2);
    
    V.push_back(vm);
    W.push_back(w);
    T.push_back(time);
    
    while (time < tmax + dt)
      {
  	w = w + dt * ms.g(w, vm);
  	vm = vm + dt * ms.Itot(w, vm, time);
      
  	time+=dt;
	
  	V.push_back(vm);
  	W.push_back(w);
  	T.push_back(time);
	
  	if (mod(time, 50) < dt)
  	  {
  	    std::cout << "Simulation at time: " << time
  		      << "ms" << std::endl;
  	  }
	
      }
    
    std::ofstream foutv("./Results/MS/vm.dat");
    std::ofstream foutw("./Results/MS/w.dat");
    
    for (int k = 0; k < T.size(); k++)
      {
  	foutv << T[k] << " " << V[k] << std::endl;
  	foutw << T[k] << " " << W[k] << std::endl;
      }
    foutv.close();
    foutw.close();
  }


  
  // Test of the Beeler Reuter model
  std::cout << "Simulating using Beeler-Reuter's model..." << std::endl;

  system("rm ./Results/BR/*.dat");
  system("mkdir -p ./Results/BR/");
  double dt = 0.01; // For i_stim to work, dt <= 0.5
  double tmax = 800, time = 0.0;
  std::vector<double> T;
    
  double vm; // Solutions at time t
  std::vector<double> w;
  std::vector<double> V, I;    // Solutions at all time
  std::vector<std::vector<double>> W, G, ions;
  BR br(w, vm);		    // Model used for simulation
  
  V.push_back(vm);
  W.push_back(w);
  T.push_back(time);
  G.push_back(br.g(w,vm));
  I.push_back(br.Itot(w,vm, time));
  
  ions.push_back({br.i_K1(vm), br.i_x1(vm, w), br.i_Na(vm,w), br.i_s(vm,w)});
  
  std::cout << "Initial state sol = [" << vm;
  for (int i = 0; i<w.size(); i++)
    std::cout << ", " << w[i];
  std::cout << std::endl;
  
  while (time < tmax + dt)
    {
      I.push_back(br.Itot(w,vm, time));
      w = w + dt * br.g(w, vm);
      vm = vm + dt * br.Itot(w, vm, time);
      
      time+=dt;
	
      V.push_back(vm);
      W.push_back(w);
      T.push_back(time);
      G.push_back(br.g(w,vm));
      ions.push_back({br.i_K1(vm), br.i_x1(vm, w), br.i_Na(vm,w), br.i_s(vm,w)});
      
      if (mod(time, 50) < dt)
	{
	  std::cout << "Simulation at time: " << time
		    << "ms" << std::endl;
	}
    }
  

  std::ofstream foutv("./Results/BR/Sol.dat");
  std::ofstream foutg("./Results/BR/Derivatives.dat");
  std::ofstream fouti("./Results/BR/Current.dat");
  
  
  for (int k = 0; k < T.size(); k++)
    {
      foutv << T[k] << " " << V[k] << ", ";
      foutg << T[k] << " " << I[k] << " ";

      for (int i = 0; i < w.size(); i++)
	{
	  foutv << W[k][i] << ", ";
	  foutg << G[k][i] << " ";
	}
      foutv << std::endl;
      foutg << std::endl;

      fouti << T[k] << " ";
      for (int i = 0; i < 4; i++)
	fouti << ions[k][i] << " ";
      fouti << br.i_stim(T[k]);
      fouti << std::endl;
    }
  
  foutv.close();
  foutg.close();
  fouti.close();

  return 0;
}
