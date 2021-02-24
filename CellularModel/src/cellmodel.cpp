/* 
   Classical models of cellular activity with parameters
   used in the 1D bidomain equation. 
   The function g(w,vm) describes the membrane dynamic 
   i.e. dw/dt = g(w,vm) where w is a collection of 
   gate-like variable (w=1 => open, w=0 => closed)
   representing the ionic channel on the cell's membrane
   and vm is the transmembrane potential vm = ui-ue
   where ui is the intracellular potential and ue the 
   extracellular potential.

   The number of gate variables p (w ∈ Rᵖ) depends on the model used.
*/

#ifndef CELLMODEL_CPP
#define CELLMODEL_CPP

#include "cellmodel.h"

/* 
   Mitchel and Schaeffer model:
   "A Two-Current Model for the Dynamics of Cardiac Membrane"
   COLLEEN C. MITCHELL AND DAVID G. SCHAEFFER
*/

MS::MS(){}

double MS::g(const double w, const double vm)
{
  if (vm < vgate)
    return (1.0/((vmax-vmin)*(vmax-vmin))-w)/tauopen;
  else
    return -w/tauclose;  
}

double MS::Iion(const double w, const double vm)
{
  return w*(vm-vmin)*(vm-vmin)*(vmax-vm)/(tauin*(vmax-vmin)) - (vm-vmin)/(tauout*(vmax-vmin));
}

double MS::i_stim(const double t)
{
  if (t > 50 and t < 50.5)
    return 1;
  else
    return 0.0;
}

double MS::Itot(const double w, const double vm, const double t)
{
  return Am*(Iion(w, vm) + i_stim(t)); 
}

/*
  Beeler Reuter model (1977):
  "RECONSTRUCTION OF THE ACTION POTENTIAL OF VENTRICULAR MYOCARDIAL FIBRES"
  G.W. BEELER and H. REUTER
*/


double BR::i_stim(const double t)
{
  if (t > 50 and t < 51.01)
    return 40;
  else
    return 0.0;
}

BR::BR(std::vector<double>& w, double& vm)
{
  InitModel(w, vm);
}

std::vector<double> BR::g(const std::vector<double>& w, const double vm)
{
  std::vector<double> g(7, 0.0); // intracellular calcium concentration +
   				 // 6 gating variables   

  g[0] =-1e-7*i_s(vm, w) + 0.07*(1e-7-w[0]); // [CA]

  for (int i = 0; i < 6; i++)
    {
      double a = alpha(vm, 2*i);
      double b = alpha(vm, 2*i+1);
      double tau = 1/(a+b);
      double y_inf = a/(a+b);    
      g[i+1] = (y_inf-w[i+1])/tau;
    }
  
  return g;
}

double BR::Iion(const std::vector<double>& w, const double vm)
{
  return i_K1(vm) + i_x1(vm, w) + i_Na(vm, w) + i_s(vm, w);
}

double BR::Itot(const std::vector<double>& w, const double vm,
		const double t)
{
  return -1./cm * (Iion(w, vm) - i_stim(t));
}

void BR::InitModel(std::vector<double>& w, double& vm)
{
  vm = -84;
  w = {1.86e-7, 0.1644, 0.01, 0.9814, 0.9673, 0.0033, 0.9884};
  
  // Init the parameters
  C.push_back({0.0005, 0.083, 50, 0, 0, 0.057, 1}); //a_x1
  C.push_back({0.0013, -0.06, 20, 0, 0, -0.04, 1}); //b_x1
  C.push_back({0, 0, 47, -1, 47, -0.1, -1}); // a_m
  C.push_back({40,-0.056, 72, 0,0,0,0}); // b_m
  C.push_back({0.126, -0.25, 77, 0,0,0,0}); // a_h
  C.push_back({1.7, 0, 22.5, 0, 0, -0.082, 1}); // b_h
  C.push_back({0.055, -0.25, 78, 0, 0, -0.2, 1}); // a_j
  C.push_back({0.3, 0, 32, 0,0, -0.1, 1}); // b_j
  C.push_back({0.095, -0.01, -5, 0, 0, -0.072, 1}); // a_d
  C.push_back({0.07, -0.017, 44, 0, 0, 0.05, 1}); // b_d
  C.push_back({0.012, -0.008, 28, 0, 0, 0.15, 1}); // a_f
  C.push_back({0.0065, -0.02, 30, 0, 0, -0.2, 1}); // b_f
    
}


// Ion currents
double BR::i_K1(const double vm)
{
  double ik1 = 0.35*( 4* (exp(0.04*(vm+85))-1)/(exp(0.08*(vm+53))+exp(0.04*(vm+53)))
  		      + 0.2*(vm+23)/(1-exp(-0.04*(vm+23))) );
  return ik1;
}

 double BR::i_x1(const double vm, const std::vector<double>& w)
{
  return w[1]*0.1973*(exp(0.04*(vm+77))-1)/exp(0.04*vm);
}
double BR::i_Na(const double vm, const std::vector<double>& w)
{
  return (gNa_bar*pow(w[2],3)*w[3]*w[4]+gNaC)*(vm-ENa);
}
double BR::i_s (const double vm, const std::vector<double>& w)
{
  double Es = -82.3 - 13.0287*log(w[0]);
  return gs*w[5]*w[6]*(vm-Es);
}

// Rate coefficients [ms⁻¹]
double BR::alpha(const double vm, int i)
{
  double a1 = C[i][0]*exp( C[i][1]*(vm+C[i][2]) );
  double a2 = C[i][3]*(vm+C[i][4]);
  double a3 = exp(C[i][5]*(vm+C[i][2]))+C[i][6];
  
  return (a1+a2)/a3;
}

/* Some useful operator overloadings and functions */
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());

  std::vector<double> result;
  result.reserve(a.size());

  std::transform(a.begin(), a.end(), b.begin(),
		 std::back_inserter(result), std::plus<double>());
  return result;
}
std::vector<double> operator*(const std::vector<double>& v, const double c)
{
  std::vector<double> result = v;
  for (unsigned int i = 0; i<result.size(); i++)
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


double mod(double a, double b)
{
  double mod;
  if (a < 0)
    mod = -a;
  else
    mod = a;
  if (b < 0)
    b = -b;

  while (mod >= b)
    mod = mod-b;
  if (a < 0)
    return -mod;
  return mod;
}

#endif
