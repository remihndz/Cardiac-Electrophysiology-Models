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

Vec MS::g(const Vec& w, const Vec& vm)
{
  int ndof = vm.size();
  Vec g(ndof);
  int k = 0;

  Vec above_thresh = -w.array()/tauclose;
  Vec under_thresh = (1.0/pow(vmax-vmin,2) - w.array())/tauopen;
  g = (vm.array() < vgate).select(under_thresh, above_thresh);

  // Use of iterator may make a loop faster?
  return g;
}

Vec MS::Iion(const Vec & w, const Vec & vm)
{
  int ndof = vm.size();
  Vec iion(ndof);
  int k = 0;

  iion = ( w.array()*(vm.array()-vmin)*(vm.array()-vmin)*(vmax-vm.array())
	   /(tauin*(vmax-vmin)) - (vm.array()-vmin)/(tauout*(vmax-vmin)) );
  return iion;
}

/*
  Beeler Reuter model (1977):
  "RECONSTRUCTION OF THE ACTION POTENTIAL OF VENTRICULAR MYOCARDIAL FIBRES"
  G.W. BEELER and H. REUTER
*/

BR::BR(MatW& w, Vec& vm)
{
  InitModel(w, vm);
}

MatW BR::g(const MatW & w, const Vec & vm)
{
  MatW g; // intracellular calcium concentration +
          // 6 current/gating variables   
  
  g.col(0) =-1e-7*(i_s(vm, w).array()) + 0.07*(1e-7-w.col(0).array()); // [CA]

  for (int i = 0; i < 6; i++)
    {
      Vec a = alpha(vm, 2*i);
      Vec b = alpha(vm, 2*i+1);
      Vec tau = 1/((a+b).array());
      Vec y_inf = a.array()/((a+b).array());    
      g.col(i+1) = (y_inf-w.col(i+1)).array()/tau.array();
    }
  
  return g;
}

Vec BR::Iion(const MatW & w, const Vec & vm)
{
  return i_K1(vm) + i_x1(vm, w) + i_Na(vm, w) + i_s(vm, w);
}

Vec BR::Itot(const MatW & w, const Vec & vm,
		const double t)
{
  return -1./cm * Iion(w, vm);
}

void BR::InitModel(MatW & w, Vec & vm)
{
  
  vm = Eigen::MatrixXd::Constant(vm.rows(), vm.cols(), -84.0);
  double w0[7]{1.86e-7, 0.1644, 0.01, 0.9814, 0.9673, 0.0033, 0.9884};
  for (int k = 0; k < 7; k++)
    w.col(k) = Eigen::MatrixXd::Constant(vm.rows(), 1, w0[k]);
  
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
Vec BR::i_K1(const Vec & vm)
{
  Vec ik1 = vm;
  ik1.unaryExpr([](double x) {return 0.35*( 4* (exp(0.04*(x+85))-1)/
					    (exp(0.08*(x+53))+exp(0.04*(x+53)))
					    + 0.2*(x+23)/(1-exp(-0.04*(x+23))) );});
  return ik1;
}

Vec BR::i_x1(const Vec & vm, const MatW & w)
{
  Vec ix1 = vm;
  ix1.unaryExpr([](double x) {return 0.1973*(exp(0.04*(x+77))-1)/exp(0.04*x);});
  return w.col(0).array() * ix1.array();
}

Vec BR::i_Na(const Vec & vm, const MatW & w)
{
  return( (gNa_bar* (w.col(2).array()).pow(3) *
	   (w.col(3).array()*w.col(4).array()+gNaC))
	   *(vm.array()-ENa) );
}
Vec BR::i_s (const Vec & vm, const MatW & w)
{
  Vec Es = -82.3 - 13.0287*(w.col(0).array()).log();
  return gs*(w.col(5)).array() * (w.col(6)).array() *(vm.array()-Es.array());
}

// Rate coefficients [ms⁻¹]
Vec BR::alpha(const Vec & vm, int i)
{
  Vec a1 = C[i][0]*exp( C[i][1]*(vm.array()+C[i][2]) );
  Vec a2 = C[i][3]*(vm.array()+C[i][4]);
  Vec a3 = exp(C[i][5]*(vm.array()+C[i][2]))+C[i][6];
  
  return (a1+a2).array()/a3.array();
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
