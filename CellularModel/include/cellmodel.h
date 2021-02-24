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

#ifndef CELLMODEL_H
#define CELLMODEL_H

#include <math.h>
#include <vector>
#include <assert.h>
#include <algorithm>

/* Some useful operator overloadings and functions */

std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);
std::vector<double> operator*(const std::vector<double>& v, const double c);
std::vector<double> operator*(const double c, const std::vector<double>& v);

double mod(double a, double b);

double i_stim(const double t);


/* Cellular models */

class MS
{
 private:
  double vgate    = -67.0;
  double tauclose = 100.0;
  double tauopen  = 300.0;
  double tauout   = 90.0;
  double tauin    = 4.0;

 public:
  double vmin     = -80.0;
  double vmax     = 20.0;
  double Am       = 200.0;

  /* 
     w is a vector with pointwise value of the variable
     i.e. w[k] = w(x_k) where w(x_k) is the value of 
     w at point x_k ∈ Ω 
  */
  
  MS();
  double g(const double w, const double vm);
  double i_stim(const double t);
  double Iion(const double w, const double vm);
  double Itot(const double w, const double vm, const double t);
};


class BR
{
 private:
  // Membrane capacitance
  double cm = 1; // [μF/cm²]

  // Parameters
  // Conductances [mmho/cm²]
  double gNa_bar = 4, gNaC = 0.003, gs = 0.09;

  // Voltages [mV]
  double ENa = 50;

  // Rate coefficients [ms⁻¹]
  double alpha(const double vm, int i);

  // Coefficients for rates [1/ms, 1/mV, 1/(mV·ms), mV, 1/mV, 1]
  std::vector<std::vector<double>> C;

 public:
  /* 
     w is a matrix with pointwise value of the
     ionic variable l in the row l
     i.e. w[l,k] = w_l(x_k) where w_l(x_k) is
     the value of w_l at point x_k ∈ Ω 
     k = 0...7 
  */

  // Ion currents
  double i_stim(const double t);
  double i_K1(const double vm);
  double i_x1(const double vm, const std::vector<double>& w);
  double i_Na(const double vm, const std::vector<double>& w);
  double i_s (const double vm, const std::vector<double>& w);

  BR(std::vector<double>& w, double& vm);
  std::vector<double> g(const std::vector<double>& w, const double vm);
  double              Iion(const std::vector<double>& w, const double vm);
  double              Itot(const std::vector<double>& w, const double vm,
			   const double t);
  void InitModel(std::vector<double>& w, double& vm);
};
  
#endif
