/* 
   Classical models of cellular activity with parameters
   used in the simulation of myocardial fibers. 
   The function g(w,vm) describes the cell membrane 
   dynamics i.e. dw/dt = g(w,vm) where w is a collection
   of currents and dimensionless gating variables
   acting on the transmembrane potential vm = ui-ue

   The number of  variables p (w ∈ Rᵖ) depends on the 
   model used (1 for MS, 7 for BR).
*/

#ifndef CELLMODEL_H
#define CELLMODEL_H

#include <vector>
#include <eigen3/Eigen/Dense>

typedef Eigen::Matrix<double, Eigen::Dynamic, 7> MatW;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vec;

double mod(double a, double b);

Vec i_stim(const double t);

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
  MS(Vec& w, Vec& vm);
  Vec g(const Vec& w, const Vec& vm);
  Vec Iion(const Vec& w, const Vec& vm);
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
  Vec alpha(const Vec& vm, int i);

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
  Vec i_K1(const Vec& vm);
  Vec i_x1(const Vec& vm, const MatW& w);
  Vec i_Na(const Vec& vm, const MatW& w);
  Vec i_s (const Vec& vm, const MatW& w);

  BR();
  BR(MatW & w, Vec& vm);
  MatW g(const MatW & w, const Vec& vm);
  Vec  Iion(const MatW & w, const Vec& vm);
  Vec  Itot(const MatW & w, const Vec& vm,
	    const double t);
  void InitModel(MatW & w, Vec& vm);
};
  
#endif
