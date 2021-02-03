The bidomain equations are frequently used to model electrical activity within the heart. They read:
<p align="center">
<img src="https://latex.codecogs.com/png.latex?%5Cinline%20%5Cdpi%7B150%7D%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20%5Cnabla%5Ccdot%28%5Csigma_i%5Ccdot%5Cnabla%20%5Cphi_e%29%20&plus;%20%5Cnabla%5Ccdot%20%28%5Csigma_i%5Ccdot%5Cnabla%20V_m%29%20%3D%20A_m%28C_m%5Cfrac%7B%5Cpartial%20V-m%7D%7B%5Cpartial%20t%7D%20&plus;%20I_%7Bion%7D%20-%20I_%7Bapp%7D%29%5C%5C%20%5C%5C%20%5Cnabla%5Ccdot%5Cbigl%28%28%5Csigma_i&plus;%5Csigma_e%29%5Cnabla%5Cphi_e%5Cbigr%29%20&plus;%20%5Cnabla%5Ccdot%28%5Csigma_i%5Ccdot%5Cnabla%20V_m%29%20%3D%200%5C%5C%20%5C%5C%20%5Cfrac%7B%5Cpartial%20w%7D%7B%5Cpartial%20t%7D%20%3D%20g%28V_m%2C%20w%29%20%5Cend%7Bmatrix%7D%5Cright.">
</p>

The last of those equations simulates the cell dynamics (excitable->action potential->refractory period->excitable). There are many choices for this set of ODEs, the simplet probably being the Mitchel and Schaeffer ([Mitchel and Schaeffer's paper](https://doi.org/10.1016/S0092-8240(03)00041-7)).

# What is in this repository
Mainly, some C++ codes to solve the bidomain equations using different methods (finite differences, finite volumes, finite elements?,...).
