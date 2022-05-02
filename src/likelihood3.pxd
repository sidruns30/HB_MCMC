# distutils: language=c
# distutils: extra_compile_args = -O3 -std=c99


cdef extern from "likelihood3.c":
     void calc_light_curve (double* times, double Nt, double*pars, double *template)
     void calc_radii_and_Teffs(double *params,  double *R1, double *R2, double *Teff1, double* Teff2)
     void calc_mags(double *params,  double D, double *Gmg, double *BminusV, 
double *VminusG, double *GminusT)
     double _getT(double)
     double _getR(double)
     double envelope_Radius(double logM)
     double envelope_Temp(double logM)
     pass
