# distutils: language=c
# distutils: extra_compile_args = -O3 -std=c99


cdef extern from "likelihood3.c":
     void calc_light_curve (double* times, double Nt, double*pars, double *template)
     void calc_radii_and_Teffs(double *params,  double *R1, double *R2, double *Teff1, double* Teff2)
     pass
