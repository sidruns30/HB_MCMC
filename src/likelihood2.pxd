# distutils: extra_compile_args = -O3 -std=c99

cdef extern from "likelihood2.c":
     void calc_light_curve "calc_light_curve2" (double* times, double Nt, double*pars, double *template)
     pass
