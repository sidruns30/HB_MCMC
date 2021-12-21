# distutils: extra_compile_args = -O3 -std=c99

cdef extern from "likelihood3.c":
     void calc_light_curve (double* times, double Nt, double*pars, double *template)
     pass
