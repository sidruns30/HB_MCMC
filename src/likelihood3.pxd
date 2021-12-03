cdef extern from "likelihood3.c":
     void calc_light_curve (double* times, double Nt, double*pars, double *template)
     pass
