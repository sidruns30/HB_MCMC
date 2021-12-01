cdef extern from "likelihood2.c":
     void calc_light_curve (double* times, double Nt, double*pars, double *template)
     pass
