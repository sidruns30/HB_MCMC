cdef extern from "likelihood2.c":
     void calc_light_curve "calc_light_curve2" (double* times, double Nt, double*pars, double *template)
     pass
