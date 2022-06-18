#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "likelihood3.h"

#define SQRT_2PI 2.5066282746
#define NCHAINS 50
#define NPAST   500
#define GAMMA (2.388/sqrt(2.*NPARS))
#define NUM_ELEMENTS(x) (size_of(x)/size_of(x[0]))
#define USE_RAND_PARS 1
#define ENABLE_OPENMP 1

/* From likelihood.c */
void calc_light_curve(double t_data[], long Nt, double P_[], double light_curve[]);
/* MCMC functions*/
void ptmcmc(int *index, double temp[], double logL[], double logP[], FILE *temp_swap_file);
double get_logP(double pars[], bounds limited[], bounds limits[], gauss_bounds gauss_pars[]);
/* Proposal distributions for MCMC */
double ran2(long *idum);
double gasdev2(long *idum);
void uniform_proposal(double *x, long *seed, bounds limits[], double *y);
void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y);
void hypercube_proposal(double *x, long *seed, double *y);
void differential_evolution_proposal(double *x, long *seed, double **history, double *y);
/* Same functions but for multithreading */
double ran2_parallel(long *idum, RNG_Vars *state);
double gasdev2_parallel(long *idum, RNG_Vars *state);
void uniform_proposal_parallel(double *x, long *seed, bounds limits[], double *y, RNG_Vars *state);
void gaussian_proposal_parallel(double *x, long *seed, double *sigma, double scale, double temp, double *y, RNG_Vars *state);
void differential_evolution_proposal_parallel(double *x, long *seed, double **history, double *y, RNG_Vars *state);
/* Other helpful functions */
void free_1d(double *arr);
void free_2d(double **arr, int size);
void free_3d(double ***arr, int size1, int size2);
double gaussian(double x, double mean, double sigma);
void TEST_lc_calc(double *pars);
int exists(const char *fname);
void LogSuspiciousJumps(FILE* logfile, long iter, int chain_id, double H, double alpha, double chain_temp, double logL_old, double logL_new,
                        double logP_old, double logP_new, double *pars_old, double *pars_new, int jump_type);