#ifndef LKL3_H

  #define LKL3_H 123456789
  #define PI 3.14159265358979323846
  #define G 6.6743e-8 //cgs
  #define C 2.998e10
  #define AU 1.496e13
  #define MSUN 1.9885e33
  #define RSUN 6.955e10
  #define SEC_DAY 86400.0
  #define USE_COLOR_INFO 0
  #define STORE_DATA 1
  #define ALPHA_FREE 1 // to set coefficitents as parameters in the model
  #define ALPHA_MORE 1 // to add even more flexible coefficients
  #define BLENDING 1 // Blending fraction (currently only works when alpha free and alpha more are enabled)
  #if ALPHA_FREE == 1
    #if ALPHA_MORE == 1
      #if BLENDING == 1
          #define NPARS 22
      #else
          #define NPARS 20
      #endif
    #else
      #define NPARS 16
    #endif
  #else
      #define NPARS 10
  #endif
  #define SAVECOMP 0 // Save the individual lightcurve components in the txt file
  #define BIG_NUM 1.e15

  #define SQR(x) ((x)*(x))
  #define CUBE(x) ((x)*(x)*(x))
  #define QUAD(x) ((x)*(x)*(x)*(x))

#endif // LKL3_H

struct bounds_old
{
  double lo;
  double hi;
};

struct gauss_bounds_old
{
  int flag;
};

#define NTAB 32
// structure to hold random number values
struct RNG_Vars_old
{
  long idum2; // some seed
  long iy;
  long iv[NTAB];
  int iset;
  double gset;
  // Number of times the rng is called
  long cts;
};


typedef struct bounds_old bounds_old;
typedef struct gauss_bounds_old gauss_bounds_old;
typedef struct RNG_Vars_old RNG_Vars_old;

double partition_old(double arr[], int low, int high);
void quickSort_old(double arr[], int low, int high);
void remove_median_old(double *arr, long begin, long end);
void traj_old(double t, double pars[], int i, double *X1, double *X2, double *Y1, double *Y2,
        double *Z1, double *Z2, double *rr, double *ff);
double get_alpha_beam_old(double logT);
double beaming_old(double P, double M1, double M2, double e, double inc,
                double omega0, double nu, double alpha_beam);
double ellipsoidal_old(double P, double M1, double M2, double e, double inc,
                double omega0, double nu, double R1, double a, double mu, double tau);
double reflection_old(double P, double M1, double M2, double e, double inc, 
                    double omega0, double nu , double R2, double alpha_ref1);
double eclipse_area_old(double R1, double R2, 
                double X1, double X2, double Y1, double Y2);
void calc_mags_old(double params[],  double D, double *Gmg, double *BminusV, 
double *VminusG, double *GminusT);
void calc_light_curve_old(double *times, long Nt, double *pars,  double *template);
int RocheOverflow_old(double *pars);
double loglikelihood_old(double time[], double lightcurve[], double noise[],
		     long N, double params[], double mag_data[], double magerr[], 
          int weight);
void set_limits_old(bounds_old limited[], bounds_old limits[], gauss_bounds_old gauss_pars[], double LC_PERIOD);