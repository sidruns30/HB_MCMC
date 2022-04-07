/*
Code that runs a smaller MCMC chain on the GAIA and color data of stars
A function (depending on radii, temperatures and masses of the stars) that
assumes black body distribution is used to compute the colors of the stars.
The color data is stored in ../data/magnitudes
MCMC function is borrowed from the MCMC wrapper2 file
*/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <omp.h>

#define MSUN 1.9885e33
#define RSUN 6.955e10
#define C 2.998e10
#define SQRT_2PI 2.5066282746
#define PI 3.14159265358979323846
#define NCHAINS 20
#define NPAST   100
#define ENABLE_OPENMP 1
#define MAGPARS 6		// Number of parameters
#define GAMMA (2.388/sqrt(2.*MAGPARS))
#define SUBN 4			// Numberof points in data
#define LEN(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))


struct bounds{
  double lo;
  double hi;
};

struct gauss_bounds{
	int flag;
};


typedef struct bounds bounds;
typedef struct gauss_bounds gauss_bounds;


/*
Function to get the temperature for a star from the log Mass (Msun). I am using tabulated
values given in the TESS portal paper. Note that the final temperature depends on an 
additional scaling parameter and boundary function. Returns temperature in log10K
*/
double _getT(double logM)
{
    // In solar masses
    double M_nodes[16] = {0.1, 0.26, 0.47, 0.59, 0.69, 0.87,
                          0.98, 1.085, 1.4, 1.65, 2.0, 2.5, 
                          3.0, 4.4, 15., 40.};
    // In log10 K
    double T_nodes[16] = {3.491, 3.531, 3.547, 3.584, 3.644, 3.712,
                          3.745, 3.774, 3.823, 3.863, 3.913, 3.991,
                          4.057, 4.182, 4.477, 4.623};

    double m = pow(10., logM);

    double T;

    // Edge cases
    if (m <= M_nodes[0])
    {
        T = T_nodes[0];
    }
    else if (m >= M_nodes[15])
    {
        T = T_nodes[15];

    }

    // Linear interp otherwise
    else
    {
        for (int j=0; j<16; j++)
        {
            if (m < M_nodes[j])
            {
                T = T_nodes[j-1] + (m - M_nodes[j-1]) * (T_nodes[j] - 
                                       T_nodes[j-1]) / (M_nodes[j] - M_nodes[j-1]);
                break;
            }
        }

    }

    return T;
}

/*
Function to get the radius of the star from its mass. Nodes taken from the TESS input catalog
paper and slightly tweaked by John Baker. Note that the final radius depends on an 
additional scaling parameter and boundary function. Returns radius in log10 Rsun
*/
double _getR(double logM)
{

    // In solar masses
    double M_nodes[10] = {0.07, 0.2, 0.356, 0.655, 0.784, 0.787, 1.377,
                          4.4, 15., 40.};
    // In log10 Rsun
    double logR_nodes[10] = {-0.953, -0.627, -0.423, -0.154, -0.082, -0.087,
                             0.295, 0.477, 0.792, 1.041};

        double m = pow(10., logM);

    // Edge cases
    if (m <= M_nodes[0])
        return logR_nodes[0];

    else if (m >= M_nodes[9])
        return logR_nodes[9];

    // Linear interp otherwise
    else
    {
        for (int j=0; j<10; j++)
        {
            if (m < M_nodes[j])
            {
                return logR_nodes[j-1] + (m - M_nodes[j-1]) * (logR_nodes[j] - 
                                       logR_nodes[j-1]) / (M_nodes[j] - M_nodes[j-1]);
            }
        }
    }
}

/*
Scaling functions for the temperature and the radius. Allows for flexibility in the mass and the
radius. Defined by John Baker; parameters tweaked by Siddhant to make the points lie within
1 sigma
*/
double envelope_Temp(double logM)
{
    /*The distribution of log10(x/y) ~ N(mu~0, std=0.02264)
    We assume that y(m) = model(m) x 10 ^ (scale x alpha)
    log10(x/y) is a normal distribution therefore we want scale x alpha
    to be a normal distribution. alpha is normally distributed around -1
    and 1 so we just rescale it by multiplying by the std of log10(x/y)
    */

    return 0.0224;
}

double envelope_Radius(double logM)
{
    double m = pow(10., logM);
    double n = 4.22;
    double slope = 15.68;
    double floor=0.01;
    double corner=1.055;
    double ceil=0.17;

    double boundary = 1/(1/ceil+1/(slope*pow((pow(m, n) + pow(corner, n)), (1/n)) - (slope*corner-floor)));

    return boundary;
}

// Standard gaussian
double gaussian(double x, double mean, double sigma)
{
  return (1 / sigma / SQRT_2PI) * exp(- pow((x - mean) / sigma, 2.));
}

// Get the log priors over all priors
double get_logP(double pars[], bounds limited[], bounds limits[], gauss_bounds gauss_pars[])
{
  double logP = 0.;
  double mean;
  double sigma;
  for (int i=0; i<MAGPARS; i++)
  {
    if (gauss_pars[i].flag == 1)
    {
      mean = 0.5 * (limits[i].lo + limits[i].hi);
      sigma = (limits[i].hi - limits[i].lo) / 2.;
      logP += log10(gaussian(pars[i], mean, sigma));
    }
  }
  return logP;
}


/*
GAIA magnitude model (Jeremy Schnittman 2021)
Assumes that the stars are two blackbodies, calculates the color spectrum from
the sum of the black body spectra
*/
void get_mags(double params[],  double D, double *model)
{
	double logM1 = params[0];
    double logM2 = params[1];

	double rr1 = params[2];
	double rr2 = params[3];

	double alpha_T1 = params[4];
	double alpha_T2 = params[5];

	double R1 = pow(10., _getR(logM1) + rr1*envelope_Radius(logM1)); 
	double R2 = pow(10., _getR(logM2) + rr2*envelope_Radius(logM2)); 

	double T1 = pow(10., _getT(logM1) + alpha_T1*envelope_Temp(logM1));
	double T2 = pow(10., _getT(logM2) + alpha_T2*envelope_Temp(logM2));
	
    // Convert Rsun to cgs
    R1 *= RSUN;
    R2 *= RSUN;

  double lam[4] = {442,540,673,750}; //wavelength in nm
  double nu[4],f_nu[4];
  double h = 6.626e-27;
  double k = 1.38e-16;
  double pc_cgs = 3.086e18;
  double blending = 0.;
  int j;


  for (j=0;j<4;j++) {
    nu[j]=C/(lam[j]*1e-7);
    f_nu[j]	=	PI*(SQR(R1)*(2.*h*CUBE(nu[j])/SQR(C)/expm1(h*nu[j]/(k*T1)))+
        		SQR(R2)*(2.*h*CUBE(nu[j])/SQR(C)/expm1(h*nu[j]/(k*T2))))
      			/(SQR(D)*SQR(pc_cgs));
    }

  double Bmag = -2.5*log10(f_nu[0])-48.6;
  double Vmag = -2.5*log10(f_nu[1])-48.6;
  double Gmag = -2.5*log10(f_nu[2])-48.6;
  double Tmag = -2.5*log10(f_nu[3])-48.6;

  // Return the differences
  double BminusV, VminusG, GminusT; 
  BminusV = Bmag - Vmag;
  VminusG = Vmag - Gmag;
  GminusT = Gmag - Tmag;

  model[0] = Gmag;
  model[1] = BminusV;
  model[2] = VminusG;
  model[3] = GminusT;
}

/*
Likelihood function
*/
double model_likelihood(double *input_data, double *data_err, 
						double *model, double *params, double Distance)
{
	get_mags(params, Distance, model);

	double chi2 = 0.;

	for (int i=0; i<SUBN; i++)
	{
		double residual = (input_data[i] - model[i]) / data_err[i];
		chi2 += residual * residual;
	}

	return (-chi2/2.0);
}

/* 
Function to check if (parameter) file exists
*/
int exists(const char *fname)
{
    FILE *file;
    if (access(fname, R_OK) == 0){
      return 1;
    }
    else {return 0;}
}

/* Functions to free memory*/
void free_1d(double *arr){
  free(arr);
}

void free_2d(double **arr, int size){
  int i;
  for (i=0;i<size;i++){free(arr[i]);}
  free(arr);
}

void free_3d(double ***arr, int size1, int size2){
  int i,j;
  for (i=0;i<size1;i++){
    for (j=0;j<size2;j++){free(arr[i][j]);}
    free(arr[i]);
  }
  free(arr);
}

/*Functions to use without openmp*/
double ran2(long *idum);
double gasdev2(long *idum);

/*Functions to use *with* openmp*/
double gsl_uniform(gsl_rng *r);
double gsl_gaussian(gsl_rng *r);

/*
Read the magnitude data and store it in an array
*/
void read_mag_data(char *TIC, double *ydata, double *yerr, double *distance)
{
	char fname[100] = "../data/magnitudes/";
	strcat(fname, TIC);
	strcat(fname, ".txt");

	FILE *data_file;
	int file_flag = 0;

	printf("Opening magnitude file %s \n", fname);

	if (exists(fname)) 
	{
		data_file = fopen(fname, "r");
		file_flag = 1;
	}
	else printf("Could not open magnitude file %s \n", fname); 

	if (file_flag == 1)
	{
		fscanf(data_file, "%lf\n", distance);
		for (int i=0; i<4; i++)
		{
			double tmp1, tmp2;
			fscanf(data_file, "%lf\t%lf\n", &tmp1, &tmp2);
			ydata[i] = tmp1;
			yerr[i] = tmp2;
		}
		fclose(data_file);
	}
}

/*
Set limits on the model parameters
*/
void set_limits(bounds limited[], bounds limits[], gauss_bounds gauss_pars[])
{
	//limits on M1, in log10 MSUN
	limited[0].lo = 1; 
	limits[0].lo = -1.5;
	limited[0].hi = 1;
	limits[0].hi = 2.0;
	gauss_pars[0].flag = 0;
	//limits on M2, in log10 MSUN
	limited[1].lo = 1; 
	limits[1].lo = -1.5;
	limited[1].hi = 1;
	limits[1].hi = 2.0;
	gauss_pars[1].flag = 0;
	//limits on rr1
	limited[2].lo = 1;
	limits[2].lo = -3.;
	limited[2].hi = 1;
	limits[2].hi = 3.0;
	gauss_pars[2].flag = 1;
	//limits on log rr2, the scale factor for R2
	limited[3].lo = 1;
	limits[3].lo = -3.;
	limited[3].hi = 1;
	limits[3].hi = 3.;
	gauss_pars[3].flag = 1;
	// limits on (log) Teff coefficient for star 1
	limited[4].lo = 1;
	limits[4].lo = -3.;
	limited[4].hi = 1;
	limits[4].hi = 3.;
	gauss_pars[4].flag = 1;
	// limits on (log) Teff coefficient for star 2
	limited[5].lo = 1;
	limits[5].lo = -3.;
	limited[5].hi = 1;
	limits[5].hi = 3.;
	gauss_pars[5].flag = 1;
}

/*
Create log Files for the MCMC
These files include the parameter file, chain file,
likelihood file and the output file
*/
void create_log_files(char *TIC, char* chain_fname, char* par_fname, 
					   char* logL_fname, char* out_fname,
					   FILE** chain_file, FILE** par_file, FILE** logL_file, FILE**
					   out_file)
{
	// Chain file
	strcat(chain_fname,"../data/chains/");
	strcat(chain_fname, TIC);
	strcat(chain_fname, "_GAIA_run.txt");
	*chain_file = fopen(chain_fname, "w");

	// Parameter file
	strcat(par_fname, "../data/subpars/");
	strcat(par_fname, TIC);
	strcat(par_fname, "_GAIA_run.txt");
	*par_file = fopen(par_fname, "w");

	// Likelihood file
	strcat(logL_fname, "../data/logL/");
	strcat(logL_fname, TIC);
	strcat(logL_fname, "_GAIA_run.txt");
	*logL_file = fopen(logL_fname, "w");;

	// Output file
	strcat(out_fname, "../data/GAIA_runs/");
	strcat(out_fname, TIC);
	strcat(out_fname, "_GAIA_run.txt");
	*out_file = fopen(out_fname, "w");
}

/*
Allocate memory for the arrays in the program
*/
void alloc_memory(double ***x, double ***P_,
				  double ****history, double **sigma, double **magdata,
				  double **magerr, double **xmap, int **index, 
				  double **logLx, int **DEtrial, double **model, double ***y)
{
	*sigma = (double *)malloc(MAGPARS*sizeof(double));

	*x = (double **)malloc(NCHAINS*sizeof(double *));
	for(int i=0;i<NCHAINS;i++) (*x)[i]=(double *)malloc(MAGPARS*sizeof(double));

	*y = (double **)malloc(NCHAINS*sizeof(double *));
	for(int i=0;i<NCHAINS;i++) (*y)[i]=(double *)malloc(MAGPARS*sizeof(double));

	*P_ = (double **)malloc(NCHAINS*sizeof(double *));
	for(int i=0;i<NCHAINS;i++) (*P_)[i]=(double *)malloc(MAGPARS*sizeof(double));

	*history = (double ***)malloc(NCHAINS*sizeof(double **));
  	for(int i=0;i<NCHAINS;i++) {
    	(*history)[i]=(double **)malloc(NPAST*sizeof(double *));
    	for(int j=0;j<NPAST;j++) (*history)[i][j]=(double *)malloc(MAGPARS*sizeof(double));
  	}	

	*magdata     = (double *)malloc(SUBN*sizeof(double));
    *magerr     = (double *)malloc(SUBN*sizeof(double));

	*xmap = (double *)malloc(MAGPARS*sizeof(double));
	*index = (int *)malloc(NCHAINS*sizeof(int));
	*logLx = (double *)malloc(NCHAINS*sizeof(double));
	*DEtrial = (int *)malloc(NCHAINS*sizeof(int));
	*model = (double *)malloc(SUBN*sizeof(double));
}

/*
Initialize Proposals - Set the sigma scale values for the various parameter distributions
*/
void init_proposals(double *sigma, double ***history)
{
	sigma[0] = 1.e-2;      // log M1 (MSUN)
	sigma[1] = 1.e-2;       // log M2 (MSUN)
	//sigma[0]  = 1.0e-2;  //log R1 (RSUN)
	//sigma[1]  = 1.0e-2;  //log R2 (RSUN)
	//sigma[2]  = 1.0e-3;  //log T1 (K)
	//sigma[3]  = 1.0e-3;  //log T2 (K)

}


/*
Initialize Chains -  parameters are picked randomly from the input parameter space
*/
void init_chain(bounds limits[], gsl_rng **r, double **x, double **P_, double *xmap,
				double *temp, int *index, double logLmap, double *magdata, 
				double *magerr, double distance, double *logLx, int NTHREADS)
{	
	for (int i=0; i<MAGPARS; i++)
	{
		for (int j=0; j<NCHAINS; j++)
		{
			int k = i % NTHREADS;
			double tmp = limits[i].lo + gsl_uniform(r[k])*(limits[i].hi - limits[i].lo);
			P_[j][i] = tmp;
			x[j][i]  = P_[j][i];
			xmap[i]  = P_[0][i];
		}
	}

	double dtemp = 1.2;
  
	temp[0] = 1.0;
	index[0] = 0;
  
	for(int i=1; i<NCHAINS; i++) 
	{
		temp[i]  = temp[i-1]*dtemp;
		index[i] = i;
	}

	double *model_ = (double *)malloc(SUBN*sizeof(double));
	
	logLmap = model_likelihood(magdata, magerr, model_, P_[0], distance);
	printf("initial chi2 %g\n",-2*logLmap);
	
	free(model_);
	for(int i=0; i<NCHAINS; i++) logLx[i] = logLmap;
}

void ptmcmc(int *index, double heat[], double logL[]);
void uniform_proposal(double *x, long *seed, bounds limits[], double *y);
void gaussian_proposal(double *x, gsl_rng *r, double *sigma, double scale, double temp, double *y);
void differential_evolution_proposal(double *x, gsl_rng *r, double **history, double *y);


void run_chain(gsl_rng *r, int iter, double **x, double *sigma, double *temp, 
				int *index, double ***history, int chain_id, int *DEtrial, bounds 
				limits[], bounds limited[], double *magdata, double *magerr, 
				double *model, double Distance, double *acc, double *DEacc, 
				double *logLx, double **y, gauss_bounds gauss_pars[])
{	
	double alpha = gsl_uniform(r);
	double jscale = pow(10., -6. + 6.*alpha);
	int jump = 0;
	if(gsl_uniform(r) < 0.5 && iter>NPAST) {jump=1;}
	//gaussian jumps along parameter directions

    if(jump==0){gaussian_proposal(x[index[chain_id]], r, sigma, jscale, 
									temp[chain_id], y[chain_id]);}
	if(jump==1) {
		if(index[chain_id]==0)DEtrial[chain_id] = 1;
		differential_evolution_proposal(x[index[chain_id]], r, history[chain_id], y[chain_id]);
		double dx_mag=0;
		for (int i=0;i<MAGPARS;i++) {
			dx_mag+=(x[index[chain_id]][i]-y[chain_id][i])*(x[index[chain_id]][i]-y[chain_id][i]);
		}
		if (dx_mag < 1e-6)
		gaussian_proposal(x[index[chain_id]], r, sigma, jscale, temp[chain_id], y[chain_id]);
      }

	// Enforce boundary conditions
	for (int i=0;i<MAGPARS;i++)
	{	//(reflecting boundary conditions)
		if ((limited[i].lo == 1)&&(y[chain_id][i] < limits[i].lo))
		y[chain_id][i] = 2.0*limits[i].lo - y[chain_id][i];
		if ((limited[i].hi == 1)&&(y[chain_id][i] > limits[i].hi)) 
		y[chain_id][i] = 2.0*limits[i].hi - y[chain_id][i];
		// (periodic boundary conditions)
		if ((limited[i].lo == 2)&&(y[chain_id][i] < limits[i].lo)) 
		y[chain_id][i] = limits[i].hi + (y[chain_id][i]-limits[i].lo);
		if ((limited[i].hi == 2)&&(y[chain_id][i] > limits[i].hi)) 
		y[chain_id][i] = limits[i].lo + (y[chain_id][i]-limits[i].hi);

	}

	double logPx, logPy;

    logPx = get_logP(x[index[chain_id]], limited, limits, gauss_pars);
    logPy = get_logP(y[chain_id], limited, limits, gauss_pars);

	logLx[index[chain_id]] = model_likelihood(magdata, magerr, model, x[index[chain_id]],
							 Distance);
	double logLy = model_likelihood(magdata, magerr, model, y[chain_id], Distance);
	/* evaluate new solution */
      //Hasting's ratio
      double H = exp( (logLy-logLx[index[chain_id]])/temp[chain_id]) * pow(10., logPy - logPx);
      //acceptance probability
      alpha = gsl_uniform(r);
      
	  //conditional acceptance of y
      if (alpha <= H) {
	      if(index[chain_id]==0) acc[chain_id] = 1;
	      for (int i=0;i<MAGPARS;i++) {
	        x[index[chain_id]][i]=y[chain_id][i];
	      }
        logLx[index[chain_id]] = logLy;
	      if(jump==1 && index[chain_id]==0) DEacc[chain_id] = 1;
      }


	  //printf("history assigned \n");
	  /* parallel tempering */
      //ptmcmc(index,temp,logLx);

	      /*  fill history  */
      //long k = iter - (iter/NPAST)*NPAST;
	  //for(int i=0; i<MAGPARS; i++) history[chain_id][k][i] = x[index[chain_id]][i];
}

/*
Write data to files
*/
void log_data(FILE** chain_file, FILE** par_file, FILE** logL_file, FILE** out_file,
				double **x, double *LogLx, double *magdata, double *model, int *index,
				double distance)
{
	// Write best chain parameters
	fprintf(*chain_file, "%.10g\t", LogLx[index[0]]);
	for (int i=0; i<MAGPARS; i++)
	{
		fprintf(*chain_file, "%.10g\t", x[index[0]][i]);
	}
	fprintf(*chain_file, "\n");

	// Write parameters
	rewind(*par_file);
	for (int i=0; i<MAGPARS; i++)
	{
		fprintf(*par_file, "%.10g\t", x[index[0]][i]);
	}

	fprintf(*par_file, "\n");

	// Write Likelihood
	for (int i=0; i<NCHAINS; i++)
	{
		fprintf(*logL_file, "%.10g\t", LogLx[index[i]]);
	}

	fprintf(*logL_file, "\n");

	// Write Magnitudes
	get_mags(x[index[0]], distance, model);

	rewind(*out_file);
	for (int i=0; i<SUBN; i++)
	{
		fprintf(*out_file, "%.10g\n", model[i]);
	}
}

/*
Close all the files
*/
void close_files(FILE** chain_file, FILE** par_file, FILE** logL_file, FILE** out_file)
{
	fclose(*chain_file);
	fclose(*logL_file);
	fclose(*par_file);
	fclose(*out_file);
}

/*
Data: Mag data
Model: get_mags(model)
Params: x
*/


/*
Main MCMC function
*/
void run_mcmc(long Niter, char TIC[20], int NTHREADS)
{
	char 	chain_fname[80] = "", par_fname[80] = "", 
			logL_fname[80] = "", out_fname[80] = "";
	// Arrays that store the data
	double **x;
	double ** P_;
	double ***history;
	double *sigma;

	double *model;
	double *magdata;
	double *magerr;
	double **y;

	double *xmap;
	double *logLx;
	double logLmap;

	long seed = Niter;

	double Distance;

	double heat[NCHAINS],temp[NCHAINS];
	
	int *index;
	int *DEtrial;
	double DEacc[NCHAINS], acc[NCHAINS];
	double DEacc_tot = 0., acc_tot = 0.;

	FILE* chain_file;
	FILE* par_file;
	FILE* logL_file;
	FILE* out_file;

	bounds limited[MAGPARS], limits[MAGPARS];
	gauss_bounds gauss_pars[MAGPARS];

    gsl_rng *r[NTHREADS];

    for (int i=0; i<NTHREADS; i++)
    {
        // Call returns pointer to the newly generated random
        // number type
        r[i] = gsl_rng_alloc (gsl_rng_ranlxs1);

        // Set the seed to each generator
        gsl_rng_set(r[i], seed + i);
    } 


	// Initialize everything
	create_log_files(TIC, chain_fname, par_fname, logL_fname, out_fname, 
						&chain_file, &par_file, &logL_file, &out_file);
	alloc_memory(&x, &P_, &history,  &sigma, &magdata, &magerr, &xmap, &index, &logLx,
				 &DEtrial, &model, &y);
	
	set_limits(limited, limits, gauss_pars);
	init_proposals(sigma, history);
	read_mag_data(TIC, magdata, magerr, &Distance);
	init_chain(limits, r, x, P_, xmap, temp, index, logLmap, magdata, magerr, 
				Distance, logLx, NTHREADS);
	int test_threads = 0;
	// Iteration loop
	for (int iter=0; iter<Niter; iter++)
	{
#pragma omp parallel for schedule(static)  //private(y)
		for (int chain_id = 0; chain_id<NCHAINS; chain_id++)
		{
			int k = chain_id % NTHREADS;
			test_threads = omp_get_num_threads();
			run_chain(r[k], iter, x, sigma, temp, 
				index, history, chain_id, DEtrial, limits,
				limited, magdata, magerr, model, 
				Distance, acc, DEacc, logLx, y, gauss_pars);
		}

		for (int chain_id = 0; chain_id<NCHAINS; chain_id++)
		{
			ptmcmc(index,temp,logLx);

			/*  fill history  */
			long k = iter - (iter/NPAST)*NPAST;
			for(int i=0; i<MAGPARS; i++) history[chain_id][k][i] = x[index[chain_id]][i];
		}

		// Store best parameters
		if (logLx[index[0]] > logLmap) {
      		for (int i=0;i<MAGPARS;i++) 
			{
	      		xmap[i]=x[index[0]][i];
      		}
      		logLmap = logLx[index[0]];
    	}

		// Write to files
		if (iter % 10 == 0)
		{
			log_data(&chain_file, &par_file, &logL_file, &out_file,
			x, logLx, magdata, model, index, Distance);
		}

		// Print results
		if (iter % 10000 == 0)
		{
			printf("Iter: %d \t Best likelihood %f \n", iter, logLx[index[0]]);
			printf("number of processes used = %d \n", test_threads);
		}
	}
	
	// Free the memory
	free_1d(magdata);
	free_1d(magerr);
	free_1d(sigma);
	free_1d(model);
	free_1d(logLx);
	//free_2d(y, NCHAINS);
	free(index);
	free(DEtrial);
	free_2d(x, NCHAINS);
	free_2d(P_, NCHAINS);
	free_3d(history, NCHAINS, NPAST);

	close_files(&chain_file, &par_file, &logL_file, &out_file);
	
	for (int i=0; i<NTHREADS; i++)
    {
        gsl_rng_free(r[i]);
    }
}

int main(int argc, char* argv[])
{	
	double Niter = (long)atoi(argv[1]);
	char TIC[20];
	strcpy(TIC,argv[2]);
	srand(Niter);
	int NTHREADS = 1;

	if (ENABLE_OPENMP)
  	{
    	NTHREADS = (long)atoi(argv[3]);
		omp_set_num_threads(NTHREADS);
		printf("number of threads = %d \n", NTHREADS);
  	}

	//for (int k=0; k<10; k++)
	//{
		//long int temp_ = 1727113733;//rand();
		//double res;
		//printf("Random number is %d \t", temp_);
		//res = ran2(&temp_);
		//printf(" value is %f \n", res);
	//}
	
	run_mcmc(Niter, TIC, NTHREADS);

	return 1;

}

/*************************Random Number Stuff*******************/


/***********************Series Functions ************************/
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define eps 1.2e-7
#define RNMX (1.0-eps)

double ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789; // some seed
	static long iy=0;
	static long iv[NTAB];
	double temp;
	
	printf("idum before: %d \n", *idum);
	//printf("idum2 start %d \n", idum2);
	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	//printf("idum2 after %d \n", idum2);
	//printf("idum after: %d \n", *idum);
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	//printf("interediate idum2: %d \n", idum2);
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	//printf("interediate idum2 (2): %d \n", idum2);
	printf("value of j is %d \n", j);
	iv[j] = *idum;
	//printf("idum2: %d \n", idum2);
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

// gaussian random number
double gasdev2(long *idum)
{
	double ran2(long *idum);
	static int iset=0;
	static double gset;
	double fac, rsq, v1, v2;
	
	if(*idum < 0) iset = 0;
	if(iset == 0){
		do{
			v1 = 2.0 * ran2(idum)-1.0;
			v2 = 2.0 * ran2(idum)-1.0;
			rsq = v1*v1+v2*v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset = v1 * fac;
		iset = 1;
		return(v2*fac);
	} else {
		iset = 0;
		return (gset);
	}
}

/***************************Parallel Functions ***********/
double gsl_uniform(gsl_rng *r)
{
	int random = gsl_rng_get(r);
	unsigned long int rmax = gsl_rng_max(r);
	double res = gsl_rng_uniform(r);
	if (res >= 1.)
	{
		printf("random number exceed 1 %f \n", res);
	}
	return res;
}

double gsl_gaussian(gsl_rng *r)
{
	double res =  gsl_ran_gaussian(r, 1.);
	return res;
}

/************************MCMC Iterations ****************/

void ptmcmc(int *index, double heat[], double logL[])
{
  int a, b;
	int olda, oldb;
	
	double heat1, heat2;
	double logL1, logL2;
	double dlogL;
	double H;
	double alpha;
	double beta;

	/*
	 Possible evidence for over-coupling using this
	 chain swapping scheme?  Neil & Laura & I see that
	 var(logL) < D/2 w/ PT, var(logL) ~ D/2 w/out.  
	 Neil just randomly chooses a pair of chains to 
	 exchange instead of marching down the line in his
	 EMRI code and doesn't see this problem.  But Joey
	 uses this in the eccentric binary code and also doesn't
	 see the problem.  WTF?  In the meantime, I'll switch to
	 randomly proposing a pair, but this is very puzzling.
	 */ 
	
  /* Siddhant: b can be -1, gives seg fault, put bounds on b*/
	b = (int)((double)rand()/(RAND_MAX)*((double)(NCHAINS-1)));
	a = b + 1;
	
	olda = index[a];
	oldb = index[b];
	heat1 = heat[a];
	heat2 = heat[b];
	logL1 = logL[olda];
	logL2 = logL[oldb];
	dlogL = logL2 - logL1;
	H  = (heat2 - heat1)/(heat2*heat1);
	alpha = exp(dlogL*H);
	beta  = (double)rand()/(RAND_MAX);
	if(alpha >= beta)
	{
		index[a] = oldb;
		index[b] = olda;
	}
}


void uniform_proposal(double *x, long *seed, bounds limits[], double *y)
{
	int n;
	double dx[MAGPARS];
	
	//compute size of jumps
	for(n=0; n<MAGPARS; n++) dx[n] = ran2(seed)*(limits[n].hi - limits[n].lo);
	
	//uniform draw on prior range
	for(n=0; n<MAGPARS; n++) y[n] = limits[n].lo + dx[n];
}

void gaussian_proposal(double *x, gsl_rng *r, double *sigma, double scale, double temp, double *y)
{
	int n;
	double gamma;
	double sqtemp;
	double dx[MAGPARS];
	//scale jumps by temperature
	sqtemp = sqrt(temp);
	//compute size of jumps
	for(n=0; n<MAGPARS; n++) dx[n] = gsl_gaussian(r)*sigma[n]*sqtemp*scale;
	//jump in parameter directions scaled by dx
	for(n=0; n<MAGPARS; n++) {
	  y[n] = x[n] + dx[n];
	}
}

void differential_evolution_proposal(double *x, gsl_rng *r, double **history, double *y)
{
	int n;
	int a;
	int b;
	double dx[MAGPARS];
	
	//choose two samples from chain history
	a = (int) (gsl_uniform(r)*NPAST);
	b = a;
	while(b==a) b = (int) (gsl_uniform(r)*NPAST);
	
	//compute vector connecting two samples
	for(n=0; n<MAGPARS; n++) {
	  dx[n] = history[b][n] - history[a][n];
	}
	//Blocks?
	
	//90% of jumps use Gaussian distribution for jump size
	if(gsl_uniform(r) < 0.9) for(n=0; n<MAGPARS; n++) dx[n] *= gsl_gaussian(r)*GAMMA;

	//jump along vector w/ appropriate scaling
	for(n=0; n<MAGPARS; n++) y[n] = x[n] + dx[n];
}
