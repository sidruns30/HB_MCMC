/*********************************************************/
//similar to mcmc_wrapper.c, but uses likelihood3.c,
//with independently varying R1,R2
//instead of using main sequence scaling relations
/*********************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include<unistd.h>
#include "likelihood3.h"

#define SQRT_2PI 2.5066282746
#define NCHAINS 50
#define NPAST   500
#define GAMMA (2.388/sqrt(2.*NPARS))
#define NUM_ELEMENTS(x) (size_of(x)/size_of(x[0]))
#define ENABLE_OPENMP 1

#if ENABLE_OPENMP == 1
  #include <omp.h>
#endif

/* From likelihood.c */
void calc_light_curve(double t_data[], long Nt, double P_[], double light_curve[]);
/* Siddhant: What are ran2 and gasdev2 */ 
void ptmcmc(int *index, double heat[], double logL[]);
double ran2(long *idum);
double gasdev2(long *idum);
/* Proposal distributions for MCMC */
void initialize_proposals(double *sigma, double ***history);
void uniform_proposal(double *x, long *seed, bounds limits[], double *y);
void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y);
void hypercube_proposal(double *x, long *seed, double *y);
void differential_evolution_proposal(double *x, long *seed, double **history, double *y);


// CHeck diff runs
// Adds John's version
// Fold on lc
// Profiler

double SIGMAP, SIGMA_BLEND;

/* Siddhant: Functions to free memory from arrays */
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

// Function to check if (parameter) file exists
int exists(const char *fname){
    FILE *file;
    if (access(fname, R_OK) == 0){
      return 1;
    }
    else {return 0;}
}


/* Set priors on parameters, and whether or not each parameter is bounded*/
/* Siddhant: Maybe just use an if condition/switch statment instead of limited*/
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
  //limits on P, in log10 days
  limited[2].lo = 1;
  limits[2].lo = -2.0;
  limited[2].hi = 1;
  limits[2].hi = 3.0;
  gauss_pars[2].flag = 0;
  //limits on e
  limited[3].lo = 1;
  limits[3].lo = 0.0;
  limited[3].hi = 1;
  limits[3].hi = 1;
  gauss_pars[3].flag = 0;
  //limits on inc, in rads
  limited[4].lo = 2;
  limits[4].lo = 0;
  limited[4].hi = 2;
  limits[4].hi = PI/2;
  gauss_pars[4].flag = 0;
  //limits on Omega, in rads
  limited[5].lo = 2;
  limits[5].lo = -PI;
  limited[5].hi = 2;
  limits[5].hi = PI;
  gauss_pars[5].flag = 0;
  //limits on omega0, in rads
  limited[6].lo = 2;
  limits[6].lo = -PI;
  limited[6].hi = 2;
  limits[6].hi = PI;
  gauss_pars[6].flag = 0;
  //limits on T0, in MDJ-2450000
  limited[7].lo = 1;
  limits[7].lo = -1000;
  limited[7].hi = 1;
  limits[7].hi = 1000.;
  gauss_pars[7].flag = 0;
  //limits on log rr1, the scale factor for R1
  limited[8].lo = 1;
  limits[8].lo = -3.;
  limited[8].hi = 1;
  limits[8].hi = 3.;
  gauss_pars[8].flag = 1.;
  //limits on log rr2, the scale factor for R2
  limited[9].lo = 1;
  limits[9].lo = -3.;
  limited[9].hi = 1;
  limits[9].hi = 3.;
  gauss_pars[9].flag = 1.;
  if (ALPHA_FREE == 1){
    // Limits of the alpha_coefficients
    // limits on limb darkening coefficient for star 1
    limited[10].lo = 1;
    limits[10].lo = 0.12;
    limited[10].hi = 1;
    limits[10].hi = 0.20;
    gauss_pars[10].flag = 1;
    // limits on gravity darkening coefficient for star 1
    limited[11].lo = 1;
    limits[11].lo = 0.3;
    limited[11].hi = 1;
    limits[11].hi = 0.38;
    gauss_pars[11].flag = 1;
    // limits on limb darkening coefficient for star 2
    limited[12].lo = 1;
    limits[12].lo = 0.12;
    limited[12].hi = 1;
    limits[12].hi = 0.20;
    gauss_pars[12].flag = 1;
    // limits on gravity darkening coefficient for star 2
    limited[13].lo = 1;
    limits[13].lo = 0.3;
    limited[13].hi = 1;
    limits[13].hi = 0.38;
    gauss_pars[13].flag = 1;
    // limits on reflection coefficients on star 1
    limited[14].lo = 1;
    limits[14].lo = 0.8;
    limited[14].hi = 1;
    limits[14].hi = 1.2;
    gauss_pars[14].flag = 1;
    // limits on reflection coefficients on star 2
    limited[15].lo = 1;
    limits[15].lo = 0.8;
    limited[15].hi = 1;
    limits[15].hi = 1.2;
    gauss_pars[15].flag = 1;
    if (ALPHA_MORE == 1){
      // limits on extra (log) beaming coefficient for star 1
      limited[16].lo = 1;
      limits[16].lo = -0.1;
      limited[16].hi = 1;
      limits[16].hi = 0.1;
      gauss_pars[16].flag = 1;
      // limits on extra (log) beaming coefficient for star 2
      limited[17].lo = 1;
      limits[17].lo = -0.1;
      limited[17].hi = 1;
      limits[17].hi = 0.1;
      gauss_pars[17].flag = 1;
      // limits on Teff coefficient for star 1
      limited[18].lo = 1;
      limits[18].lo = -3.;
      limited[18].hi = 1;
      limits[18].hi = 3.;
      gauss_pars[18].flag = 1.;
      // limits on Teff coefficient for star 2
      limited[19].lo = 1;
      limits[19].lo = -3.;
      limited[19].hi = 1;
      limits[19].hi = 3.;
      gauss_pars[19].flag = 1.;
      if (BLENDING == 1){
        // Blending coefficient in the flux
        limited[20].lo = 1;
        limits[20].lo = 0.;
        limited[20].hi = 1;
        limits[20].hi = 1.;
        gauss_pars[20].flag = 0;
        // FLux tune coefficient
        limited[21].lo = 1;
        limits[21].lo = 0.99;
        limited[21].hi = 1;
        limits[21].hi = 1.01;
        gauss_pars[21].flag = 0;
      }
    }
  }
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
  for (int i=0; i<NPARS; i++)
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



int main(int
 argc, char* argv[])
{

  /* Siddhant: Adding small variable descriptions */
  long Niter;                     // Chain iterations
  double **P_;                    // Contains parameters - check why different from x - seems uneccessary
  double P_0[NPARS];            
  double **x;                     // Parameter chains
  double ***history;              // Chain history
  double xmap[NPARS];           // Contains the chain with best paramters
  double dx[NPARS];             // Parameter step
  double dx_mag;                  
  double logLx[NCHAINS];          // Log likelihood for all chains
  double logLmap;                 
  double logLy;                   // Likelihood for new step
  double tmp,tmp1,tmp2,tmp3,tmp4;
  double true_err,Nerr,suberr;    // True error contains error on each LC point; others not used
  double scale;                   
  double dtemp;                   // Ask Jeremy
  double heat[NCHAINS],temp[NCHAINS]; // More PTMCMC stuff
  double *rdata;                  
  double *a_data,*a_model,*e_data,*t_data,*q_data;
  double *mag_data, *mag_err;
  double *light_curve;            
  double *light_curve2;
  double blending;        
  long i,j,k,m,subi;              
  long Nt,Nstart,Nstop,subN,obs_id,nch;
  long acc;                       
  long iter,TICid;                
  int *index, burn_in;
  int run;
  int nthreads;
  double ls_period;  
  //double SIGMAP;        

  bounds limited[NPARS], limits[NPARS];
  gauss_bounds gauss_pars[NPARS];
  FILE *param_file, *data_file, *chain_file, *logL_file, *mag_file;
  /* Siddhant: Initializing empty arrays help avoid memory leaks*/
  char  pfname[80]="", dfname[80]="",  parname[80]="",
        subparname[80]="",  outname[80]="", chainname[80]="",
        logLname[80]="", RUN_ID[11]="", mag_name[80]="";
  //characteristic magnitude of MCMC step in each parameter
  double  *sigma;

  Niter = (long)atoi(argv[1]);
  strcpy(RUN_ID,argv[2]);
  true_err = (double)atof(argv[3]);
  burn_in = atoi(argv[4]);
  ls_period = (double)atof(argv[5]);
  run = atoi(argv[6]);
  blending = (double)atof(argv[7]);
  nthreads = (long)atoi(argv[8]);

  if (ENABLE_OPENMP)
  {
    omp_set_num_threads(nthreads);
  }
  //double weight = (double)atof(argv[7]);
  // Set period search
  if (burn_in == 2) 
  {
    SIGMAP = 1.e-10;
  }
  else              {
    SIGMAP = 1.e-1;
    }

  strcat(subparname,"../data/subpars/subpar.");
  strcat(subparname,RUN_ID);
  strcat(subparname,".dat");
  strcat(parname,"../data/pars/par.");
  strcat(parname,RUN_ID);
  strcat(parname,".dat");
  strcat(chainname,"../data/chains/chain.");
  strcat(chainname,RUN_ID);
  strcat(chainname,".dat");
  strcat(logLname,"../data/logL/logL.");
  strcat(logLname,RUN_ID);
  strcat(logLname,".dat");
  strcat(outname,"../data/lightcurves/mcmc_lightcurves/");
  strcat(outname,RUN_ID);
  strcat(outname,".out");
  strcpy(mag_name,"../data/magnitudes/");
  strcat(mag_name,RUN_ID);
  strcat(mag_name,".txt");
  printf("%s\n",parname);

  //double weight = (double)atof(argv[7]);
  char suffix[30];
  sprintf(suffix, "%d", run);
  strcat(subparname, suffix);
  strcat(parname, suffix);
  strcat(chainname, suffix);
  strcat(logLname, suffix);
  strcat(outname, suffix);

  printf("%s\n", subparname);
  printf("%s\n", parname);
  printf("%s\n", chainname);
  printf("%s\n", logLname);
  printf("%s\n", outname);

  /*Use binned lc if period is found, otherwise use the full lightcurve*/
  if (burn_in == 2){
  strcpy(dfname,"../data/lightcurves/folded_lightcurves/");
  strcat(dfname,RUN_ID);
  strcat(dfname,"_new.txt");
  }
  else {
    strcpy(dfname,"../data/lightcurves/original/");
    strcat(dfname,RUN_ID);
    strcat(dfname,".txt");
  }


  P_ = (double **)malloc(NCHAINS*sizeof(double));
  for(i=0;i<NCHAINS;i++) P_[i]=(double *)malloc(NPARS*sizeof(double));
  
  x = (double **)malloc(NCHAINS*sizeof(double));
  for(i=0;i<NCHAINS;i++) x[i]=(double *)malloc(NPARS*sizeof(double));


  //y = (double *)malloc(NPARS*sizeof(double));
  
  history = (double ***)malloc(NCHAINS*sizeof(double));
  for(i=0;i<NCHAINS;i++) {
    history[i]=(double **)malloc(NPAST*sizeof(double));
    for(j=0;j<NPAST;j++)history[i][j]=(double *)malloc(NPARS*sizeof(double));
  }

  sigma = (double *)malloc(NPARS*sizeof(double));
  scale  = 1.0/((double)(NPARS));
  //

  //true_err = 1.0;
  srand(Niter);
  long seed = run;
  long chain_seeds[NCHAINS];

  for (int i=0; i<NCHAINS; i++)
  {
    chain_seeds[i] = run + i;
  }

  //read in starting parameters near best solution
  printf("Opening parameter file \n");
  if (exists(parname)){param_file = fopen(parname,"r");}
  printf("Parameter file: %s\n",parname);
  //initialize priors
  set_limits(limited,limits, gauss_pars);
  //set up proposal distribution
  initialize_proposals(sigma, history);
  sigma[2] = SIGMAP;

  int param_file_flag = 0;
  /* Siddhant: Initializing the chain parameters*/
  for (i=0;i<NPARS;i++) {
    if (exists(parname)){
      fscanf(param_file,"%lf\n", &tmp); 
      printf("Par %d value: %.5e \n",i, tmp);
      param_file_flag = 1;
      }
  
    else {

      tmp = limits[i].lo + ran2(&seed)*(limits[i].hi - limits[i].lo);
      printf("\t assinging random pars \n");
      }

    P_[0][i] = tmp;
    P_0[i]   = tmp;
    x[0][i]  = P_[0][i];
    xmap[i]  = P_[0][i];
    for(j=1; j<NCHAINS; j++) {
      P_[j][i] = P_[0][i];
      x[j][i]  = P_[0][i];
    }
    // burn_in == 1 seems redundant?
    if (burn_in == 1) {
      for(j=0; j<NCHAINS; j++) {
        if (param_file_flag != 1){
          P_[j][i] = limits[i].lo + ran2(&seed)*(limits[i].hi - limits[i].lo);
          x[j][i]  = P_[j][i];
        }
        
      }
    }

    if (burn_in == 2) {
      for(j=0; j<NCHAINS; j++) {
        //if (param_file_flag != 1){
          P_[j][i] = limits[i].lo + ran2(&seed)*(limits[i].hi - limits[i].lo);
          if (i == 2) P_[j][i] = ls_period;  // Keep the period
          if (i == 20) P_[j][i] = blending;
          x[j][i]  = P_[j][i];
        // }
      }
    }
  }

  if (burn_in == 1){
      //printf("%s", dfname);
    if (exists(parname)){fclose(param_file);}
    /* Read binned data file */
    printf("Opening data file %s", dfname);
    data_file = fopen(dfname,"r");
    tmp = 0;
    fscanf(data_file,"%ld\n", &Nt);
    subN         = 0;
    rdata        = (double *)malloc(4*Nt*sizeof(double));
    index        = (int *)malloc(NCHAINS*sizeof(int));

    for (i=0;i<Nt;i++) {
      fscanf(data_file,"%lf %lf %lf %lf\n", &tmp1, &tmp2, &tmp3, &tmp4);
      rdata[i*4]=tmp1;
      rdata[i*4+1]=tmp2;
      rdata[i*4+2]=tmp3;
      rdata[i*4+3]=tmp4;
      if (tmp4 == 0) subN++;
    }
    printf("Closing data file \n");
    fclose(data_file);
    
    a_data       = (double *)malloc(subN*sizeof(double));
    a_model      = (double *)malloc(subN*sizeof(double));
    t_data       = (double *)malloc(subN*sizeof(double));
    e_data       = (double *)malloc(subN*sizeof(double));
    q_data       = (double *)malloc(subN*sizeof(double));
    light_curve  = (double *)malloc(subN*sizeof(double));
    light_curve2 = (double *)malloc(subN*sizeof(double));
    mag_data     = (double *)malloc(5*sizeof(double));
    mag_err     = (double *)malloc(4*sizeof(double));
    
    subi = 0;
    for (i=0;i<Nt;i++) {
      tmp4 = rdata[i*4+3];
      if (tmp4 == 0) {
      t_data[subi] = rdata[i*4]; 
      a_data[subi] = rdata[i*4+1]; 
      e_data[subi] = sqrt(rdata[i*4+1])*true_err;
      subi++;
      }
  }
    printf("Opening magnitude file");
    mag_file = fopen(mag_name, "r");
    double tmp1, tmp2;
    fscanf(mag_file, "%lf\n", &tmp1);
    mag_data[0] = tmp1;
    for (i=0;i<4;i++){
      fscanf(mag_file, "%lf\t%lf\n", &tmp1, &tmp2);
      mag_data[i+1] = tmp1;
      mag_err[i] = tmp2;
    }
    fclose(mag_file);      
  }

  else if (burn_in == 2){
    if (exists(parname)){fclose(param_file);}
    printf("Opening folded lc data file %s \n", dfname);
    /* Read binned data file */
    data_file = fopen(dfname,"r");
    tmp = 0;
    fscanf(data_file,"%ld\n", &Nt);
    Nt = Nt;
    subN         = 0;
    rdata        = (double *)malloc(3*Nt*sizeof(double));
    index        = (int *)malloc(NCHAINS*sizeof(int));

    for (i=0;i<Nt;i++) {
      fscanf(data_file,"%lf\t%lf\t%lf\n", &tmp1, &tmp2, &tmp3);
      rdata[i*3]=tmp1;
      rdata[i*3+1]=tmp2;
      rdata[i*3+2]=tmp3;
      subN++;
    }
    printf("Closing data file \n");
    fclose(data_file);

    a_data       = (double *)malloc(subN*sizeof(double));
    a_model      = (double *)malloc(subN*sizeof(double));
    t_data       = (double *)malloc(subN*sizeof(double));
    e_data       = (double *)malloc(subN*sizeof(double));
    q_data       = (double *)malloc(subN*sizeof(double));
    light_curve  = (double *)malloc(subN*sizeof(double));
    light_curve2 = (double *)malloc(subN*sizeof(double));
    mag_data     = (double *)malloc(5*sizeof(double));
    mag_err     = (double *)malloc(4*sizeof(double));
    
    for (subi=0;subi<Nt;subi++) {
      t_data[subi] = rdata[subi*3]; 
      a_data[subi] = rdata[subi*3+1]; 
      e_data[subi] = rdata[subi*3+2];
  }

    printf("Opening magnitude file");
    mag_file = fopen(mag_name, "r");
    double tmp1, tmp2;
    fscanf(mag_file, "%lf\n", &tmp1);
    mag_data[0] = tmp1;
    for (i=0;i<4;i++){
      fscanf(mag_file, "%lf\t%lf\n", &tmp1, &tmp2);
      mag_data[i+1] = tmp1;
      mag_err[i] = tmp2;
    }   
  fclose(mag_file); 
  }


	
  // initialize parallel tempering
  /* Siddhant: why is it set to an arbitrary value like this */
  dtemp = 1.2;
  if (burn_in >= 1) dtemp = 1.3;
  
  temp[0] = 1.0;
  index[0] = 0;
  
  for(i=1; i<NCHAINS; i++) {
    temp[i]  = temp[i-1]*dtemp;
    index[i] = i;
  }

  //initialize likelihood
  // read the likelihood from file
  //char *test_fname = "../debug/input_tests.txt";
  //FILE *outfile = fopen(test_fname, "r");
  //double tst_db;
  //for (int i=0; i<NPARS; i++)
  //{
    //fscanf(outfile, "\t%lf", &tst_db);
    //P_[0][i] = tst_db;
  //}

  //fclose(outfile);
  //TEST_lc_calc(P_[0]);

  logLmap = loglikelihood(t_data,a_data,e_data,subN,P_[0], mag_data, mag_err, 0);

  printf("initial chi2 and likelihood %g \t %g\n",-2*logLmap, logLmap);

  for(i=0; i<NCHAINS; i++) logLx[i] = logLmap;

  acc=0;
  printf("Creating chain and log files \n");
  chain_file = fopen(chainname,"w");
  logL_file  = fopen(logLname,"w");

  
  int atrial=0;
  int DEtrial=0;
  int DEacc=0;
  
  clock_t begin, end;
  /* Main MCMC Loop*/
  for (iter=0; iter<Niter; iter++) {
    if(iter % 10 == 0) {begin = clock();}
    //loop over chains
    #pragma omp parallel for schedule(static)
    for(j=0; j<NCHAINS; j++) {
      int jump;

      double *y = (double *)malloc(NPARS*sizeof(double));
      double alpha = ran2(&chain_seeds[j]);  
      double jscale = pow(10.,-6.+6.*alpha);
      /* propose new solution */
      jump=0;
      /* DE proposal; happens after 500 cycles */
      if(ran2(&chain_seeds[j])<0.5 && iter>NPAST) {jump=1;}
      //gaussian jumps along parameter directions
      if(jump==0){gaussian_proposal(x[index[j]], &chain_seeds[j], sigma, jscale, temp[j], y);}
      //jump along correlations derived from chain history
      if(jump==1) {
	      if(index[j]==0)DEtrial++;
	      differential_evolution_proposal(x[index[j]], &chain_seeds[j], history[j], y);
	      dx_mag=0;
	        for (i=0;i<NPARS;i++) {
	          dx_mag+=(x[index[j]][i]-y[i])*(x[index[j]][i]-y[i]);
	        }
	      if (dx_mag < 1e-6)
	        gaussian_proposal(x[index[j]], &chain_seeds[j], sigma, jscale, temp[j], y);
      }

      for (i=0;i<NPARS;i++) {

	      //enforce priors (reflecting boundary conditions)
	      if ((limited[i].lo == 1)&&(y[i] < limits[i].lo))
	        y[i] = 2.0*limits[i].lo - y[i];
	      if ((limited[i].hi == 1)&&(y[i] > limits[i].hi)) 
	        y[i] = 2.0*limits[i].hi - y[i];
	      //enforce priors (periodic boundary conditions)
	      if ((limited[i].lo == 2)&&(y[i] < limits[i].lo)) 
	        y[i] = limits[i].hi + (y[i]-limits[i].lo);
	      if ((limited[i].hi == 2)&&(y[i] > limits[i].hi)) 
	        y[i] = limits[i].lo + (y[i]-limits[i].hi);

        // Enforce limits on x too
	      //enforce priors (reflecting boundary conditions)
	      if ((limited[i].lo == 1)&&(x[index[j]][i] < limits[i].lo))
	        x[index[j]][i] = 2.0*limits[i].lo - x[index[j]][i];
	      if ((limited[i].hi == 1)&&(x[index[j]][i] > limits[i].hi)) 
	        x[index[j]][i] = 2.0*limits[i].hi - x[index[j]][i];
	      //enforce priors (periodic boundary conditions)
	      if ((limited[i].lo == 2)&&(x[index[j]][i] < limits[i].lo)) 
	        x[index[j]][i] = limits[i].hi + (x[index[j]][i]-limits[i].lo);
	      if ((limited[i].hi == 2)&&(x[index[j]][i] > limits[i].hi)) 
	        x[index[j]][i] = limits[i].lo + (x[index[j]][i]-limits[i].hi);
        
      }
      // Gaussian log priors
      double logPx, logPy;

      logPx = get_logP(x[index[j]], limited, limits, gauss_pars);
      logPy = get_logP(y, limited, limits, gauss_pars);

      //compute trial signal
      logLx[index[j]] = loglikelihood(t_data,a_data,e_data,subN,x[index[j]], mag_data, mag_err, 0);
      logLy = loglikelihood(t_data,a_data,e_data,subN,y, mag_data, mag_err, 0);

      /* evaluate new solution */
      //Hasting's ratio
      double H = exp( (logLy-logLx[index[j]])/temp[j] ) * pow(10., logPy - logPx);
      //acceptance probability
      alpha = ran2(&chain_seeds[j]);
      //printf("Log(Ly), Log(Lx): %12.5e %12.5e\n", logLy,logLx[index[j]]);
      //conditional acceptance of y
      if (alpha <= H) {
	      if(index[j]==0) acc++;
	      for (i=0;i<NPARS;i++) {
	        x[index[j]][i]=y[i];
	      }
        logLx[index[j]] = logLy;
	      if(jump==1 && index[j]==0) DEacc++;
      }

      /* parallel tempering */
      #pragma omp critical
      ptmcmc(index,temp,logLx);

      /*  fill history  */
      k = iter - (iter/NPAST)*NPAST;
      for(i=0; i<NPARS; i++) history[j][k][i] = x[index[j]][i];

      free_1d(y); 
      //#pragma omp barrier
    }
    /********Chain Loop ends**********/
    //update map parameters
    if (logLx[index[0]] > logLmap) {
      for (i=0;i<NPARS;i++) {
	      xmap[i]=x[index[0]][i];
      }
      logLmap = logLx[index[0]];
    }
    // index[0]: index of coldest chain
    //down-sample chains until I get a better proposal
    if(iter%10==0) {
      //print parameter chains
      fprintf(chain_file,"%ld %.12g ",iter/10,logLx[index[0]]);
      for(i=0; i<NPARS; i++) fprintf(chain_file,"%.12g ",x[index[0]][i]);
      fprintf(chain_file,"\n");
      
      //print log likelihood chains
      fprintf(logL_file,"%ld ",iter/10);
      for(i=0; i<NCHAINS; i++) fprintf(logL_file,"%.12g ",logLx[index[i]]);
      fprintf(logL_file,"\n"); 
    }
    
    //update progress to screen
    if(iter%100==0) {
	    printf("%ld/%ld logL=%.10g acc=%.3g DEacc=%.3g",iter,Niter,logLx[index[0]],
	    (double)(acc)/((double)atrial),
	    (double)DEacc/(double)DEtrial);
	    printf("\n");
      printf("Npars = %d \n", NPARS);
      }
    atrial++;
    
    if(iter%1000==0) acc=atrial=0;
    if(iter%100==0) {
      acc=atrial=0;
      data_file = fopen(outname,"w");
      calc_light_curve(t_data,subN,xmap,a_model);
      fprintf(data_file,"%ld\n",subN);
      for (i=0;i<subN;i++) {
	      fprintf(data_file,"%12.5e %12.5e %12.5e\n",t_data[i],a_data[i],a_model[i]);
      }
      fclose(data_file);  
      
      param_file = fopen(subparname,"w");
      for (int z=0; z<NPARS; z++){
        fprintf(param_file, "%12.5e ", x[index[0]][z]);
      }
      fprintf(param_file,"\n");
      fclose(param_file);
    }
    
  }
  /***** MCMC Loop ends *****/

  data_file = fopen(outname,"w");
  calc_light_curve(t_data,subN,xmap,a_model);
  fprintf(data_file,"%ld\n",subN);
  for (i=0;i<subN;i++) {
    fprintf(data_file,"%12.5e %12.5e %12.5e\n",t_data[i],a_data[i],a_model[i]);
  }
  fclose(data_file);  
  param_file = fopen(parname,"w");
  for (int z=0; z<NPARS; z++){
    fprintf(param_file, "%12.5e ", x[index[0]][z]);
  }
  fprintf(param_file,"\n");
  fclose(param_file);

  // Free the memory from arrays
  free_2d(P_, NCHAINS);
  free_2d(x, NCHAINS);
  free_3d(history, NCHAINS, NPAST);
  free_1d(sigma);
  free_1d(rdata);
  free(index);
  free_1d(a_data);
  free_1d(a_model);
  free_1d(t_data);
  free_1d(e_data);
  free_1d(q_data);
  free_1d(light_curve);
  free_1d(light_curve2);
  //free_1d(y);
  return(0);
}


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
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
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


void uniform_proposal(double *x, long *seed, bounds limits[], double *y)
{
	int n;
	double dx[NPARS];
	
	//compute size of jumps
	for(n=0; n<NPARS; n++) dx[n] = ran2(seed)*(limits[n].hi - limits[n].lo);
	
	//uniform draw on prior range
	for(n=0; n<NPARS; n++) y[n] = limits[n].lo + dx[n];
}

void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y)
{
	int n;
	double gamma;
	double sqtemp;
	double dx[NPARS];
	
	//scale jumps by temperature
	sqtemp = sqrt(temp);
	
	//compute size of jumps
	for(n=0; n<NPARS; n++) dx[n] = gasdev2(seed)*sigma[n]*sqtemp*scale;
	
	//jump in parameter directions scaled by dx
	for(n=0; n<NPARS; n++) {
	  y[n] = x[n] + dx[n];
    
	  //printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",
	  //	 scale,sqtemp,sigma[n],x[n],dx[n]);
	}
}

void hypercube_proposal(double *x, long *seed, double *y)
{
	
}
void differential_evolution_proposal(double *x, long *seed, double **history, double *y)
{
	int n;
	int a;
	int b;
	double dx[NPARS];
	
	//choose two samples from chain history
	a = ran2(seed)*NPAST;
	b = a;
	while(b==a) b = ran2(seed)*NPAST;
	
	//compute vector connecting two samples
	for(n=0; n<NPARS; n++) {
	  dx[n] = history[b][n] - history[a][n];
	}
	//Blocks?
	
	//90% of jumps use Gaussian distribution for jump size
	if(ran2(seed) < 0.9) for(n=0; n<NPARS; n++) dx[n] *= gasdev2(seed)*GAMMA;

	//jump along vector w/ appropriate scaling
	for(n=0; n<NPARS; n++) y[n] = x[n] + dx[n];
}

void initialize_proposals(double *sigma, double ***history)
{
	int n,i,j;
	double junk;
	double x[NPARS];
	FILE *hfile;
	
	
	/*********************/    
	/* Gaussian Proposal */
	/*********************/    
	
	//jump sizes are set by hand based on what works
	sigma[0]  = 1.0e-1;  //log M1 (MSUN)
	sigma[1]  = 1.0e-1;  //log M2 (MSUN)
	sigma[2]  = 1.0e-1;  //log P (days)
	sigma[3]  = 1.0e-2;  //e
	sigma[4]  = 1.0e-2;  //inc (rad)
	sigma[5]  = 1.0e-2;  //Omega (rad)
	sigma[6]  = 1.0e-2;  //omega0 (rad)
	sigma[7]  = 1.0e+0;  //T0 (day)
	sigma[8]  = 1.0e-2;  //log rr1 normalization
	sigma[9]  = 1.0e-2;  //log rr2 normalization
  sigma[10] = 1.0e-3;   // mu 1
  sigma[11] = 1.0e-3;   // tau 1
  sigma[12] = 1.0e-3;   // mu 2
  sigma[13] = 1.0e-3;   // tau 2
  sigma[14] = 1.0e-2;   // ref 1
  sigma[15] = 1.0e-2;   // ref 2
  sigma[16] = 1.0e-3;   // extra alph 1
  sigma[17] = 1.0e-3;   // extra alph 2
  sigma[18] = 1.0e-2;   // temp 1
  sigma[19] = 1.0e-2;   // temp 2
  sigma[20] = 1.0e-2;   // blending
  sigma[21] = 1.0e-3;   // flux tune
	//if (burn_in == 0) for(i=0; i<NPARS; i++) sigma[i] /= 100;
	
	//for(i=0; i<NPARS; i++) sigma[i] /= ((double)(NPARS));
		
	
	/**********************/    
	/* Hypercube Proposal */
	/**********************/    
	//use sigma to set bin size?
		
		
	/**************************/    
	/* Differential Evolution */
	/**************************/    
		/*
	//read initial "history" samples from saved chain	
	hfile = fopen("history.dat","r");
	
	for(i=0; i<NPAST; i++)
	{
		//read in parameter from history file
		fscanf(hfile,"%i %lg",&n,&junk);
		for(n=0; n<NPARS; n++)fscanf(hfile,"%lg",&x[n]);
		
		//copy parameters across all chains
		for(j=0; j<NCHAINS; j++) for(n=0; j<NPARS; n++) history[j][i][n] = x[n];
	}
	
	fclose(hfile);
		 */
}

