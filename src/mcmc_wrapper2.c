/*********************************************************/
//similar to mcmc_wrapper.c, but uses likelihood3.c,
//with independently varying R1,R2
//instead of using main sequence scaling relations
/*********************************************************/
#include "mcmc_wrapper2.h"

int main(int argc, char* argv[])
{
  // Run variables
  long NITER;
  long iter;
  
  double logLy;
  double logPy;
  double logLmap; 
  double  *sigma;
  double **x;
  double ***history;
  double xmap[NPARS];  
  double logLx[NCHAINS];
  double logPx[NCHAINS];

  bounds limited[NPARS], limits[NPARS];
  gauss_bounds gauss_pars[NPARS];
  RNG_Vars states[NCHAINS];
  
  int acc = 0;
  int atrial = 0;
  int DEtrial = 0;
  int DEacc = 0;

  int DEtrial_arr[NCHAINS];
  int acc_arr[NCHAINS];
  int DEacc_arr[NCHAINS];

  long seeds[NCHAINS];
 
  // Parallel tempering
  double dtemp;
  double temp[NCHAINS];
  int *index;

  // Arrays and files to store and read lightcurve data
  double *rdata, *a_data, *a_model, *e_data, *t_data;
  double *mag_data, *mag_err;
  FILE *param_file, *data_file, *chain_file, *logL_file, *mag_file, *logfile;
  char  pfname[150] = "";
  char  dfname[150] = "";
  char  parname[150] = "";
  char  subparname[150] = "";
  char  outname[150] = "";
  char  chainname[150] = "";
  char  logLname[150] = "";
  char  RUN_ID[50] = "";
  char  mag_name[150] = "";
  char logname[150] = "";
  
  // Temporary variables
  double tmp, tmp1, tmp2, tmp3, tmp4;              
  long Nt, subN;

  // Miscellaneous
  int run;
  int weight;
  double log_LC_PERIOD;
  char prefix[100] = "";
  char suffix[100] = "";
  char run_num[15] = "";

  NITER = (long)atoi(argv[1]);
  strcpy(RUN_ID,argv[2]);;
  log_LC_PERIOD = (double)atof(argv[3]);
  run = atoi(argv[4]);


  const double LC_PERIOD = pow(10., log_LC_PERIOD);

  if (ENABLE_OPENMP)
  {
    // 25 threads give optimum speed for a 1e4 likelihood evaluations
    printf("Using OPENMP (nthreads = 25)\n");
    omp_set_num_threads(25);
  }

  // Initialze seeds and other vairables
  srand(NITER);

  for (int i=0; i<NCHAINS; i++)
  {
    // To differentiate the starting seeds in openmp runs
    seeds[i] = i + run;
    DEtrial_arr[i] = 0;
    DEacc_arr[i] = 0;
    acc_arr[i] = 0;

    // for parallel random number generation
    states[i].idum2 = 123456789;
    states[i].iy =  0;
    for (int l=0; l<NTAB; l++)
    {
      states[i].iv[l] = 0;
    }
    states[i].iset =  0;
    states[i].gset = 0.;
    states[i].cts = 0;
  }


  // Initialize the file names
  sprintf(prefix, "/scratch/ssolanski/HB_MCMC/data");
  sprintf(run_num, "_%d", run);
  strcat(suffix,RUN_ID);
  
  if (USE_COLOR_INFO)
  {
    strcat(suffix, "_color");
  }

  if (ENABLE_OPENMP)
  {
    strcat(suffix, "_OMP");
  }

  strcat(suffix, run_num);

  strcat(subparname, prefix);
  strcat(parname, prefix);
  strcat(chainname, prefix);
  strcat(logLname, prefix);
  strcat(logname, prefix);
  strcat(outname, prefix);
  strcat(mag_name, prefix);
  strcat(dfname, prefix);

  strcat(subparname,"/subpars/subpar.");
  strcat(parname,   "/pars/par.");
  strcat(chainname, "/chains/chain.");
  strcat(logLname,  "/logL/logL.");
  strcat(logname,   "/log/log.");
  strcat(outname,   "/lightcurves/mcmc_lightcurves/");
  strcat(mag_name,  "/magnitudes/");
  strcat(dfname,    "/lightcurves/folded_lightcurves/");

  strcat(subparname, suffix);
  strcat(parname, suffix);
  strcat(chainname, suffix);
  strcat(logLname, suffix);
  strcat(logname, suffix);
  strcat(outname, suffix);

  strcat(subparname, ".dat");
  strcat(parname, ".dat");
  strcat(chainname, ".dat");
  strcat(logLname, ".dat");
  strcat(logname, ".dat");
  strcat(outname, ".out");

  strcat(mag_name,RUN_ID);
  strcat(mag_name,".txt");
  strcat(dfname,RUN_ID);
  strcat(dfname,"_new.txt");

  printf("Parfile: %s\n", parname);
  printf("Subparfile: %s\n", subparname);
  printf("Chainfile: %s\n", chainname);
  printf("logLfile: %s\n", logLname);
  printf("outfile: %s\n", outname);
  printf("logfile: %s\n", logname);

  // Allocate memory for the mcmc arrays
  x = (double **)malloc(NCHAINS*sizeof(double));
  for(int i=0;i<NCHAINS;i++) 
  {
    x[i]=(double *)malloc(NPARS*sizeof(double));
  }
  
  history = (double ***)malloc(NCHAINS*sizeof(double));
  for(int i=0;i<NCHAINS;i++)
  {
    history[i]=(double **)malloc(NPAST*sizeof(double));
    for(int j=0;j<NPAST;j++)
    {
      history[i][j]=(double *)malloc(NPARS*sizeof(double));
    }
  }

  sigma = (double *)malloc(NPARS*sizeof(double));

  // Initialize priors
  set_limits(limited, limits, gauss_pars, LC_PERIOD);
  
  // Set up proposal distribution
  initialize_proposals(sigma, history);

  // Read data from the parameter file
  printf("Opening parameter file \n");

  if ((exists(parname)) && (!USE_RAND_PARS))
  {
    printf("Reading contents from parameter file \n");
    param_file = fopen(parname, "r");
    for (int i=0;i<NPARS;i++) 
    {
      fscanf(param_file,"%lf\n", &tmp); 
      printf("Par %d value: %.5e \n",i, tmp);

      if (i == 2) 
      {
          tmp = log_LC_PERIOD;  // Keep the period from input
      }
      if (i == 7)
      {
          tmp = fmod(tmp, LC_PERIOD);
      }

      x[0][i]  = tmp;
      xmap[i]  = tmp;

      for(int j=1; j<NCHAINS; j++) 
      {
        x[j][i]  = tmp;
      }
    }

    fclose(param_file);
  }

  else 
  {
    printf("Parameter file not found, assigning random pars \n");
    for (int i=0;i<NPARS;i++) 
    {
      for(int j=0; j<NCHAINS; j++) 
      {
        tmp1 = ran2_parallel(&seeds[0], &states[0]);
        x[j][i] = (limits[i].lo + tmp1*(limits[i].hi - limits[i].lo));
        if (i == 2) 
        {
          x[j][i] = log_LC_PERIOD;  // Keep the period from input
        }
        if (i == 7)
        {
          x[j][i] = fmod(x[j][i], LC_PERIOD);
        }
      }
    }
  }

  // Open the folded lightcurve file
  printf("Opening folded lc data file %s \n", dfname);

  if (exists(dfname))
  {
    data_file = fopen(dfname,"r");
    fscanf(data_file,"%ld\n", &Nt);
    subN = 0;
    rdata        = (double *)malloc(3*Nt*sizeof(double));
    index        = (int *)malloc(NCHAINS*sizeof(int));

    for (int i=0; i<Nt; i++) 
    {
      fscanf(data_file,"%lf\t%lf\t%lf\n", &tmp1, &tmp2, &tmp3);
      rdata[i*3]=tmp1;
      rdata[i*3+1]=tmp2;
      rdata[i*3+2]=tmp3;
      subN++;
    }

    printf("Closing folded lc data file \n");
    fclose(data_file);

  }

  else
  {
    printf("Lightcurve datafile not found; terminating program \n");
    exit(0);
  }

  // Assign memory for the other arrays
  a_data       = (double *)malloc(subN*sizeof(double));
  a_model      = (double *)malloc(subN*sizeof(double));
  t_data       = (double *)malloc(subN*sizeof(double));
  e_data       = (double *)malloc(subN*sizeof(double));
  mag_data     = (double *)malloc(5*sizeof(double));
  mag_err     = (double *)malloc(4*sizeof(double));
  
  for (int subi= 0; subi<Nt; subi++) 
  {
    t_data[subi] = rdata[subi*3]; 
    a_data[subi] = rdata[subi*3+1]; 
    e_data[subi] = rdata[subi*3+2];
  }

  // Read the magnitude data
  printf("Opening magnitude file");
  if (exists(mag_name) && (USE_COLOR_INFO))
  {
    printf("Using color information \n");
    weight = 1;
    mag_file = fopen(mag_name, "r");
    fscanf(mag_file, "%lf\n", &tmp1);
    mag_data[0] = tmp1;
    for (int i=0;i<4;i++)
    {
      fscanf(mag_file, "%lf\t%lf\n", &tmp1, &tmp2);
      mag_data[i+1] = tmp1;
      mag_err[i] = tmp2;
    }   
  
  fclose(mag_file); 

  }

  else
  {
    printf("Magnitude file not found/used; assigning infinite error to mag data \n");
    weight = 0;
    mag_data[0] = 1000.;
    for (int i=0;i<4;i++)
    {
      mag_data[i+1] = 1.;
      mag_err[i] = BIG_NUM;
    } 
  }

  // Initialize parallel tempering
  dtemp = 1.4;
  temp[0] = 1.0;
  index[0] = 0;
  
  for(int i=1; i<NCHAINS; i++) 
  {
    temp[i]  = temp[i-1]*dtemp;
    index[i] = i;
  }

  // Perform first likelihood evaluation
  logLmap = loglikelihood(t_data,a_data,e_data,subN, x[0], mag_data, mag_err, weight);
  
  for(int i=0; i<NCHAINS; i++) 
  {
    logLx[i] = logLmap;
  }
  
  printf("initial chi2 and likelihood %lf \t %lf\n",-2*logLmap, logLmap);

  if (STORE_DATA)
  {
    // Make chain and Log files
    printf("Creating chain and log files %s and %s \n", chainname, logLname);
    chain_file = fopen(chainname,"w");
    logL_file  = fopen(logLname,"w");
    logfile    = fopen(logname, "w");
  }

  FILE *temp_log_files[NCHAINS];
  char temp_log_fname[NCHAINS][150];

  for (int j=0; j<NCHAINS; j++)
  {
      sprintf(temp_log_fname[j], "/scratch/ssolanski/HB_MCMC/debug/temp_%d_log.txt", j);
      temp_log_files[j] = fopen(temp_log_fname[j], "w");
      fclose(temp_log_files[j]);

      // Store the values of parameters in the chains
      temp_log_files[j] = fopen(temp_log_fname[j], "a");
  }

  // To monitor temperature swaps
  FILE *temp_swap_file = fopen("/scratch/ssolanski/HB_MCMC/debug/temp_swap_file.txt", "w");

  printf("Begining main mcmc loop \n");
  /* Main MCMC Loop*/
  for (iter=0; iter<NITER; iter++) 
  {

    int k = iter - (iter / NPAST) * NPAST;

    #pragma omp parallel for schedule(static) if(ENABLE_OPENMP) private(logLy, logPy)
    for(int j=0; j<NCHAINS; j++) 
    {
      // Test parameters
      double y[NPARS];

      // Random number
      double alpha = ran2_parallel(&seeds[j], &states[j]);

      // Jump scale and steps
      double jscale = pow(10.,-6.+6.*alpha);

      double dx[NPARS];
      double dx_mag = 0;
      int jump = 0;
      int chain_id = index[j];
      int jump_type = 0;


      // Take steps in parameter space
      if((ran2_parallel(&seeds[j], &states[j]) < 0.5) && (iter > NPAST)) 
      {
        jump = 1;
      }

      //gaussian jumps
      if(jump == 0)
      {
        gaussian_proposal_parallel(x[chain_id], &seeds[j], sigma, jscale, temp[j], y, &states[j]);
        jump_type = 1;
      }

      //jump along correlations derived from chain history
      if(jump == 1)
      {
        /* DE proposal; happens after 500 cycles */
        if (chain_id == 0)
        {
          DEtrial_arr[j]++;
        }

        differential_evolution_proposal_parallel(x[chain_id], &seeds[j], history[j], y, &states[j]);
        jump_type = 2;

        for (int i=0;i<NPARS;i++) 
        {
          dx_mag += (x[chain_id][i] - y[i]) * (x[chain_id][i] - y[i]);
        }

        if (dx_mag < 1e-6)
        {
          gaussian_proposal_parallel(x[chain_id], &seeds[j], sigma, jscale, temp[j], y, &states[j]);
          jump_type = 1;
        }
      }

      // Enforce priors
      for (int i=0;i<NPARS;i++) 
      {
	      // Reflecting boundary conditions
	      while (((limited[i].lo == 1) && (y[i] < limits[i].lo)) || ((limited[i].hi == 1) && (y[i] > limits[i].hi)))
        {
          if (y[i] < limits[i].lo)
          {
            y[i] = 2.0*limits[i].lo - y[i];
          }

          else
          {
            y[i] = 2.0*limits[i].hi - y[i]; 
          }
        }

	      // Periodic boundary conditions
	      while ((limited[i].lo == 2) && (y[i] < limits[i].lo))
        {
	        y[i] = limits[i].hi + (y[i]-limits[i].lo);
        }

	      while ((limited[i].hi == 2) && (y[i] > limits[i].hi))
        {
	        y[i] = limits[i].lo + (y[i]-limits[i].hi);
        }

      }

      // Order the masses
      if (y[1] > y[0])
      {
        double tmp = y[1];
        y[1] = y[0];
        y[0] = y[1];
      }

      // Fix the period
      y[2] = log_LC_PERIOD;

      // Make the phase modulo period
      y[7] = fmod(y[7], LC_PERIOD);

      // Gaussian priors
      logPx[chain_id] = get_logP(x[chain_id], limited, limits, gauss_pars);
      logPy = get_logP(y, limited, limits, gauss_pars);

      //compute current and trial likelihood
      logLx[chain_id] = loglikelihood(t_data, a_data, e_data, subN, x[chain_id], mag_data, mag_err, weight);
      logLy           = loglikelihood(t_data, a_data, e_data, subN, y, mag_data, mag_err, weight);

      /* evaluate new solution */
      alpha = ran2_parallel(&seeds[j], &states[j]);

      //Hasting's ratio
      double H = exp((logLy-logLx[chain_id])/temp[j] + (logPy-logPx[chain_id]));

      //
      //for (int i=0; i<NPARS; i++)
      //{
        // More file printing 
      //  fprintf(temp_log_files[j], "%lf\t%lf\t", x[chain_id][i], y[i]);
      //}

      //conditional acceptance of y
      if (alpha <= H) 
      {

        // print if the step was accepted
        //fprintf(temp_log_files[j], "%ld\n", jump_type);

        // Warn if the best likelihood evaluation looks weird for the 5 coldest chains
        if ( (logLx[chain_id] / logLy <= 0.5)  && (iter > 10000) && (j <= 5))
        {
          LogSuspiciousJumps(logfile, iter, chain_id, H, alpha, temp[j], logLx[chain_id], logLy,
                             logPx[chain_id], logPy, x[chain_id], y, jump_type);
        }

	      if (chain_id == 0) 
        {
          acc_arr[j]++;
        }

	      for (int i=0;i<NPARS;i++) 
        {
	        x[chain_id][i] = y[i];
	      }

        logLx[chain_id] = logLy;

	      if((jump == 1) && (chain_id == 0)) 
        {
          DEacc_arr[j]++;
        }

      }

    //else
    //  {
        // Print if the step was accepted or not
        //fprintf(temp_log_files[j], "%ld\n", 0);
    //  }

    for(int i=0; i<NPARS; i++) 
    {
      history[j][k][i] = x[chain_id][i];
    }

    //#pragma omp critical
    //fprintf(logfile, "Iter: %d Thread id: %d Loop iter: %d index[j]: %d seed: %d cts: %d \n", iter, omp_get_thread_num(), j, index[j], seeds[j], states[j].cts);
  }
  
    /********Chain Loop ends**********/

    for (int i=0; i<NCHAINS; i++)
    {
      acc += acc_arr[i];
      DEacc += DEacc_arr[i];
      DEtrial += DEtrial_arr[i];
      acc_arr[i] = 0;

      /* parallel tempering */
      ptmcmc(index, temp, logLx, logPx, temp_swap_file);
    }
    //update map parameters
    if (logLx[index[0]] > logLmap)
    {
      for (int i=0;i<NPARS;i++) 
      {
	      xmap[i]=x[index[0]][i];
      }
      logLmap = logLx[index[0]];
    }
    
    //update progress to screen
    if(iter%1000==0) 
    {
	    printf("%ld/%ld logL=%.10g acc=%.3g DEacc=%.3g",iter,NITER,logLx[index[0]],
	    (double)(acc)/((double)atrial),
	    (double)DEacc/(double)DEtrial);
	    printf("\n");

      printf("Parameter values: \n");
      // Print first few parameters
      for (int i=0; i<5; i++)
      {
        printf("%lf\t", x[index[10]][i]);
      }
      printf("\n");
    }

    atrial++;

    if((iter%1000==0) && (STORE_DATA)) 
    {
      //print parameter chains
      fprintf(chain_file,"%ld %.12g ",iter/10,logLx[index[0]]);

      for(int i=0; i<NPARS; i++) 
      {
        fprintf(chain_file,"%.12g ",x[index[0]][i]);
      }

      fprintf(chain_file,"\n");
      
      //print log likelihood chains
      fprintf(logL_file,"%ld ",iter/10);

      for(int i=0; i<NCHAINS; i++)
      {
        fprintf(logL_file,"%.12g ",logLx[index[i]]);
        
        for(int j=0; j<NPARS; j++)
        {
          fprintf(temp_log_files[i], "%lf\t", x[index[i]][j]);
        }
        
        fprintf(temp_log_files[i], "\n");
      }

      fprintf(logL_file,"\n");

      acc=atrial=0;

      for (int i=0; i<NCHAINS; i++)
      {
        DEtrial_arr[i] = 0;
        DEacc_arr[i] = 0;
        acc_arr[i] = 0;
      }

      data_file = fopen(outname,"w");
      calc_light_curve(t_data,subN,xmap,a_model);
      fprintf(data_file,"%ld\n",subN);

      for (int i=0;i<subN;i++) 
      {
	      fprintf(data_file,"%12.5e %12.5e %12.5e\n",t_data[i],a_data[i],a_model[i]);
      }

      fclose(data_file);  
  
      param_file = fopen(subparname,"w");
      for (int z=0; z<NPARS; z++)
      {
        fprintf(param_file, "%12.5e ", x[index[0]][z]);
      }
      fprintf(param_file,"\n");
      fclose(param_file);
    }
  }

  fclose(temp_swap_file);
  /***** MCMC Loop ends *****/

  if (STORE_DATA)
  {
    // Print computed lightcurve
    data_file = fopen(outname,"w");
    calc_light_curve(t_data,subN,xmap,a_model);
    fprintf(data_file,"%ld\n",subN);

    for (int i=0;i<subN;i++) 
    {
      fprintf(data_file,"%12.5e %12.5e %12.5e\n",t_data[i],a_data[i],a_model[i]);
    }

    // Print best parameters
    param_file = fopen(parname,"w");
    for (int z=0; z<NPARS; z++)
    {
      fprintf(param_file, "%12.5e ", x[index[0]][z]);
    }

    fprintf(param_file,"\n");

    fclose(logfile);
    fclose(data_file);  
    fclose(param_file);
    fclose(chain_file);
    fclose(logL_file);
  }

  for (int j=0; j<NCHAINS; j++)
  {
    fclose(temp_log_files[j]);
  }

  // Free the memory from arrays
  free_2d(x, NCHAINS);
  free_3d(history, NCHAINS, NPAST);
  free_1d(sigma);
  free_1d(rdata);
  free(index);
  free_1d(a_data);
  free_1d(a_model);
  free_1d(t_data);
  free_1d(e_data);
  return(0);
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
      sigma = (limits[i].hi - limits[i].lo) / 3.;
      logP += log(gaussian(pars[i], mean, sigma));
    }
  }
  return logP;
}


void ptmcmc(int *index, double temp[], double logL[], double logP[], FILE *temp_swap_file)
{
  int a, b;
	int olda, oldb;
	
	double heat1, heat2;
	double logL1, logL2;
  double logP1, logP2;
	double dlogL;
  double dlogP;
  double dlogE;
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
	b = (int) (((double)rand()/(RAND_MAX))*((double)(NCHAINS-1)));
	a = b + 1;
	
	olda = index[a];
	oldb = index[b];
	heat1 = temp[a];
	heat2 = temp[b];
	logL1 = logL[olda];
	logL2 = logL[oldb];
  logP1 = logP[olda];
  logP2 = logP[oldb];
	dlogL = logL2 - logL1;
  dlogP = logP1 - logP2;
	H  = (heat2 - heat1)/(heat2*heat1);
	alpha = exp(dlogL*H);
	beta  = ((double)rand()/(RAND_MAX));
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


double ran2_parallel(long *idum, RNG_Vars* state)
{
  if (ENABLE_OPENMP)
  {
    int j;
    long k;
    double temp;

    // Update counter
    state-> cts += 1;
    
    if (*idum <= 0) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      state->idum2=(*idum);
      for (j=NTAB+7;j>=0;j--) {
        k=(*idum)/IQ1;
        *idum=IA1*(*idum-k*IQ1)-k*IR1;
        if (*idum < 0) *idum += IM1;
        if (j < NTAB) state->iv[j] = *idum;
      }
      state->iy=state->iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=state->idum2/IQ2;
    state->idum2=IA2*(state->idum2-k*IQ2)-k*IR2;
    if (state->idum2 < 0) state->idum2 += IM2;
    j=state->iy/NDIV;
    state->iy=state->iv[j]-state->idum2;
    state->iv[j] = *idum;
    if (state->iy < 1) state->iy += IMM1;

    if ((temp=AM*state->iy) > RNMX) 
    {
      return RNMX;
    }
    else 
    {
      return temp;
    }
  }

  else
  {
    double temp = ran2(idum);
    return temp;
  }
}


// gaussian random number
double gasdev2_parallel(long *idum, RNG_Vars *state)
{
  if (ENABLE_OPENMP)
  {
    double fac, rsq, v1, v2;
    
    if(*idum < 0) state->iset = 0;
    if(state->iset == 0){
      do{
        v1 = 2.0 * ran2_parallel(idum, state)-1.0;
        v2 = 2.0 * ran2_parallel(idum, state)-1.0;
        rsq = v1*v1+v2*v2;
      } while (rsq >= 1.0 || rsq == 0.0);
      fac=sqrt(-2.0*log(rsq)/rsq);
      state->gset = v1 * fac;
      state->iset = 1;
      return(v2*fac);
    } else {
      state->iset = 0;
      return (state->gset);
    }
  }

  else
  {
    gasdev2(idum);
  }
}


void uniform_proposal(double *x, long *seed, bounds limits[], double *y)
{
	int n;
	double dx[NPARS];
	
	//compute size of jumps
	for(n=0; n<NPARS; n++) 
  {
    dx[n] = ran2(seed)*(limits[n].hi - limits[n].lo);
  }
	
	//uniform draw on prior range
	for(n=0; n<NPARS; n++) 
  {
    y[n] = limits[n].lo + dx[n];
  }
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
	for(n=0; n<NPARS; n++)
  {
	  y[n] = x[n] + dx[n];
	}
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


void uniform_proposal_parallel(double *x, long *seed, bounds limits[], double *y, RNG_Vars *state)
{
  if (ENABLE_OPENMP)
  {
    int n;
    double dx[NPARS];
    
    //compute size of jumps
    for(n=0; n<NPARS; n++) dx[n] = ran2_parallel(seed, state)*(limits[n].hi - limits[n].lo);
    
    //uniform draw on prior range
    for(n=0; n<NPARS; n++) y[n] = limits[n].lo + dx[n];
  }

  else
  {
    uniform_proposal(x, seed, limits, y);
  }
}

void gaussian_proposal_parallel(double *x, long *seed, double *sigma, double scale, double temp, double *y, RNG_Vars *state)
{
  if (ENABLE_OPENMP)
  {
    int n;
    double gamma;
    double sqtemp;
    double dx[NPARS];
    
    //scale jumps by temperature
    sqtemp = sqrt(temp);
    
    //compute size of jumps
    for(n=0; n<NPARS; n++) dx[n] = gasdev2_parallel(seed, state)*sigma[n]*sqtemp*scale;
    
    //jump in parameter directions scaled by dx
    for(n=0; n<NPARS; n++) 
    {
      y[n] = x[n] + dx[n];
    }
  }

  else
  {
    gaussian_proposal(x, seed, sigma, scale, temp, y);
  }
}


void differential_evolution_proposal_parallel(double *x, long *seed, double **history, double *y, RNG_Vars *state)
{
  if (ENABLE_OPENMP)
  {
    int n;
    int a;
    int b;
    int c;
    double dx[NPARS];
    double epsilon[NPARS];
    
    //choose two samples from chain history
    a = ran2_parallel(seed, state)*NPAST;
    a = ran2_parallel(seed, state);
    b = a;
    while(b==a) 
    {
      b = ran2_parallel(seed, state)*NPAST;
    }

    //compute vector connecting two samples
    for(n=0; n<NPARS; n++) 
    {
      dx[n] = history[b][n] - history[a][n];
      epsilon[n] = dx[n] * (gaussian(c, 0, 1.e-4) - 0.5);
    }
    //Blocks?
    
    //90% of jumps use Gaussian distribution for jump size
    if(ran2_parallel(seed, state) < 0.9) 
    {
      for(n=0; n<NPARS; n++) 
      {
        dx[n] *= gasdev2_parallel(seed, state) * GAMMA;
      }
    }

    //jump along vector w/ appropriate scaling
    for(n=0; n<NPARS; n++) 
    {
      dx[n] += epsilon[n];
      y[n] = x[n] + dx[n];
    }
  }

  else
  {
    differential_evolution_proposal(x, seed, history, y);
  }
}

/* Set priors on parameters, and whether or not each parameter is bounded*/
/* Siddhant: Maybe just use an if condition/switch statment instead of limited*/
void set_limits(bounds limited[], bounds limits[], gauss_bounds gauss_pars[], double LC_PERIOD)
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
  limited[4].lo = 1;
  limits[4].lo = 0;
  limited[4].hi = 1;
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
  limits[7].lo = 0.;
  limited[7].hi = 1;
  limits[7].hi = LC_PERIOD;
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
      // limits on (log) Teff coefficient for star 1
      limited[18].lo = 1;
      limits[18].lo = -3.;
      limited[18].hi = 1;
      limits[18].hi = 3.;
      gauss_pars[18].flag = 1.;
      // limits on (log) Teff coefficient for star 2
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

void initialize_proposals(double *sigma, double ***history)
{
	int n,i,j;
	double junk;
	double x[NPARS];
	
	
	/*********************/    
	/* Gaussian Proposal */
	/*********************/    
	
	//jump sizes are set by hand based on what works
	sigma[0]  = 1.0e-2;  //log M1 (MSUN)
	sigma[1]  = 1.0e-2;  //log M2 (MSUN)
	sigma[2]  = 1.0e-8;  //log P (days)
	sigma[3]  = 1.0e-2;  //e
	sigma[4]  = 1.0e-3;  //inc (rad)
	sigma[5]  = 1.0e-3;  //Omega (rad)
	sigma[6]  = 1.0e-3;  //omega0 (rad)
	sigma[7]  = 1.0e-5;  //T0 (day)
	sigma[8]  = 1.0e-1;  //log rr1 normalization
	sigma[9]  = 1.0e-1;  //log rr2 normalization
  sigma[10] = 1.0e-2;   // mu 1
  sigma[11] = 1.0e-2;   // tau 1
  sigma[12] = 1.0e-2;   // mu 2
  sigma[13] = 1.0e-2;   // tau 2
  sigma[14] = 1.0e-2;   // ref 1
  sigma[15] = 1.0e-2;   // ref 2
  sigma[16] = 1.0e-2;   // extra alph 1
  sigma[17] = 1.0e-2;   // extra alph 2
  sigma[18] = 1.0e-1;   // temp 1
  sigma[19] = 1.0e-1;   // temp 2
  sigma[20] = 1.0e-3;   // blending
  sigma[21] = 1.0e-5;   // flux tune
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

/* Other functions*/
/* Free memory from arrays */
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

// Standard gaussian
double gaussian(double x, double mean, double sigma)
{
  return (1 / sigma / SQRT_2PI) * exp(- pow((x - mean) / sigma, 2.));
}

/*
Function to test the lightcurve calculation
*/
void TEST_lc_calc(double *pars)
{
	// create a synthetic times array
	double times[100];
	for (int i=0; i<100; i++)
	{
		times[i] = (double) (i);
	}

	double lc_template[100];
	calc_light_curve(times, 100, pars, lc_template);

	// Write to a file
	FILE* outfile;
	char *outname = "../debug/lightcurve_calc_original.txt";
	outfile = fopen(outname, "w");

  // parameters
  for(int i=0; i<NPARS; i++)
  {
    fprintf(outfile, "%.10g\n", pars[i]);
  }

  // lightcurve
	for (int i=0; i<100; i++)
	{
		fprintf(outfile, "%.10g\n", lc_template[i]);
	}

  // magnitude calculation
  double D = 100.;
  double m1, m2, m3, m4;
  calc_mags(pars, D, &m1, &m2, &m3, &m4);
  fprintf(outfile, "Blending %.10g\n", pars[20]);
  fprintf(outfile, "%.10g\n", m1);
  fprintf(outfile, "%.10g\t", m2);
  fprintf(outfile, "%.10g\t", m3);
  fprintf(outfile, "%.10g\t", m4);

	fclose(outfile);
	
	return;
}


// log chain data on the suspicious big likelihood jumps
// chain id corresponds to index[j] and temp corresponds to temp[j], or 1.3**(j+1)
void LogSuspiciousJumps(FILE* logfile, long iter, int chain_id, double H, double alpha, double chain_temp, double logL_old, double logL_new,
                        double logP_old, double logP_new, double *pars_old, double *pars_new, int jump_type)
{
  fprintf(logfile, "Big jump in likelihood detected on iternation: %ld and chain id: %ld ;", iter, chain_id);
  fprintf(logfile, " temperature of the chain: %f \n", chain_temp);
  fprintf(logfile, "Old log prior: %f new log prior: %f old log likelihood: %f new log likelihood: %f \n", logP_old, logP_new, 
          logL_old, logL_new);
  fprintf(logfile, "Hastings ratio [exp((logLy-logLx[chain_id])/temp[j]) * pow(10., logPy - logPx)] %f, alpha %f \
          and jump type %d\n", H, alpha, jump_type);
  fprintf(logfile, "Printing old and new parameters \n");
  
  for (int i=0; i<NPARS; i++)
  {
    fprintf(logfile, "%lf \t", pars_old[i]);
  }
  
  fprintf(logfile, "\n");
  
  for (int i=0; i<NPARS; i++)
  {
    fprintf(logfile, "%lf \t", pars_new[i]);
  }
  
  fprintf(logfile, "\n ************************************************ \n");
  return;
}