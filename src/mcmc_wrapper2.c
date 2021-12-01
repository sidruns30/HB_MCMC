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
#include <omp.h>
#include<unistd.h>
#include "likelihood3.c"

#define NCHAINS 50
#define NPAST   500
#define GAMMA (2.388/sqrt(2.*NPARS))
#define NUM_ELEMENTS(x) (size_of(x)/size_of(x[0]))


/*Change SIGMAP if necessary*/
double SIGMAP;

/* From likelihood.c */
void calc_light_curve(double t_data[], long Nt, double P_[], double light_curve[]);
double likelihood(double t_data[], double a_data[], double e_data[], long Nt, double P_[]);
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
    if ((file = fopen(fname, "r"))){
        fclose(file);
        return 1;
    }
    return 0;
}

int main(int
 argc, char* argv[])
{
  //omp_set_num_threads(1);
  /* Siddhant: Adding small variable descriptions */
  long Niter;                     // Chain iterations
  double **P_;                    // Contains parameters - check why different from x - seems uneccessary
  double P_0[NPARS];            
  double *y;                      // Updated Parameter
  double **x;                     // Parameter chains
  double ***history;              // Chain history
  double xmap[NPARS];           // Contains the chain with best paramters
  double dx[NPARS];             // Parameter step
  double dx_mag;                  
  double alpha;                   // Random number stuff
  double logLx[NCHAINS];          // Log likelihood for all chains
  double logLmap;                 
  double logLy;                   // Likelihood for new step
  double tmp,tmp1,tmp2,tmp3,tmp4;
  //double Tmin;                    // Never used!
  //double t_range[2];              // Never used!
  double true_err,Nerr,suberr;    // True error contains error on each LC point; others not used
  double scale;                   
  double dtemp;                   // Ask Jeremy
  double heat[NCHAINS],temp[NCHAINS]; // More PTMCMC stuff - heat never used
  double H;                       // Hasting's ratio
  double *rdata;                  
  double *a_data,*a_model,*e_data,*t_data,*q_data;
  double *light_curve;            
  double *light_curve2;           
  long i,j,k,m,subi;              
  long Nt,Nstart,Nstop,subN,obs_id,nch;
  long acc;                       
  long iter,TICid;                
  int *index, burn_in;
  double ls_period;  
  double SIGMAP;        

  bounds limited[NPARS], limits[NPARS];
  FILE *param_file, *data_file, *chain_file, *logL_file;
  /* Siddhant: Initializing empty arrays help avoid memory leaks*/
  char  pfname[80]="", dfname[80]="",  parname[80]="",
        subparname[80]="",  outname[80]="", chainname[80]="",
        logLname[80]="", RUN_ID[11]="";
  //characteristic magnitude of MCMC step in each parameter
  double  *sigma;

  Niter = (long)atoi(argv[1]);
  strcpy(RUN_ID,argv[2]);
  true_err = (double)atof(argv[3]);
  burn_in = atoi(argv[4]);
  ls_period = (double)atof(argv[5]);

  // Set period search
  if (burn_in == 2) {SIGMAP = 1.e-10;}
  else              {SIGMAP = 1.e-1;}

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
  printf("%s\n",parname);

  /*Use binned period if period is found, otherwise use the full lightcurve*/
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
  
  history = (double ***)malloc(NCHAINS*sizeof(double));
  for(i=0;i<NCHAINS;i++) {
    history[i]=(double **)malloc(NPAST*sizeof(double));
    for(j=0;j<NPAST;j++)history[i][j]=(double *)malloc(NPARS*sizeof(double));
  }

  y = (double *)malloc(NPARS*sizeof(double));
  sigma = (double *)malloc(NPARS*sizeof(double));
  scale  = 1.0/((double)(NPARS));
  //

  //true_err = 1.0;
  srand(Niter);
  long seed = Niter;
  //read in starting parameters near best solution
  printf("Opening parameter file \n");
  if (exists(parname)){param_file = fopen(parname,"r");}
  printf("Parameter file: %s\n",parname);
  //initialize priors
  set_limits(limited,limits);
  //set up proposal distribution
  initialize_proposals(sigma, history);
  
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
        if (param_file_flag != 1){
          P_[j][i] = limits[i].lo + ran2(&seed)*(limits[i].hi - limits[i].lo);
          if (i == 2) P_[j][i] = ls_period;  // Keep the period
          if (i == 7) P_[j][i] = 0;          // Set T0 = 0
          x[j][i]  = P_[j][i];
        }
      }
    }
  }

  if (burn_in == 1){
      printf("Test6 \n");
      printf("Dfname:");
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

    printf("Test7 \n");
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
    
    for (subi=0;subi<Nt;subi++) {
      t_data[subi] = rdata[subi*3]; 
      a_data[subi] = rdata[subi*3+1]; 
      e_data[subi] = rdata[subi*3+2];
  }
  }


	
  // initialize parallel tempering
  /* Siddhant: what is dtemp and why is it set to an arbitrary value like this */
  dtemp = 1.2;
  if (burn_in >= 1) dtemp = 1.3;
  
  temp[0] = 1.0;
  index[0] = 0;
  
  for(i=1; i<NCHAINS; i++) {
    temp[i]  = temp[i-1]*dtemp;
    index[i] = i;
  }

  //initialize likelihood
  logLmap = loglikelihood(t_data,a_data,e_data,subN,P_[0]);

  printf("initial chi2 %g\n",-2*logLmap);
 
  //for(i=0; i<NCHAINS; i++) logLx[i] = logLmap;

  acc=0;
  printf("Creating chain and log files \n");
  chain_file = fopen(chainname,"w");
  logL_file  = fopen(logLname,"w");

  int atrial=0;
  int DEtrial=0;
  int DEacc=0;
  int jump;
  double jscale;
  
  clock_t begin, end;
  /* Main MCMC Loop*/
  for (iter=0; iter<Niter; iter++) {
    if(iter % 10 == 0) {begin = clock();}
    //loop over chains
    //#pragma omp for schedule(static)
    for(j=0; j<NCHAINS; j++) {
      alpha = ran2(&seed);  
      jscale = pow(10.,-6.+6.*alpha);
      /* propose new solution */
      jump=0;
      /* DE proposal; happens after 500 cycles */
      if(ran2(&seed)<0.5 && iter>NPAST) {jump=1;}
      //gaussian jumps along parameter directions
      if(jump==0){gaussian_proposal(x[index[j]], &seed, sigma, jscale, temp[j], y);}
      //jump along correlations derived from chain history
      if(jump==1) {
	      if(index[j]==0)DEtrial++;
	      differential_evolution_proposal(x[index[j]], &seed, history[j], y);
	      dx_mag=0;
	        for (i=0;i<NPARS;i++) {
	          dx_mag+=(x[index[j]][i]-y[i])*(x[index[j]][i]-y[i]);
	        }
	      if (dx_mag < 1e-6)
	        gaussian_proposal(x[index[j]], &seed, sigma, jscale, temp[j], y);
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
      }

      //compute trial signal
      logLx[index[j]] = loglikelihood(t_data,a_data,e_data,subN,x[index[j]]);
      logLy = loglikelihood(t_data,a_data,e_data,subN,y);
      
      /* evaluate new solution */
      //Hasting's ratio
      H = exp( (logLy-logLx[index[j]])/temp[j] );
      //acceptance probability
      alpha = ran2(&seed);
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
      ptmcmc(index,temp,logLx);

      /*  fill history  */
      k = iter - (iter/NPAST)*NPAST;
      for(i=0; i<NPARS; i++) history[j][k][i] = x[index[j]][i];
    }
    /********Chain Loop ends**********/
    // Call timer after 10 iterations
    //if (iter%10 == 9){
    //  end = clock();
    //  double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    //  printf("time spent in the last 10 iterations: %f \n", time_spent);
    //  }

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
      fprintf(param_file,"%12.5e %12.5e %12.5e %12.5e %12.5e\n",
	    x[index[0]][0],x[index[0]][1],x[index[0]][2],x[index[0]][3],x[index[0]][4]);
      fprintf(param_file,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	    x[index[0]][5],x[index[0]][6],x[index[0]][7],x[index[0]][8],x[index[0]][9],x[index[0]][10]);
      if (ALPHA_FREE == 1){
         fprintf(param_file,"%12.5e %12.5e %12.5e %12.5e %12.5e\n",
	       x[index[0]][11],x[index[0]][12],x[index[0]][13],x[index[0]][14],x[index[0]][15]);
      }
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
  fprintf(param_file,"%12.5e %12.5e %12.5e %12.5e %12.5e\n",
	x[index[0]][0],x[index[0]][1],x[index[0]][2],x[index[0]][3],x[index[0]][4]);
  fprintf(param_file,"%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
	x[index[0]][5],x[index[0]][6],x[index[0]][7],x[index[0]][8],x[index[0]][9],x[index[0]][10]);
  if (ALPHA_FREE == 1){
    fprintf(param_file,"%12.5e %12.5e %12.5e %12.5e %12.5e\n",
	  x[index[0]][11],x[index[0]][12],x[index[0]][13],x[index[0]][14],x[index[0]][15]);
  }
  fclose(param_file);

  // Free the memory from arrays
  free_2d(P_, NCHAINS);
  free_2d(x, NCHAINS);
  free_3d(history, NCHAINS, NPAST);
  free_1d(y);
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
	sigma[2]  = SIGMAP;  //log P (days)
	sigma[3]  = 1.0e-1;  //e
	sigma[4]  = 1.0e+0;  //inc (deg)
	sigma[5]  = 1.0e+1;  //Omega (deg)
	sigma[6]  = 1.0e+1;  //omega0 (deg)
	sigma[7]  = 1.0e+0;  //T0 (day)
	sigma[8]  = 1.0e-1;  //log flux normalization
	sigma[9]  = 1.0e-2;  //log rr1 normalization
	sigma[10]  = 1.0e-2;  //log rr2 normalization
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
