/*********************************************************/
//similar to mcmc_wrapper.c, but uses likelihood2.c,
//with independently varying R1,R2
//instead of using main sequence scaling relations
/*********************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "likelihood2.c"
//#define Niter 100000
//#define burn_in 1
#define Nparams 11
#define Nchains 50
#define Npast   500
#define sigmaP 1.0e-12
#define GAMMA (2.388/sqrt(2.*Nparams))

typedef struct {
  double lo;
  double hi;
} bounds;

void calc_light_curve(double t_data[], long Nt, double P_[], double light_curve[]);
double likelihood(double t_data[], double a_data[], double e_data[], long Nt, double P_[]);

void ptmcmc(int *index, double heat[], double logL[]);
double ran2(long *idum);
double gasdev2(long *idum);

void initialize_proposals(double *sigma, double ***history);
void uniform_proposal(double *x, long *seed, bounds limits[], double *y);
void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y);
void hypercube_proposal(double *x, long *seed, double *y);
void differential_evolution_proposal(double *x, long *seed, double **history, double *y);

/*Set priors on parameters, and whether or not each parameter is bounded*/
void set_limits(bounds limited[], bounds limits[])
{
  //limits on M1, in log10 Msun
  limited[0].lo = 1; 
  limits[0].lo = -1.5;
  limited[0].hi = 1;
  limits[0].hi = 2.0;
  
  //limits on M2, in log10 Msun
  limited[1].lo = 1;
  limits[1].lo = -1.5;
  limited[1].hi = 1;
  limits[1].hi = 2.0;
	
  //limits on P, in log10 days
  limited[2].lo = 1;
  limits[2].lo = -2.0;
  limited[2].hi = 1;
  limits[2].hi = 3.0;
  
  //limits on e
  limited[3].lo = 1;
  limits[3].lo = 0.0;
  limited[3].hi = 1;
  limits[3].hi = 1.;
  
  //limits on inc, in degrees
  limited[4].lo = 2; //reflecting BCs
  limits[4].lo = -90.;
  limited[4].hi = 2;
  limits[4].hi = 90.;
  
  //limits on Omega, in degrees
  limited[5].lo = 2;
  limits[5].lo = -180.;
  limited[5].hi = 2;
  limits[5].hi = 180.;
  
  //limits on omega0, in degrees
  limited[6].lo = 2;
  limits[6].lo = -180;
  limited[6].hi = 2;
  limits[6].hi = 180.;
  
  //limits on T0, in MDJ-2450000
  limited[7].lo = 1;
  limits[7].lo = -1000;
  limited[7].hi = 1;
  limits[7].hi = 1000.;
  
  //limits on log Flux0, in TESS counts
  limited[8].lo = 1;
  limits[8].lo = -5.0;
  limited[8].hi = 1.;
  limits[8].hi = 5.0;

  //limits on log rr1, the scale factor for R1
  limited[9].lo = 1;
  limits[9].lo = -2.0;
  limited[9].hi = 1.;
  limits[9].hi = 2.0;

  //limits on log rr2, the scale factor for R2
  limited[10].lo = 1;
  limits[10].lo = -2.0;
  limited[10].hi = 1.;
  limits[10].hi = 2.0;
}

int main(int argc, char* argv[])
{
  long Niter;
  double **P_;
  double P_0[Nparams];
  double *y;
  double **x;
  double ***history;
  double xmap[Nparams];
  double dx[Nparams];
  double dx_mag;
  double alpha;
  double logLx[Nchains];
  double logLmap;
  double logLy;
  double tmp,tmp1,tmp2,tmp3,tmp4;
  double Tmin;
  double t_range[2];
  double true_err,Nerr,suberr;
  double scale;
  double dtemp;
  double heat[Nchains],temp[Nchains];
  double H;
  double *rdata;
  double *a_data,*a_model,*e_data,*t_data,*q_data;
  double *light_curve;
  double *light_curve2;
  long i,j,k,m,subi;
  long Nt,Nstart,Nstop,subN,obs_id,nch;
  long acc;
  long iter,TICid;
  int *index, burn_in;
  bounds limited[Nparams], limits[Nparams];
  FILE *param_file, *data_file, *chain_file, *logL_file;
  char pfname[80],dfname[80],parname[80],subparname[80],
    outname[80],chainname[80],logLname[80];
  char RUN_ID[11];
  //characteristic magnitude of MCMC step in each parameter
  double  *sigma;

  Niter = (long)atoi(argv[1]);
  strcpy(RUN_ID,argv[2]);
  true_err = (double)atof(argv[3]);
  burn_in = atoi(argv[4]);

  strcat(subparname,"subpar.");
  strcat(subparname,RUN_ID);
  strcat(subparname,".dat");
  strcat(parname,"par.");
  strcat(parname,RUN_ID);
  strcat(parname,".dat");
  strcat(chainname,"chain.");
  strcat(chainname,RUN_ID);
  strcat(chainname,".dat");
  strcat(logLname,"logL.");
  strcat(logLname,RUN_ID);
  strcat(logLname,".dat");
  strcat(outname,RUN_ID);
  strcat(outname,".out");
  printf("%s\n",parname);
  P_ = (double **)malloc(Nchains*sizeof(double));
  for(i=0;i<Nchains;i++) P_[i]=(double *)malloc(Nparams*sizeof(double));
  x = (double **)malloc(Nchains*sizeof(double));
  for(i=0;i<Nchains;i++) x[i]=(double *)malloc(Nparams*sizeof(double));
  
  history = (double ***)malloc(Nchains*sizeof(double));
  for(i=0;i<Nchains;i++) {
    history[i]=(double **)malloc(Npast*sizeof(double));
    for(j=0;j<Npast;j++)history[i][j]=(double *)malloc(Nparams*sizeof(double));
  }
  y = (double *)malloc(Nparams*sizeof(double));
  sigma = (double *)malloc(Nparams*sizeof(double));

  
  scale  = 1.0/((double)(Nparams));
  //true_err = 1.0;
  srand(Niter);

  //TESS data file
  long seed = Niter;
  strcpy(dfname,"");
  strcat(dfname,RUN_ID);
  strcat(dfname,".txt"); //complete data set
  //strcat(dfname,".bin"); //binned data

  //read in starting parameters near best solution

  param_file = fopen(parname,"r");
  printf("%s\n",parname);
  //initialize priors
  set_limits(limited,limits);
  
  //set up proposal distribution
  initialize_proposals(sigma, history);
  
  for (i=0;i<Nparams;i++) {
    fscanf(param_file,"%lf\n", &tmp);
    P_[0][i] = tmp;
    P_0[i] = tmp;

    //if no map file, draw from prior
    //P_[0][i]   = limits[i].lo + 0.5*(limits[i].hi - limits[i].lo);//tmp;
		
    x[0][i] = P_[0][i];
    xmap[i] = P_[0][i];

    for(j=1; j<Nchains; j++) {
      P_[j][i]   = P_[0][i];
      x[j][i]  = P_[0][i];
    }
    if (burn_in == 1) {
      for(j=0; j<Nchains; j++) {
	P_[j][i]   = limits[i].lo + ran2(&seed)*(limits[i].hi - limits[i].lo);
	x[j][i]  = P_[j][i];
      }
    }
    if (burn_in == 2) {
      for(j=0; j<Nchains; j++) {
	P_[j][i]   = limits[i].lo + ran2(&seed)*(limits[i].hi - limits[i].lo);
	if (i == 2) P_[j][i] = P_0[i];
	x[j][i]  = P_[j][i];
      }
    }
    printf("%12.5e\n",tmp);
  }
  fclose(param_file);

  //read data file
  data_file = fopen(dfname,"r");
  printf("%s\n",dfname);
  tmp = 0;
  fscanf(data_file,"%ld\n", &Nt);
  printf("%ld\n",Nt);
  
  subN         = 0;
  rdata        = (double *)malloc(4*Nt*sizeof(double));
  index        = (int *)malloc(Nchains*sizeof(int));
  
  for (i=0;i<Nt;i++) {
    fscanf(data_file,"%lf %lf %lf %lf\n", &tmp1, &tmp2, &tmp3, &tmp4);
    rdata[i*4]=tmp1;
    rdata[i*4+1]=tmp2;
    rdata[i*4+2]=tmp3;
    rdata[i*4+3]=tmp4;
    if (tmp4 == 0) subN++;
    //printf("%12.5e %12.5e %12.5e %12.5e\n",tmp1,tmp2,tmp3,tmp4);
  }
  fclose(data_file);

  a_data       = (double *)malloc(subN*sizeof(double));
  a_model       = (double *)malloc(subN*sizeof(double));
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
      //e_data[subi] = rdata[i*4+2]; //binned data
      //printf("%ld %12.5e %12.5e %12.5e\n",subi,t_data[subi],a_data[subi],e_data[subi]);
      subi++;
    }
  }      
	
  //initialize parallel tempering
  dtemp = 1.2; //1.2;
  if (burn_in >= 1) dtemp = 1.3; //1.3;
  temp[0] = 1.0;
  for(i=1; i<Nchains; i++) {
    temp[i]  = temp[i-1]*dtemp;
    index[i] = i;
  }

  //initialize likelihood
  logLmap = loglikelihood(t_data,a_data,e_data,subN,P_[0]);
  printf("initial chi2 %g\n",-2*logLmap);
  for(i=0; i<Nchains; i++) logLx[i] = logLmap;

  //mcmc loop
  acc=0;
  chain_file = fopen(chainname,"w");
  logL_file  = fopen(logLname,"w");
  int atrial=0;
  int DEtrial=0;
  int DEacc=0;
  int jump;
  double jscale;
  for (iter=0;iter<Niter;iter++) {
    //simulated annealing
    //for(j=0; j<Nchains; j++) heat[j]=(1.0);// + 100000.0*exp(-(double)iter/1000.0));
    
    //loop over chains
    for(j=0; j<Nchains; j++) {
      alpha = ran2(&seed);
      //jscale = pow(10.,-8.+6.*alpha);
      jscale = pow(10.,-6.+6.*alpha);
      /* propose new solution */
      jump=0;
      if(ran2(&seed)<0.5 && iter>Npast) jump=1;
      //gaussian jumps along parameter directions
      if(jump==0)
	gaussian_proposal(x[index[j]], &seed, sigma, jscale, temp[j], y);
      //uniform draw on priors
      //uniform_proposal(x[index[j]], &seed, limits, y);
      
      //jump along correlations derived from chain history
      if(jump==1) {
	if(index[j]==0)DEtrial++;
	differential_evolution_proposal(x[index[j]], &seed, history[j], y);
	dx_mag=0;
	for (i=0;i<Nparams;i++) {
	  dx_mag+=(x[index[j]][i]-y[i])*(x[index[j]][i]-y[i]);
	}
	if (dx_mag < 1e-6)
	  gaussian_proposal(x[index[j]], &seed, sigma, jscale, temp[j], y);
      }
      
      for (i=0;i<Nparams;i++) {
	
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
      //printf("%d %d %12.5e %12.5e\n",j,index[j],logLy,logLx[index[j]]);
      //conditional acceptance of y
      if (alpha <= H)  {
	if(index[j]==0)acc++;
	for (i=0;i<Nparams;i++) {
	  x[index[j]][i]=y[i];
	}
	logLx[index[j]] = logLy;
	
	if(jump==1 && index[j]==0)DEacc++;
      }
      
      /* parallel tempering */
      ptmcmc(index,temp,logLx);
      
      /*  fill history  */
      k = iter - (iter/Npast)*Npast;
      for(i=0; i<Nparams; i++) history[j][k][i] = x[index[j]][i];
    }//end loop over chains
    
    //update map parameters
    if (logLx[index[0]] > logLmap) {
      for (i=0;i<Nparams;i++) {
	xmap[i]=x[index[0]][i];
      }
      logLmap = logLx[index[0]];
    }
    
    //down-sample chains until I get a better proposal
    if(iter%10==0) {
      //print parameter chains
      fprintf(chain_file,"%ld %.12g ",iter/10,logLx[index[0]]);
      for(i=0; i<Nparams; i++) fprintf(chain_file,"%.12g ",x[index[0]][i]);
      fprintf(chain_file,"\n");
      
      //print log likelihood chains
      fprintf(logL_file,"%ld ",iter/10);
      for(i=0; i<Nchains; i++) fprintf(logL_file,"%.12g ",logLx[index[i]]);
      fprintf(logL_file,"\n");
      
    }
    
    //update progress to screen
    if(iter%10==0) 
      {
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
      fclose(param_file);
    }
    
  }//end mcmc loop
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
  fclose(param_file);
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
	
	b = (int)((double)rand()/(RAND_MAX)*((double)(Nchains)))-1;
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
		//printf("%d %d %12.5e %12.5e\n",a,b,heat1,heat2);
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
	static long idum2=123456789;
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
	double dx[Nparams];
	
	//compute size of jumps
	for(n=0; n<Nparams; n++) dx[n] = ran2(seed)*(limits[n].hi - limits[n].lo);
	
	//uniform draw on prior range
	for(n=0; n<Nparams; n++) y[n] = limits[n].lo + dx[n];
}

void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y)
{
	int n;
	double gamma;
	double sqtemp;
	double dx[Nparams];
	
	
	//scale jumps by temperature
	sqtemp = sqrt(temp);
	
	//compute size of jumps
	for(n=0; n<Nparams; n++) dx[n] = gasdev2(seed)*sigma[n]*sqtemp*scale;
	
	//jump in parameter directions scaled by dx
	for(n=0; n<Nparams; n++) {
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
	double dx[Nparams];
	
	//choose two samples from chain history
	a = ran2(seed)*Npast;
	b = a;
	while(b==a) b = ran2(seed)*Npast;
	
	//compute vector connecting two samples
	for(n=0; n<Nparams; n++) {
	  dx[n] = history[b][n] - history[a][n];
	}
	//Blocks?
	
	//90% of jumps use Gaussian distribution for jump size
	if(ran2(seed) < 0.9) for(n=0; n<Nparams; n++) dx[n] *= gasdev2(seed)*GAMMA;

	//jump along vector w/ appropriate scaling
	for(n=0; n<Nparams; n++) y[n] = x[n] + dx[n];
}

void initialize_proposals(double *sigma, double ***history)
{
	int n,i,j;
	double junk;
	double x[Nparams];
	FILE *hfile;
	
	
	/*********************/    
	/* Gaussian Proposal */
	/*********************/    
	
	//jump sizes are set by hand based on what works
	sigma[0]  = 1.0e-1;  //log M1 (Msun)
	sigma[1]  = 1.0e-1;  //log M2 (Msun)
	sigma[2]  = sigmaP;  //log P (days)
	sigma[3]  = 1.0e-1;  //e
	sigma[4]  = 1.0e+0;  //inc (deg)
	sigma[5]  = 1.0e+1;  //Omega (deg)
	sigma[6]  = 1.0e+1;  //omega0 (deg)
	sigma[7]  = 1.0e+0;  //T0 (day)
	sigma[8]  = 1.0e-1;  //log flux normalization
	sigma[9]  = 1.0e-2;  //log rr1 normalization
	sigma[10]  = 1.0e-2;  //log rr2 normalization
	//if (burn_in == 0) for(i=0; i<Nparams; i++) sigma[i] /= 100;
	
	//for(i=0; i<Nparams; i++) sigma[i] /= ((double)(Nparams));
		
	
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
	
	for(i=0; i<Npast; i++)
	{
		//read in parameter from history file
		fscanf(hfile,"%i %lg",&n,&junk);
		for(n=0; n<Nparams; n++)fscanf(hfile,"%lg",&x[n]);
		
		//copy parameters across all chains
		for(j=0; j<Nchains; j++) for(n=0; j<Nparams; n++) history[j][i][n] = x[n];
	}
	
	fclose(hfile);
		 */
}
