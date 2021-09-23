/*********************************************************/
//similar to likelihood.C, but with independently varying R1,R2
//instead of using main sequence scaling relations
/*********************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PI 3.14159265358979323846
#define G 6.6743e-8 //cgs
#define C 2.998e10
#define AU 1.496e13
#define MSUN 1.9885e33
#define RSUN 6.955e10
#define SEC_DAY 86400.0
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

static inline void swap(double *x, double *y){
  double temp = *x;
  *x = *y;
  *y = temp;
}

// Old functions
double A_rh_old(double R, double h);
double overlap_old(double r1, double r2, double d);
void traj_old(double t, double pars[], double pos[],
	  double *zdot, double *rE, double *theta,
	  double *rr, double *ff);
void calc_light_curve_old(double *times, long Nt, double *pars, double *template);
double loglikelihood_old(double time[], double data[], double noise[],
		     long N, double params[]);
// New functions
double A_rh(double R, double h);
double overlap(double r1, double r2, double d);
void traj(double t, double pars[], double pos[],
	  double *zdot, double *rE, double *theta,
	  double *rr, double *ff);
void calc_light_curve(double *times, long Nt, double *pars, double *template);
double loglikelihood(double time[], double data[], double noise[],
		     long N, double params[]);
void calc_light_curve_old_old(double *times, long Nt, double *pars, double *template);

/*testing*/
void test(){
    for (int i=0;i<10;i++){
        double R1 = (double)rand()/(RAND_MAX);
        double R2 = (double)rand()/(RAND_MAX);
        double R3 = (double)rand()/(RAND_MAX);
        double R4 = (double)rand()/(RAND_MAX);
        double R5 = (double)rand()/(RAND_MAX);
        double R6 = (double)rand()/(RAND_MAX);
        double R7 = (double)rand()/(RAND_MAX);
        double R8 = (double)rand()/(RAND_MAX);
        double R9 = (double)rand()/(RAND_MAX);
        double R10 = (double)rand()/(RAND_MAX);
        double R11 = (double)rand()/(RAND_MAX);
        double R12 = (double)rand()/(RAND_MAX);

        // Test A_rh:
        double f1_old = A_rh_old(R1,R2);
        double f1_new = A_rh(R1,R2);
        printf("Test (%d of 10) for A_rh: f_old - f_new = %f \n",i, (f1_old-f1_new));

        // Testing overlap
        double f2_old = overlap_old(R1, R2, R3);
        double f2_new = overlap(R1, R2, R3);
        printf("Test (%d of 10) for overlap: f_old - f_new = %f \n",i, (f2_old-f2_new));

        // Testing traj
        double pars_[] = {R1,R2,R3,R4,R5,R6,R7,R8};
        double pos_[] = {R1,R2,R3,R4,R5,R6};
        double pars[] = {R1,R2,R3,R4,R5,R6,R7,R8};
        double pos[] = {R1,R2,R3,R4,R5,R6};
        double zdot_ = R9, rE_ = R10, theta_ = R12, rr_ = R11, ff_ = R2;
        double zdot = R9, rE = R10, theta = R12, rr = R11, ff = R2;
        traj_old(R2, pars_, pos_, &zdot_, &rE_, &theta_, &rr_, &ff_);
        traj(R2, pars, pos, &zdot, &rE, &theta, &rr, &ff);
        printf("Test (%d out of 10) for traj: \n",i);
        printf("pos_old - pos_new = (%f, %f, %f, %f, %f, %f) \n", (pos_[0] - pos[0]), (pos_[1] - pos[1]),(pos_[2] - pos[2]),(pos_[3] - pos[3]),(pos_[4] - pos[4]),(pos_[5] - pos[5]));
        printf("zdot, rE, theta, rr, ff: (%f, %f, %f, %f, %f) \n", (zdot_-zdot), (rE_-rE), (theta_-theta), (rr_-rr), (ff_-ff));

        // Testing lightcurve
        double times_[] = {0.,0.1,0.2,0.5,0.9,10.,20.,100.};
        long Nt_ = 8;
        double pars1_[] = {R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11};
        double template_[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        double times[] = {0.,0.1,0.2,0.5,0.9,10.,20.,100.};
        long Nt = 8;
        double pars1[] = {R1,R2,R3,R4,R5,R6,R7,R8,R9,R10,R11};
        double template[] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
        calc_light_curve(times_,Nt_, pars1_, template_);
        calc_light_curve_old_old(times_,Nt, pars1_, template);
        printf("Test (%d out of 10) for calc_lc \n", i);
        printf("Delta lc values: \t");
        for (int j=0;j<10;j++){printf("%f, \t", (template_[j] - template[j]));}

        // End test
        printf(" \n ***************************************\n");
    }   
}

int main(){
    test();
    return 0;
}

/*********************************************************/
double A_rh_old(double R, double h)
{
  /*Commented stuff with extra overhead*/
  //double area;
  //area = R*R*asin(h/R)-h*sqrt(R*R-h*h);
  //return area;
  return R*R*asin(h/R)-h*sqrt(R*R-h*h);
}
/*********************************************************/
double overlap_old(double r1, double r2, double d)
{
  double h,r,dc,area,h_sq;
  if (r2 > r1) swap(&r1, &r2);
  d = fabs(d);
  if (d >= (r1+r2)) area = 0.;
  else if (d < (r1-r2)) area = PI*r2*r2;
  dc = sqrt(r1*r1-r2*r2);
  h_sq = (4.*d*d*r1*r1- SQR(d*d-r2*r2+r1*r1))
          /(4.*d*d);
  h = sqrt(h);
  if ((d > dc)&(d < (r1+r2))) { area = A_rh_old(r1,h)+A_rh_old(r2,h);}
  if ((d <= dc)&(d >= (r1-r2))) { area = PI*r2*r2-(A_rh_old(r2,h)-A_rh_old(r1,h));}
  return area;
}

/***********************************************************/
void traj_old(double t, double pars[], double pos[],
	  double *zdot, double *rE, double *theta,
	  double *rr, double *ff)
{
  int i;
  double P,M1,M2,a,e,inc,Omega,omega0,T0;
  double Mtot,r,f;
  t = t*SEC_DAY;
  M1 = pow(10.,pars[0])*MSUN;
  M2 = pow(10.,pars[1])*MSUN;
  P = pow(10.,pars[2])*SEC_DAY;
  e = pars[3];
  inc = pars[4]*(PI/180.);
  Omega = pars[5]*(PI/180.);
  omega0 = pars[6]*(PI/180.);
  T0 = pars[7]*SEC_DAY;
  Mtot = M1+M2;
  a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  //P = sqrt(4.0*PI*PI*a*a*a/(G*Mtot));

  // Check per likelihood time evaluation
  // Get rid of extra data in folded lc
  // Change parameters in lc
  // Fit a model and remove high std dev points
  double M = 2.*PI * (t-T0)/P; 
  M = fmod(M,2*PI);
  double EE = M;

  if(sin(M) != 0.0){ 
  	 EE = M + 0.85*e*sin(M)/fabs(sin(M));
  }
  for(i=0;i<6;i++){
	EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE));
  }
  r = a*(1-e*cos(EE));
  f = 2.*atan(sqrt((1.+e)/(1.-e))*tan(EE/2.));

  double cos_Omega = cos(Omega);
  double cos_omega0_f = cos(omega0+f);
  double sin_Omega = sin(Omega);
  double sin_omega0_f = sin(omega0+f);
  double cos_inc = cos(inc);
  double sin_inc = sin(inc);
  double cos_omega0 = cos(omega0);

  double XX = r*(cos_Omega*cos_omega0_f - sin_Omega*sin_omega0_f*cos_inc);
  double YY = r*(sin_Omega*cos_omega0_f + cos_Omega*sin_omega0_f*cos_inc);
  double ZZ = r*sin_omega0_f*sin_inc;

  *zdot = 1./(sqrt(1-e*e))*(cos_omega0_f + e*cos_omega0);
  pos[0] = XX*(M2/Mtot);
  pos[1] = YY*(M2/Mtot);
  pos[2] = ZZ*(M2/Mtot);
  pos[3] = -XX*(M1/Mtot);
  pos[4] = -YY*(M1/Mtot);
  pos[5] = -ZZ*(M1/Mtot);
  
  //printf("t,T0,EE,a,zdot: %g,%g,%g,%g,%g\n",t,P,EE,a,*zdot);
  *rE = sqrt(4.*G*M1/(C*C)*fabs(ZZ)); //approx D_L = D_S
  //printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",G,M1,ZZ,a,r);
  //theta = acos(sin(ff+omega0)*sin(inc));
  *theta = omega0+f-PI/2.;
  *rr = r;
  *ff = f;
}


/***********************************************************/
// Guts of the main code, adapted for more general interface
// John Baker
void calc_light_curve_old(double *times, long Nt, double *pars, double *template)
{
  //times: input array of time sample points
  double R1 = 1.0; //units of RSUN
  double R2 = 1.0; //units of RSUN
  double Teff1 = 1.0; //units of Tsun
  double Teff2 = 1.0; //units of Tsun
  double Flux1 = 1.0; //units of Fsun
  double Flux2 = 1.0; //units of Fsun
  double Flux_TESS = 1.0; //TESS flux units
  double alphabeam = 1.;
  double alphaev = 1.;
  double M1,M2,a,e,inc,Omega,omega0,T0;
  double t,zdot=0.,rE=0.,theta=0.,P,Pdays,rr=0.,ff=0.,Mtot,aR=0.;
  double Amag_limb;
  double Adoppler;
  double Aellipse_phi, Aellipse_mean;
  double result1, result2, result3;
  double lowerval = 0.0;
  double area,d,Am = 0.;
  double rr1,rr2;
  double Amag1[Nt], Amag2[Nt];
  double pos[] = {1,0,0,-1,0,0};
  double Tcoeff[] = {3.74677,0.557556,0.184408,-0.0640800,-0.0359547};
  double Rcoeff[] = {0.00158766,0.921233,-0.155659,-0.0739842,0.0581150};
  int j,itime=0;

  M1 = pow(10.,pars[0]);
  M2 = pow(10.,pars[1]);
  P = pow(10.,pars[2])*SEC_DAY;
  e = pars[3];
  inc = pars[4]*(PI/180.);
  Omega = pars[5]*(PI/180.);
  omega0 = pars[6]*(PI/180.);
  T0 = pars[7]*SEC_DAY;
  Flux_TESS = pars[8];
  Flux_TESS = pow(10.,Flux_TESS);
  rr1 = pow(10.,pars[9]);
  rr2 = pow(10.,pars[10]);
  Mtot = (M1+M2)*MSUN;
  a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  //P = sqrt(4.0*PI*PI*a*a*a/(G*Mtot));
  Pdays = P/SEC_DAY;
  //printf("%12.5e %12.5e %12.5e\n",a,P,Pdays);
  R1 = 0.0;
  Teff1 = 0.0;
  R2 = 0.0;
  Teff2 = 0.0;
  M1 = log10(M1);
  M2 = log10(M2);
  for (j=0;j<=4;j++) {
    R1 += Rcoeff[j]*pow(M1,j);
    Teff1 += Tcoeff[j]*pow(M1,j);
    R2 += Rcoeff[j]*pow(M2,j);
    Teff2 += Tcoeff[j]*pow(M2,j);
  }
  M1 = pow(10.,M1);
  M2 = pow(10.,M2);
  R1 = pow(10.,R1)*rr1;
  R2 = pow(10.,R2)*rr2;
  Teff1 = pow(10.,Teff1)/5580.;
  Teff2 = pow(10.,Teff2)/5580.;

  Flux1 = PI*R1*R1*pow(Teff1,4.);
  Flux2 = PI*R2*R2*pow(Teff2,4.);
  printf("Lightcurve variables are (M1,M2,R1,R2,Teff1,Teff): (%f,%f,%f,%f,%f,%f) \n",
        M1,M2,R1,R2,Teff1,Teff2);
  for (itime=0; itime<Nt; itime++){
    t=times[itime];
    traj_old(t,pars,pos,&zdot,&rE,&theta,&rr,&ff);
    rr = rr/RSUN; //units of RSUN
    //projected separation in units of RSUN
    d = sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+(pos[4]-pos[1])*(pos[4]-pos[1]))/RSUN;
    Amag_limb=1.;

    //convert back to Agnieszka units
    Mtot = M1+M2;
    aR = a/RSUN;
    Adoppler = 2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M1*zdot;
    Aellipse_phi = -alphaev*(M1/M2)*sin(inc)*sin(inc)*cos(2.*theta)*R2*R2*R2/(rr*rr*rr);
    Am = (1./9.)*alphaev*(2.+5.*M1/M2)*(2.-3.*sin(inc)*sin(inc))*R2*R2*R2/(aR*aR*aR);
    Aellipse_mean = Am*(e*cos(ff)*((3.+3.*e*cos(ff)+e*e*cos(ff)*cos(ff)))
			+3.*e-3.*e*e+e*e*e);
    printf("\t **Am, e, cos_ff, Acmean: (%f,%f,%f,%f) ** \n"
            ,Am, e, cos(ff), Aellipse_mean);
    Aellipse_mean = Aellipse_mean/((1.0-e*e)*(1.0-e*e)*(1.0-e*e));
    //printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
    //	   t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    Amag2[itime] = Flux2*(Amag_limb+Adoppler+1*Aellipse_phi+1*Aellipse_mean);
    if (pos[5] > pos[2]) {
      area = overlap_old(R1,R2,d);
      Amag2[itime]-=area*pow(Teff2,4.);
    }
    Adoppler = -2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M2*zdot;
    Aellipse_phi = -alphaev*(M2/M1)*sin(inc)*sin(inc)*cos(2.*theta)*R1*R1*R1/(rr*rr*rr);
    Am = (1./9.)*alphaev*(2.+5.*M2/M1)*(2.-3.*sin(inc)*sin(inc))*R1*R1*R1/(aR*aR*aR);
    Aellipse_mean = Am*(e*cos(ff)*((3.+3.*e*cos(ff)+e*e*cos(ff)*cos(ff)))
			+3.*e-3.*e*e+e*e*e);
    Aellipse_mean = Aellipse_mean/((1.0-e*e)*(1.0-e*e)*(1.0-e*e));
    //printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
    //	   t,zdot,rr,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    Amag1[itime] = Flux1*(Amag_limb+Adoppler+1*Aellipse_phi+1*Aellipse_mean);
    if (pos[5] < pos[2]) {
      area = overlap_old(R1,R2,d);
      Amag1[itime]-=area*pow(Teff1,4.);
    }
    template[itime]=(Amag1[itime]+Amag2[itime])*Flux_TESS;
  }
}

/*given an array of observed magnitudes a_data at times t_data, with  measurment errors
 e_data, calculates the goodness-of-fit for the set of parameters P_*/

double loglikelihood_old(double time[], double data[], double noise[],
		     long N, double params[])
{
  double *template;
  double residual;
  double chi2;
  long i;
	
  //allocate memory for light curve
  template = (double *)malloc(N*sizeof(double));
	
  //compute template light curve
  calc_light_curve_old(time,N,params,template);
	
  //sum square of residual to get chi-squared
  chi2 = 0;
  for (i=0;i<N;i++) 
    {
      residual = (template[i]-data[i])/noise[i];
      chi2    += residual*residual;
      //printf("%ld %12.5e %12.5e %12.5e\n",i,data[i],template[i],noise[i]);
    }
  //free memory
  free(template);
  
  //return log likelihood
  return(-chi2/2.0);
}

/*********************************************************************************************/
/*********************************************************/
//similar to likelihood.C, but with independently varying R1,R2
//instead of using main sequence scaling relations
/*********************************************************/

/*********************************************************/
double A_rh(double R, double h)
{
  return R*R*asin(h/R)-h*sqrt(R*R-h*h);
}
/*********************************************************/
double overlap(double r1, double r2, double d)
{
  double h,r,dc,area,h_sq;
  if (r2 > r1) swap(&r1, &r2);
  d = fabs(d);
  if (d >= (r1+r2)) area = 0.;
  else if (d < (r1-r2)) area = PI*r2*r2;
  dc = sqrt(r1*r1-r2*r2);
  h_sq = (4.*d*d*r1*r1- SQR(d*d-r2*r2+r1*r1))
          /(4.*d*d);
  h = sqrt(h);
  if ((d > dc)&(d < (r1+r2))) { area = A_rh(r1,h)+A_rh(r2,h);}
  if ((d <= dc)&(d >= (r1-r2))) { area = PI*r2*r2-(A_rh(r2,h)-A_rh(r1,h));}
  return area;
}

/***********************************************************/
void traj(double t, double pars[], double pos[],
	  double *zdot, double *rE, double *theta,
	  double *rr, double *ff)
{
  int i;
  double P,M1,M2,a,e,inc,Omega,omega0,T0;
  double Mtot,r,f;
  t = t*SEC_DAY;
  M1 = pow(10.,pars[0])*MSUN;
  M2 = pow(10.,pars[1])*MSUN;
  P = pow(10.,pars[2])*SEC_DAY;
  e = pars[3];
  inc = pars[4]*(PI/180.);
  Omega = pars[5]*(PI/180.);
  omega0 = pars[6]*(PI/180.);
  T0 = pars[7]*SEC_DAY;
  Mtot = M1+M2;
  a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  //P = sqrt(4.0*PI*PI*a*a*a/(G*Mtot));

  // Check per likelihood time evaluation
  // Get rid of extra data in folded lc
  // Change parameters in lc
  // Fit a model and remove high std dev points
  double M = 2.*PI * (t-T0)/P; 
  M = fmod(M,2*PI);
  double EE = M;
  double sin_M = sin(M);

  if(sin_M != 0.0){ 
  	 EE = M + 0.85*e*sin_M/fabs(sin_M);
  }
  for(i=0;i<6;i++){
	EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE));
  }
  r = a*(1-e*cos(EE));
  f = 2.*atan(sqrt((1.+e)/(1.-e))*tan(EE/2.));

  double cos_Omega = cos(Omega);
  double cos_omega0_f = cos(omega0+f);
  double sin_Omega = sin(Omega);
  double sin_omega0_f = sin(omega0+f);
  double cos_inc = cos(inc);
  double sin_inc = sin(inc);
  double cos_omega0 = cos(omega0);

  double XX = r*(cos_Omega*cos_omega0_f - sin_Omega*sin_omega0_f*cos_inc);
  double YY = r*(sin_Omega*cos_omega0_f + cos_Omega*sin_omega0_f*cos_inc);
  double ZZ = r*sin_omega0_f*sin_inc;

  *zdot = 1./(sqrt(1-e*e))*(cos_omega0_f + e*cos_omega0);
  pos[0] = XX*(M2/Mtot);
  pos[1] = YY*(M2/Mtot);
  pos[2] = ZZ*(M2/Mtot);
  pos[3] = -XX*(M1/Mtot);
  pos[4] = -YY*(M1/Mtot);
  pos[5] = -ZZ*(M1/Mtot);
  
  //printf("t,T0,EE,a,zdot: %g,%g,%g,%g,%g\n",t,P,EE,a,*zdot);
  *rE = sqrt(4.*G*M1/(C*C)*fabs(ZZ)); //approx D_L = D_S
  //printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",G,M1,ZZ,a,r);
  //theta = acos(sin(ff+omega0)*sin(inc));
  *theta = omega0+f-PI/2.;
  *rr = r;
  *ff = f;
}


/***********************************************************/
// Guts of the main code, adapted for more general interface
// John Baker
void calc_light_curve(double *times, long Nt, double *pars, double *template)
{
  //times: input array of time sample points
  double R1 = 1.0; //units of RSUN
  double R2 = 1.0; //units of RSUN
  double Teff1 = 1.0; //units of Tsun
  double Teff2 = 1.0; //units of Tsun
  double Flux1 = 1.0; //units of Fsun
  double Flux2 = 1.0; //units of Fsun
  double Flux_TESS = 1.0; //TESS flux units
  double alphabeam = 1.;
  double alphaev = 1.;
  double M1,M2,a,e,inc,Omega,omega0,T0;
  double t,zdot=0.,rE=0.,theta=0.,P,Pdays,rr=0.,ff=0.,Mtot,aR=0.;
  double Amag_limb;
  double Adoppler;
  double Aellipse_phi, Aellipse_mean;
  double result1, result2, result3;
  double lowerval = 0.0;
  double area,d,Am = 0.;
  double rr1,rr2;
  double Amag1[Nt], Amag2[Nt];
  double pos[] = {1,0,0,-1,0,0};
  double Tcoeff[] = {3.74677,0.557556,0.184408,-0.0640800,-0.0359547};
  double Rcoeff[] = {0.00158766,0.921233,-0.155659,-0.0739842,0.0581150};
  int j,itime=0;

  // This removes 4 log10 calls:
  double logM1 = pars[0];
  double logM2 = pars[1]; 
  R1 = 0.0;
  Teff1 = 0.0;
  R2 = 0.0;
  Teff2 = 0.0;
  for (j=0;j<=4;j++) {
    R1 += Rcoeff[j]*pow(logM1,j);
    Teff1 += Tcoeff[j]*pow(logM1,j);
    R2 += Rcoeff[j]*pow(logM2,j);
    Teff2 += Tcoeff[j]*pow(logM2,j);
  }

  M1 = pow(10.,pars[0]);
  M2 = pow(10.,pars[1]);
  P = pow(10.,pars[2])*SEC_DAY;
  e = pars[3];
  inc = pars[4]*(PI/180.);
  Omega = pars[5]*(PI/180.);
  omega0 = pars[6]*(PI/180.);
  T0 = pars[7]*SEC_DAY;
  Flux_TESS = pars[8];
  Flux_TESS = pow(10.,Flux_TESS);
  rr1 = pow(10.,pars[9]);
  rr2 = pow(10.,pars[10]);
  Mtot = (M1+M2)*MSUN;
  a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  //P = sqrt(4.0*PI*PI*a*a*a/(G*Mtot));
  Pdays = P/SEC_DAY;
  //printf("%12.5e %12.5e %12.5e\n",a,P,Pdays);
  R1 = pow(10.,R1)*rr1;
  R2 = pow(10.,R2)*rr2;
  Teff1 = pow(10.,Teff1)/5580.;
  Teff2 = pow(10.,Teff2)/5580.;

  Flux1 = PI*R1*R1*QUAD(Teff1);
  Flux2 = PI*R2*R2*QUAD(Teff2);

  /*Siddhant: taking out all the redundant function calls from the loop*/
  //convert back to Agnieszka units
  Mtot = M1+M2;
  aR = a/RSUN;
  Amag_limb=1.;

  double sin_inc = sin(inc);
  double iPdays_to_one_third = pow(Pdays, -1./3);
  double iMtot_to_two_thirds = pow(Mtot, -2./3);

  double cos_ff;
  double cos_2theta;
  printf("Lightcurve variables are (M1,M2,R1,R2,Teff1,Teff): (%f,%f,%f,%f,%f,%f) \n",
        M1,M2,R1,R2,Teff1,Teff2);
  for (itime=0; itime<Nt; itime++){
    t=times[itime];
    traj(t,pars,pos,&zdot,&rE,&theta,&rr,&ff);
    rr = rr/RSUN; //units of RSUN
    cos_ff = cos(ff);
    cos_2theta = cos(2.*theta);
    //projected separation in units of RSUN
    d = sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+(pos[4]-pos[1])*(pos[4]-pos[1]))/RSUN;
    /******************COMPUTE COEFFICIENTS************************/
    /* For Amag2 */
    Adoppler = 2.8e-3 * alphabeam * sin_inc * iPdays_to_one_third * iMtot_to_two_thirds * M1 * zdot;
    Aellipse_phi = -alphaev * (M1/M2) * SQR(sin_inc) * cos_2theta * CUBE(R2) / CUBE(rr);
    Am = (1./9.) * alphaev * (2.+5.*M1/M2) * (2.-3.* SQR(sin_inc)) *CUBE(R2) / CUBE(aR);
    Aellipse_mean = Am * (e*cos_ff * (3. + 3.*e*cos_ff + e*e*cos_ff*cos_ff)
			                    + 3.*e - 3.*e*e + e*e*e);
    //printf("\t **Am, e, cos_ff, Acmean: (%f,%f,%f,%f) ** \n"
    //,Am, e, cos(ff), Aellipse_mean);
    Aellipse_mean = Aellipse_mean/(CUBE(1.0-e*e));
    //printf("Cube test: %f \n", (CUBE(1.0-e*e)));
    //printf("\t **#2 Am, e, cos_ff, Acmean: (%f,%f,%f,%f) ** \n"
    //,Am, e, cos(ff), Aellipse_mean);
    //printf("New (1): t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean, %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n"
    //,t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    Amag2[itime] = Flux2 * (Amag_limb + Adoppler + Aellipse_phi + Aellipse_mean);
    
    /*Siddhant: why are we computing the same thing again? I might remove this*/
    /* For Amag 1*/
    Adoppler = -1 * Adoppler * M2 / M1; //-2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M2*zdot;
    Aellipse_phi = -alphaev * (M2/M1) * SQR(sin_inc) * cos_2theta * CUBE(R1) / CUBE(rr);
    Am = (1./9.) * alphaev * (2.+5.*M2/M1) * (2.-3.* SQR(sin_inc)) * CUBE(R1) / CUBE(aR);
    Aellipse_mean = Am * (e*cos_ff * (3. + 3.*e*cos_ff + e*e*cos_ff*cos_ff)
			                    + 3.*e - 3.*e*e + e*e*e);
    Aellipse_mean = Aellipse_mean/CUBE(1.0-e*e);
    Amag1[itime] = Flux1 * (Amag_limb+Adoppler + Aellipse_phi + Aellipse_mean);
    //printf("New (2): t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean, %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n"
    //,t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    /*Siddhant: what is this for? */
    area = overlap(R1,R2,d);
    if (pos[5] > pos[2]) {  Amag2[itime]-=area*QUAD(Teff2);}
    else if (pos[5] < pos[2]) { Amag1[itime]-=area*QUAD(Teff1);}

    template[itime]=(Amag1[itime]+Amag2[itime])*Flux_TESS;
  }
}

/*given an array of observed magnitudes a_data at times t_data, with  measurment errors
 e_data, calculates the goodness-of-fit for the set of parameters P_*/

double loglikelihood(double time[], double data[], double noise[],
		     long N, double params[])
{
  double *template;
  double residual;
  double chi2;
  long i;
	
  //allocate memory for light curve
  template = (double *)malloc(N*sizeof(double));
	
  //compute template light curve
  calc_light_curve(time,N,params,template);
	
  //sum square of residual to get chi-squared
  chi2 = 0;
  for (i=0;i<N;i++) 
    {
      residual = (template[i]-data[i])/noise[i];
      chi2    += residual*residual;
      //printf("%ld %12.5e %12.5e %12.5e\n",i,data[i],template[i],noise[i]);
    }
  //free memory
  free(template);
  
  //return log likelihood
  return(-chi2/2.0);
}


void calc_light_curve_old_old(double *times, long Nt, double *pars, double *template)
{
  //times: input array of time sample points
  double R1 = 1.0; //units of RSUN
  double R2 = 1.0; //units of RSUN
  double Teff1 = 1.0; //units of Tsun
  double Teff2 = 1.0; //units of Tsun
  double Flux1 = 1.0; //units of Fsun
  double Flux2 = 1.0; //units of Fsun
  double Flux_TESS = 1.0; //TESS flux units
  double alphabeam = 1.;
  double alphaev = 1.;
  double M1,M2,a,e,inc,Omega,omega0,T0;
  double t,zdot=0.,rE=0.,theta=0.,P,Pdays,rr=0.,ff=0.,Mtot,aR=0.;
  double Amag_limb;
  double Adoppler;
  double Aellipse_phi, Aellipse_mean;
  double result1, result2, result3;
  double lowerval = 0.0;
  double area,d,Am = 0.;
  double rr1,rr2;
  double Amag1[Nt], Amag2[Nt];
  double pos[] = {1,0,0,-1,0,0};
  double Tcoeff[] = {3.74677,0.557556,0.184408,-0.0640800,-0.0359547};
  double Rcoeff[] = {0.00158766,0.921233,-0.155659,-0.0739842,0.0581150};
  int j,itime=0;
  M1 = pow(10.,pars[0]);
  M2 = pow(10.,pars[1]);
  P = pow(10.,pars[2])*SEC_DAY;
  e = pars[3];
  inc = pars[4]*(PI/180.);
  Omega = pars[5]*(PI/180.);
  omega0 = pars[6]*(PI/180.);
  T0 = pars[7]*SEC_DAY;
  Flux_TESS = pars[8];
  Flux_TESS = pow(10.,Flux_TESS);
  rr1 = pow(10.,pars[9]);
  rr2 = pow(10.,pars[10]);
  Mtot = (M1+M2)*MSUN;
  a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
  //P = sqrt(4.0*PI*PI*a*a*a/(G*Mtot));
  Pdays = P/SEC_DAY;
  //printf("%12.5e %12.5e %12.5e\n",a,P,Pdays);
  R1 = 0.0;
  Teff1 = 0.0;
  R2 = 0.0;
  Teff2 = 0.0;
  M1 = log10(M1);
  M2 = log10(M2);
  for (j=0;j<=4;j++) {
    R1 += Rcoeff[j]*pow(M1,j);
    Teff1 += Tcoeff[j]*pow(M1,j);
    R2 += Rcoeff[j]*pow(M2,j);
    Teff2 += Tcoeff[j]*pow(M2,j);
  }
  M1 = pow(10.,M1);
  M2 = pow(10.,M2);
  R1 = pow(10.,R1)*rr1;
  R2 = pow(10.,R2)*rr2;
  Teff1 = pow(10.,Teff1)/5580.;
  Teff2 = pow(10.,Teff2)/5580.;
  Flux1 = PI*R1*R1*pow(Teff1,4.);
  Flux2 = PI*R2*R2*pow(Teff2,4.);
  printf("Old Lightcurve variables are (M1,M2,R1,R2,Teff1,Teff): (%f,%f,%f,%f,%f,%f) \n",
        M1,M2,R1,R2,Teff1,Teff2);
  for (itime=0; itime<Nt; itime++){
    t=times[itime];
    traj(t,pars,pos,&zdot,&rE,&theta,&rr,&ff);
    rr = rr/RSUN; //units of RSUN
    //projected separation in units of RSUN
    d = sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+(pos[4]-pos[1])*(pos[4]-pos[1]))/RSUN;
    Amag_limb=1.;
    //convert back to Agnieszka units
    Mtot = M1+M2;
    aR = a/RSUN;
    Adoppler = 2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M1*zdot;
    Aellipse_phi = -alphaev*(M1/M2)*sin(inc)*sin(inc)*cos(2.*theta)*R2*R2*R2/(rr*rr*rr);
    Am = (1./9.)*alphaev*(2.+5.*M1/M2)*(2.-3.*sin(inc)*sin(inc))*R2*R2*R2/(aR*aR*aR);
    Aellipse_mean = Am*(e*cos(ff)*((3.+3.*e*cos(ff)+e*e*cos(ff)*cos(ff)))
			+3.*e-3.*e*e+e*e*e);
    //printf("\t **Am, e, cos_ff, Acmean: (%f,%f,%f,%f) ** \n"
    //,Am, e, cos(ff), Aellipse_mean);
    //printf("2,3,4 cubed: (%f,%f,%f), \t (%f,%f,%f)", CUBE(2.), CUBE(3.), CUBE(4.)
    //                    ,2.*2.*2.,3.*3.*3.,4.*4.*4.);
    Aellipse_mean = Aellipse_mean/((1.0-e*e)*(1.0-e*e)*(1.0-e*e));
    //printf("Cube test: %f, \t %f \n", CUBE(1.0-e*e), (1.0-e*e)*(1.0-e*e)*(1.0-e*e));
    //printf("\t **# 2 Am, e, cos_ff, Acmean: (%f,%f,%f,%f) ** \n"
    //,Am, e, cos(ff), Aellipse_mean);
    //printf("Old (1): t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean, %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n"
    //	   ,t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    Amag2[itime] = Flux2*(Amag_limb+Adoppler+1*Aellipse_phi+1*Aellipse_mean);
    if (pos[5] > pos[2]) {
      area = overlap(R1,R2,d);
      Amag2[itime]-=area*pow(Teff2,4.);
    }
    Adoppler = -2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M2*zdot;
    Aellipse_phi = -alphaev*(M2/M1)*sin(inc)*sin(inc)*cos(2.*theta)*R1*R1*R1/(rr*rr*rr);
    Am = (1./9.)*alphaev*(2.+5.*M2/M1)*(2.-3.*sin(inc)*sin(inc))*R1*R1*R1/(aR*aR*aR);
    Aellipse_mean = Am*(e*cos(ff)*((3.+3.*e*cos(ff)+e*e*cos(ff)*cos(ff)))
			+3.*e-3.*e*e+e*e*e);
    Aellipse_mean = Aellipse_mean/((1.0-e*e)*(1.0-e*e)*(1.0-e*e));
    //printf("Old (2): t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean, %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n"
    //,t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    //printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
    //	   t,zdot,rr,Adoppler,Aellipse_phi,Am,Aellipse_mean);
    Amag1[itime] = Flux1*(Amag_limb+Adoppler+1*Aellipse_phi+1*Aellipse_mean);
    if (pos[5] < pos[2]) {
      area = overlap(R1,R2,d);
      Amag1[itime]-=area*pow(Teff1,4.);
    }
    template[itime]=(Amag1[itime]+Amag2[itime])*Flux_TESS;
  }
}