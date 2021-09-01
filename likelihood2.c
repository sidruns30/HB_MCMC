/*********************************************************/
//similar to likelihood.c, but with independently varying R1,R2
//instead of using main sequence scaling relations
/*********************************************************/
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define pi 3.14159265358979323846
#define G 6.6743e-8 //cgs
#define c 2.998e10
#define AU 1.496e13
#define Msun 1.9885e33
#define Rsun 6.955e10
#define sec_day 86400.0

/*********************************************************/
double A_rh(double R, double h)
{
  double area;
  area = R*R*asin(h/R)-h*sqrt(R*R-h*h);
  return area;
}
/*********************************************************/
double overlap(double r1, double r2, double d)
{
  double h,r,dc,area;
  if (r2 > r1) {
    r = r1;
    r1 = r2;
    r2 = r;
  }
  d = fabs(d);
  if (d > (r1+r2)) area = 0.;
  if (d < (r1-r2)) area = pi*r2*r2;
  dc = sqrt(r1*r1-r2*r2);
  if ((d > dc)&(d < (r1+r2))) {
    h = sqrt((4.*d*d*r1*r1-pow((d*d-r2*r2+r1*r1),2.))/(4.*d*d));
    area = A_rh(r1,h)+A_rh(r2,h);
  }
  if ((d <= dc)&(d >= (r1-r2))) {
    //     d = dc+(d-r2)
    h = sqrt((4.*d*d*r1*r1-pow((d*d-r2*r2+r1*r1),2.))/(4.*d*d));    
    area = pi*r2*r2-(A_rh(r2,h)-A_rh(r1,h));
  }
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
  t = t*sec_day;
  M1 = pow(10.,pars[0])*Msun;
  M2 = pow(10.,pars[1])*Msun;
  P = pow(10.,pars[2])*sec_day;
  e = pars[3];
  inc = pars[4]*(pi/180.);
  Omega = pars[5]*(pi/180.);
  omega0 = pars[6]*(pi/180.);
  T0 = pars[7]*sec_day;
  Mtot = M1+M2;
  a = pow(G*Mtot*P*P/(4.0*pi*pi),1./3.);
  //P = sqrt(4.0*pi*pi*a*a*a/(G*Mtot));
  
  double M = 2.*pi * (t-T0)/P; 
  M = fmod(M,2*pi);
  double EE = M;
  if(sin(M) != 0.0){ 
  	 EE = M + 0.85*e*sin(M)/fabs(sin(M));
  }
  for(i=0;i<6;i++){
	EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE));
  }
  r = a*(1-e*cos(EE));
  f = 2.*atan(sqrt((1.+e)/(1.-e))*tan(EE/2.));

  double XX = r*(cos(Omega)*cos(omega0+f)-sin(Omega)*sin(omega0+f)*cos(inc));
  double YY = r*(sin(Omega)*cos(omega0+f)+cos(Omega)*sin(omega0+f)*cos(inc));
  double ZZ = r*sin(omega0+f)*sin(inc);
  *zdot = 1./(sqrt(1-e*e))*(cos(omega0+f)+e*cos(omega0));
  pos[0] = XX*(M2/Mtot);
  pos[1] = YY*(M2/Mtot);
  pos[2] = ZZ*(M2/Mtot);
  pos[3] = -XX*(M1/Mtot);
  pos[4] = -YY*(M1/Mtot);
  pos[5] = -ZZ*(M1/Mtot);
  
  //printf("t,T0,EE,a,zdot: %g,%g,%g,%g,%g\n",t,P,EE,a,*zdot);
  *rE = sqrt(4.*G*M1/(c*c)*fabs(ZZ)); //approx D_L = D_S
  //printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",G,M1,ZZ,a,r);
  //theta = acos(sin(ff+omega0)*sin(inc));
  *theta = omega0+f-pi/2.;
  *rr = r;
  *ff = f;
}


/***********************************************************/
// Guts of the main code, adapted for more general interface
// John Baker
void calc_light_curve(double *times, long Nt, double *pars, double *template)
{
  //times: input array of time sample points
  double R1 = 1.0; //units of Rsun
  double R2 = 1.0; //units of Rsun
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
  P = pow(10.,pars[2])*sec_day;
  e = pars[3];
  inc = pars[4]*(pi/180.);
  Omega = pars[5]*(pi/180.);
  omega0 = pars[6]*(pi/180.);
  T0 = pars[7]*sec_day;
  Flux_TESS = pars[8];
  Flux_TESS = pow(10.,Flux_TESS);
  rr1 = pow(10.,pars[9]);
  rr2 = pow(10.,pars[10]);
  Mtot = (M1+M2)*Msun;
  a = pow(G*Mtot*P*P/(4.0*pi*pi),1./3.);
  //P = sqrt(4.0*pi*pi*a*a*a/(G*Mtot));
  Pdays = P/sec_day;
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

  Flux1 = pi*R1*R1*pow(Teff1,4.);
  Flux2 = pi*R2*R2*pow(Teff2,4.);
  for (itime=0; itime<Nt; itime++){
    t=times[itime];
    traj(t,pars,pos,&zdot,&rE,&theta,&rr,&ff);
    rr = rr/Rsun; //units of Rsun
    //projected separation in units of Rsun
    d = sqrt((pos[3]-pos[0])*(pos[3]-pos[0])+(pos[4]-pos[1])*(pos[4]-pos[1]))/Rsun;
    Amag_limb=1.;

    //convert back to Agnieszka units
    Mtot = M1+M2;
    aR = a/Rsun;
    Adoppler = 2.8e-3*alphabeam*sin(inc)*pow(Pdays,-1./3)*pow(Mtot,-2./3.)*M1*zdot;
    Aellipse_phi = -alphaev*(M1/M2)*sin(inc)*sin(inc)*cos(2.*theta)*R2*R2*R2/(rr*rr*rr);
    Am = (1./9.)*alphaev*(2.+5.*M1/M2)*(2.-3.*sin(inc)*sin(inc))*R2*R2*R2/(aR*aR*aR);
    Aellipse_mean = Am*(e*cos(ff)*((3.+3.*e*cos(ff)+e*e*cos(ff)*cos(ff)))
			+3.*e-3.*e*e+e*e*e);
    Aellipse_mean = Aellipse_mean/((1.0-e*e)*(1.0-e*e)*(1.0-e*e));
    //printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
    //	   t,Pdays,Mtot,Adoppler,Aellipse_phi,Am,Aellipse_mean);
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

