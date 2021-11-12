/*
New version of likelihood that uses results from Engel et al 2020 to
calculate the lightcurve. Equations are initally derived in Kopal 1958
chapter 4, refined by Morris and then by Engel. The model includes
tidal distortion, rotational flattening, reflection and eclipses.
*/

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
#define NPARS 11

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

/*
Trajectory calculator - recursively solve Kepler's equation using orbital
parameters to find the coordinates of the body from observer's frame.
Parameters:
    t:          time
    pars:       relevant orbital parameters (M1, M2, period, e, inc, Omega, omega0, T0, ...)
    i:          index of time array where time = t
    X1:         observed X position of star 1 (array)
    X2:         observed X position of star 2 (array)
    Y1:         observed Y position of star 1 (array)
    Y2:         observed Y position of star 2 (array)
    Z1:         observed Z position of star 1 (array)
    Z2:         observed Z position of star 2 (array)
    rr:         radial positon (array)
    ff:         true anomaly (array)
Note: Xi, Yi, Zi are in cgs units
*/
void traj(double t, double pars[], int i, double *X1, double *X2, double *Y1, double *Y2,
        double *Z1, double *Z2, double *rr, double *ff){

    double P, M1, M2, a, e, inc, Omega, omega0, T0;
    double Mtot, r, f;

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
    
    double M = 2.*PI * (t-T0)/P; 
    M = fmod(M,2*PI);
    double EE = M;
    double sin_M = sin(M);

    if(sin_M != 0.0)    EE = M + 0.85*e*sin_M/fabs(sin_M);

    // Kepler's Equation
    for(i=0;i<4;i++)    EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE));

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

    *X1 = XX*(M2/Mtot);
    *Y1 = YY*(M2/Mtot);
    *Z1 = ZZ*(M2/Mtot);
    *X2 = -XX*(M1/Mtot);
    *Y2 = -YY*(M1/Mtot);
    *Z2 = -ZZ*(M1/Mtot);
    *rr = r;
    *ff = f;

    //printf("X1: %f \t Y1 %f \t Z1 %f \t X2 %f \t Y2 %f \t Z2 %f \t r %f \t nu %f \n",
    //        *X1/RSUN, *Y1/RSUN, *Z1/RSUN, *X2/RSUN, *Y2/RSUN, *Z2/RSUN, *rr/RSUN, *ff);
    return;
}

/*
Beaming function. Taken from Engel et. al 2020
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    nu:         true anamoly (called 'f' in traj()) (rad)
    alpha_beam  Constant
*/
double beaming(double P, double M1, double M2, double e, double inc,
                double omega0, double nu, double alpha_beam){
    double q = M2/M1;
    double fac1 = q / pow(1 + q, 2/3);
    double fac2 = pow(M1, 1./3);
    double fac3 = pow(P, -1./3);
    double fac4 = sin(inc) * cos(omega0 + nu) / sqrt(1 - SQR(e));
    double ppm = 1.e-6;

    return -2830. * alpha_beam * fac1 * fac2 * fac3 * fac4 * ppm;
}



/*
Ellipsoidal variations function. Taken from Engel et. al 2020
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    nu:         true anamoly (called 'f' in traj()) (rad)
    R1:         radius of primary (Rsun)
    a:          semi-major axis (Rsun)
Note that the gravity (tau) and limb (mu) darkenening coefficients
are set in the function
*/
double ellipsoidal(double P, double M1, double M2, double e, double inc,
                double omega0, double nu, double R1, double a){

    double mu = 1.;
    double tau = 1.;

    double alpha_11 = 15 * mu * (2 + tau) / (32 * (32 - mu));
    double alpha_21 = 3 * (15 + mu) * (1 + tau) / (20 * (3 - mu));
    double alpha_2b1 = 15 * (1 - mu) * (3 + tau) / (64 * (3 - mu));
    double alpha_01 = alpha_21 / 9;
    double alpha_0b1 = 3 * alpha_2b1 / 20;
    double alpha_31 = 5 * alpha_11 / 3;
    double alpha_41 = 7*alpha_2b1 / 4;

    double beta = (1 + e * cos(nu)) / (1 - SQR(e));
    double q = M2 / M1;
    // Angular velocity at the periapse
    double Prot = P * pow(1 - e, 3./2);

    // Constant terms (can be turned off with the flag below)
    double CONS1, CONS2, CONS3;
    // Sin terms 
    double S1, S3;
    // Cosine terms
    double C2_1, C2_2, C4;

    double ppm = 1.e-6;
    // Flag to include const terms in flux
    int enable_const = 1;
    double tot_flux = 0.;

    if (enable_const == 1){
        double SMALL = 1.e-5;
        CONS1 = 13435 * 2 * alpha_01 * (2 - 3*SQR(sin(inc))) * CUBE(R1) / (M1 * SQR(Prot));
        CONS2 = 13435 * 2 * alpha_01 * (2 - 3*SQR(sin(inc))) * q / (1 - q + SMALL) * CUBE(beta*R1) / (M1 * SQR(Prot));
        CONS3 = 759 * alpha_0b1 * (8 - 40*SQR(sin(inc)) + 35*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3)
                * pow(P, -10./3) * pow(beta * R1, 5);

        tot_flux += (CONS1 + CONS2 + CONS3) * ppm;
    }

    S1 = 3194 * alpha_11 * (4 * sin(inc) - 5 * CUBE(sin(inc))) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3)
            * QUAD(beta * R1) * sin(omega0 + nu);
    C2_1 = 13435 * alpha_21 * SQR(sin(inc)) * q / (1 + q) * CUBE(beta * R1) / (M1 * SQR(P)) * cos(2*(omega0 + nu));
    C2_2 =  759 * alpha_2b1 * (6*SQR(sin(inc)) - 7*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3)
            * (beta * R1) * QUAD(beta * R1) * cos(2*(omega0 + nu));
    S3 = 3194 * alpha_31 * CUBE(sin(inc)) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3) * QUAD(beta * R1) 
            * sin(3*(omega0 + nu));
    C4 = 759 * alpha_41 * QUAD(sin(inc)) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3) * (beta * R1) *
            QUAD(beta * R1) * cos(4*(omega0 + nu));

    tot_flux += (S1 + C2_1 + C2_2 + S3 + C4) * ppm;
    
    //printf("Cons1 %f \t Cons2 %f \t Cons3 %f \t S1 %f \t C21 %f \t C22 %f \t S3 %f \t C4 %f \n", 
    //        CONS1, CONS2, CONS3, S1, C2_1, C2_2, S3, C4);
    return tot_flux; 
}

/*
Reflection function. Taken from Engel et. al (who refer to Faigler & Mazeh 2015)
Parameters:
    Period:     period (days)
    M1:         mass of primary (Msun)
    M2:         mass of secondary (Msun)
    e:          eccentricity
    inc:        inclination (rad)
    omega0:     argument of periastron (rad)
    nu:         true anamoly (called 'f' in traj()) (rad)
    R2:         radius of secondary (Rsun)
    alpha_ref1: free paramter for reflection
*/
double reflection(double P, double M1, double M2, double e, double inc, 
                    double omega0, double nu , double R2, double alpha_ref1){

    double q = M2 / M1;
    double beta = (1 + e * cos(nu)) / (1 - SQR(e));
    double ppm = 1.e-6;

    double fac1 = pow(1 + q, -2./3);
    double fac2 = pow(M1, -2./3);
    double fac3 = pow(P, -4./3);
    double fac4 = SQR(beta *R2);
    double fac5 = 0.64 - sin(inc) * sin(omega0 + nu) + 0.18 * SQR(sin(inc)) * (1 - cos(2*(omega0 + nu)));

    return 56514 * alpha_ref1 * fac1 * fac2 * fac3 * fac4 * fac5 * ppm;

}


/*
Eclipse function - taken from the old likelihood2.c file. Written by Jeremy (?)
Calcualtes the overlapping area of the two stars and fits the radius
Parameters:
    R1:           radius of star 1 (rsun)
    R2:           radius of star 2 (rsun)
    X1:           observed X position of star 1 (cgs)
    X2:           observed X position of star 2 (cgs)
    Y1:           observed Y position of star 1 (cgs)
    Y2:           observed Y position of star 2 (cgs)
IMPORTANT: Make sure to compare Z1 and Z2 before subtracting the relevant flux from the sources
Area returned is in units of solar radius^2. Make sure Xi, Yi, Zi and Ri all have same units!
*/
double eclipse_area(double R1, double R2, 
                double X1, double X2, double Y1, double Y2){

    double h, r, dc, area, h_sq;
    // Overlap function borrowed from likelihood2.c
    double d = sqrt(SQR(X2-X1) + SQR(Y2-Y1))/RSUN;
    
    if (R2 > R1) {
        double temp_ = R1;
        R1 = R2;
        R2 = temp_;
    }
    
    area = 0.;
    d = fabs(d);
    dc = sqrt(R1*R1-R2*R2);
    // Now find the observed overlapping area between the two stars
    if (d >= (R1+R2)) area = 0.;
    if (d < (R1-R2)) area = PI*R2*R2; // Shouldn't it be R1^2

    if ((d > dc)&(d < (R1+R2))) { 
        h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d);
        h = sqrt(h_sq);
        // Ask Jeremy about Arh
        double Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h);
        double Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h);
        area = Arh1 + Arh2;
        }

    if ((d <= dc)&(d >= (R1-R2))) { 
        h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d);
        h = sqrt(h_sq);
        double Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h);
        double Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h);
        area = PI*R2*R2-(Arh1 + Arh2);
        }

    return area;
}

/*
Full lightcurve calculation
Parameters:
    times:      Time array
    Nt:         Number of points in the time array
    pars:       Parameter array containing the following parameters in order:
        M1:         Mass of primary (log Msun)
        M2:         Mass of secondary (log Msun)
        P:          Period (log days)
        e:          Eccentricity
        inc:        Inclination (deg)
        Omega:      Long of ascending node (deg)
        omega0:     Angle or periastron (deg)
        T0:         Inital Time (days)
        Flux_TESS:  Flux scaling factor
        rr1:        Radius scaling factor of primary (log)
        rr2:        Radius scaling factor of secondary (log)
    template:   Array to store the lightcurve
Note that the reflection coefficients are set to 1 for now
*/
void calc_light_curve(double *times, long Nt, double *pars, double *template){
    
    // Extract the paramters
    double logM1 = pars[0];
    double logM2 = pars[1];
    // Period in seconds
    double P = pow(10., pars[2])*SEC_DAY;
    double Pdays = pow(10., pars[2]);
    double e = pars[3];
    double inc = pars[4]*(PI/180);
    double Omega = pars[5]*(PI/180);
    double omega0 = pars[6]*(PI/180);
    double T0 = pars[7]*SEC_DAY;
    double Flux_TESS = pow(10., pars[8]);
    double rr1 = pow(10., pars[9]);
    double rr2 = pow(10., pars[10]);

    double M1 = pow(10., logM1);
    double M2 = pow(10., logM2);

    // Compute effective temperature and radius
    double Tcoeff[] = {3.74677,0.557556,0.184408,-0.0640800,-0.0359547};
    double Rcoeff[] = {0.00158766,0.921233,-0.155659,-0.0739842,0.0581150};
    double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

    for (int j=0;j<=4;j++) {
        R1 += Rcoeff[j]*pow(logM1,j);
        Teff1 += Tcoeff[j]*pow(logM1,j);
        R2 += Rcoeff[j]*pow(logM2,j);
        Teff2 += Tcoeff[j]*pow(logM2,j);
    }
    // [Units are Rsun and K respectively]
    R1 = pow(10.,R1)*rr1;
    R2 = pow(10.,R2)*rr2;
    Teff1 = pow(10.,Teff1)/5580.;
    Teff2 = pow(10.,Teff2)/5580.;

    // Semi majot axis (cgs calculation)
    double Mtot = (M1+M2)*MSUN;
    double a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
    double ar = a / RSUN;

    // Fluxes from the stars stored here
    double Amag1[Nt];
    double Amag2[Nt];

    // Positon arrays for the stars
    double X1_arr[Nt];
    double X2_arr[Nt];
    double Y1_arr[Nt];
    double Y2_arr[Nt];
    double Z1_arr[Nt];
    double Z2_arr[Nt];

    // Radial separation
    double r_arr[Nt];
    // True anomaly
    double nu_arr[Nt];

    // Calculate trajectory and store results in arrays + unitialize fluxes to 0
    for (int i = 0; i<Nt; i++){
        double t = times[i];
        double X1, X2, Y1, Y2, Z1, Z2, radius, nu;
        traj(t, pars, i, &X1, &X2, &Y1, &Y2, &Z1, &Z2, &radius, &nu);
        // Set arrays
        X1_arr[i] = X1;
        Y1_arr[i] = Y1;
        Z1_arr[i] = Z1;
        X2_arr[i] = X2;
        Y2_arr[i] = Y2;
        Z2_arr[i] = Z2;
        r_arr[i] = radius;
        nu_arr[i] = nu;

        Amag1[i] = 0.;
        Amag2[i] = 0.;

        // Beaming and reflection coefficients
        double alpha_beam = 1.;
        double alpha_ref = 1.;

        double beam1 = beaming(Pdays, M1, M2, e, inc, omega0, nu_arr[i], alpha_beam);
        double ellip1 = ellipsoidal(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R1, ar);
        double ref1 = reflection(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R2, alpha_ref);

        Amag1[i] = beam1 + ellip1 + ref1;

        double beam2 = beaming(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], alpha_beam);
        double ellip2 = ellipsoidal(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R1, ar);
        double ref2 = reflection(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R2, alpha_ref);

        Amag2[i] = beam2 + ellip2 + ref2;

        // Eclipse contribution
        double area = eclipse_area(R1, R2, X1_arr[i], X2_arr[i], Y1_arr[i], Y2_arr[i]);
        if (Z2_arr[i] > Z1_arr[i]) Amag2[i] -= area*QUAD(Teff2);
        else if (Z2_arr[i] < Z1_arr[i]) Amag1[i] -= area*QUAD(Teff1);

        // Full lightcurve, scaled by both the stars
        template[i] = (Amag1[i] + Amag2[i]) * Flux_TESS;
        //printf("r value %f \t X1 value %f \t nu value %f \n", r_arr[i], X1_arr[i],  nu_arr[i]);
        //printf("beam1 %f \t ellip 1 %f \t ref 1 %f \t beam 2 %f \t ellip2 %f \t ref 2 %f \t",
        //        beam1, ellip1, ref1, beam2, ellip2, ref2);
        //if (Z2_arr[i] > Z1_arr[i]) printf("Ecl 1 %f \n ", area*QUAD(Teff2));
        //else if (Z2_arr[i] < Z1_arr[i]) printf("Ecl 2 %f \n", area*QUAD(Teff1));
        //else printf("\n");
    }   
}

/*
Likelihood Calculator
given an array of observed magnitudes a_data at times t_data, with  measurment errors
e_data, calculates the goodness-of-fit for the set of parameters P_
Parameters:
    time:   Time array
    data:   Lightcurve data
    noise:  Errors
    N:      Number of points in data
    params: Model parameters
 */

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
  chi2 = 0.;
  for (i=0;i<N;i++) 
    { // bound on the noise:
      if (noise[i] < 1.e-5) {
          printf("Old noise: %f\n" ,noise[i]);
          noise[i] = 1.e-5;}
      residual = (template[i]-data[i])/noise[i];
      chi2    += residual*residual;
    }
  //printf("chain chi2 is: %.10e\n", chi2);
  //free memory
  free(template);
  
  //return log likelihood
  return(-chi2/2.0);
}

