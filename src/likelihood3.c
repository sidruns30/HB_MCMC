/*
New version of likelihood that uses results from Engel et al 2020 to
calculate the lightcurve. Equations are initally derived in Kopal 1958
chapter 4, refined by Morris and then by Engel. The model includes
tidal distortion, rotational flattening, reflection and eclipses.
The parameters and their units are:
M1:             Mass of star 1 (log Msun)
M2:             Mass of star 2 (log Msun)
P:              Period (log days)
e:              Eccentricity
inc:            Inclination (radians)
Omega:          Longitude of ascending node (radians)
omega0:         Argument of periapse (radians)
T0:             Phase of the lightcurve (days)
rr1:            Radius scaling factor for star 1; R1 ~ R1_i x 10^(rr1)
rr2:            Radius scaling factor for star 2
*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <gsl/gsl_rng.h>
#include "likelihood3.h"



/*
Auxiliary functions to remove median from an array
The first 3 functions are necessary for quicksort
Borrowed from 
Parameters:
    arr1        orginal array
    Nt          Number of points in array
*/
void swap(double* a, double* b) 
{ 
    double t = *a; 
    *a = *b; 
    *b = t; 
} 

/* This function takes last element as pivot, places 
the pivot element at its correct position in sorted 
array, and places all smaller (smaller than pivot) 
to left of pivot and all greater elements to right 
of pivot */
double partition (double arr[], int low, int high) 
{ 
    double pivot = arr[high]; // pivot 
    int i = (low - 1); // Index of smaller element and indicates the right position of pivot found so far
  
    for (int j = low; j <= high - 1; j++) 
    { 
        // If current element is smaller than the pivot 
        if (arr[j] < pivot) 
        { 
            i++; // increment index of smaller element 
            swap(&arr[i], &arr[j]); 
        } 
    } 
    swap(&arr[i + 1], &arr[high]); 
    return (i + 1); 
} 

/* The main function that implements QuickSort 
arr[] --> Array to be sorted, 
low --> Starting index, 
high --> Ending index */
void quickSort(double arr[], int low, int high) 
{ 
    if (low < high) 
    { 
        /* pi is partitioning index, arr[p] is now 
        at right place */
        int pi = partition(arr, low, high); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(arr, low, pi - 1); 
        quickSort(arr, pi + 1, high); 
    } 
} 
  

void remove_median(double *arr, long begin, long end)
{
    // First sort the orignal array
    double sorted_arr[end-begin];

    long Nt = end - begin;

    for (int i=0; i<Nt; i++) {sorted_arr[i] = arr[begin + i];}

    quickSort(sorted_arr, 0, Nt - 1); 

    int mid;
    if (Nt % 2 == 0) mid = (int) Nt/2;
    else mid = (int) Nt/2 + 1;
    
    double median = sorted_arr[mid];

    for (int i=0; i<Nt; i++) {arr[begin + i] -= median;}
    return;
}


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
void traj(double *times, double *traj_pars, double *d_arr, 
        double *Z1_arr, double *Z2_arr, double *rr_arr, double *ff_arr, int Nt){

    double M1_cgs = traj_pars[0];
    double M2_cgs = traj_pars[1];
    double P_cgs = traj_pars[2];
    double e = traj_pars[3];
    double inc = traj_pars[4];
    double omega0 = traj_pars[5];
    double T0_cgs = traj_pars[6];

    if (M2_cgs > M1_cgs)
    {
        swap(&M1_cgs, &M2_cgs);
    }

    double Mtot_cgs = M1_cgs + M2_cgs;
    double a = pow(G * Mtot_cgs * SQR(P_cgs) / SQR(2 * PI), 1./3.);

    //printf("*******Main traj input parameters***** \n");
    //printf("%f %f %f %f %f %f %f \n", M1_cgs, M2_cgs, P_cgs, e, inc, omega0, T0_cgs);
    
    for (int ind=0; ind<Nt; ind++)
    {
        double t_cgs = times[ind] * SEC_DAY;
        //printf("time (cgs) %f \n", t_cgs);

        double M = 2.*PI * (t_cgs - T0_cgs) / P_cgs; 
        M = fmod(M,2*PI);
        double EE = M;
        double sin_M = sin(M);

        if(sin_M != 0.0)    EE = M + 0.85*e*sin_M/fabs(sin_M);

        // Kepler's Equation
        for(int j=0; j<5; j++)    EE = EE - (EE-e*sin(EE)-M)/(1-e*cos(EE));

        //printf("M, EE, and sinM %f %f %f\n", M, EE, sin_M);
        //printf("a is %f \n", a);
        
        rr_arr[ind] = a * (1 - e * cos(EE));
        ff_arr[ind] = 2.* atan(sqrt((1. + e)/(1. - e))* tan(EE/2.));

        double cos_omega0_f = cos(omega0 + ff_arr[ind]);
        double sin_omega0_f = sin(omega0 + ff_arr[ind]);
        double cos_inc = cos(inc);
        double sin_inc = sin(inc);

        double ZZ = rr_arr[ind] * sin_omega0_f * sin_inc;
        double d_factor = sqrt(SQR(cos_omega0_f) + SQR(sin_omega0_f * cos_inc));

        d_arr[ind] =  rr_arr[ind] * d_factor;
        Z1_arr[ind] = ZZ * (M2_cgs / Mtot_cgs);
        Z2_arr[ind] = -ZZ * (M1_cgs / Mtot_cgs);

        //printf("r_arr, nu_arr, d_arr, ZZ, Z1_arr, Z2_arrs %f %f %f %f %f %f\n", 
        //rr_arr[ind], ff_arr[ind], d_arr[ind], ZZ, Z1_arr[ind], Z2_arr[ind]);
    }

    return;
}

/*
Function to calculate alpha beam
Values taken from Fig 5 (Claret et. al 2020: Doppler beaming factors for white dwarfs
                        main sequence stars and giant stars)
Using values for g=5
Input: log Temperature in Kelvin [NOT NORMALIZED BY SOLAR TEMP]
*/
double get_alpha_beam(double logT)
{
    // Initializing the alpha and temperature values
    double alphas[4] = {6.5, 4.0, 2.5, 1.2};
    double logTs[4] = {3.5, 3.7, 3.9, 4.5};

    // Return endpoints if temperature is outside the domain
    if (logT >= logTs[3]) return 1.2/4;
    if (logT < logTs[0]) return 6.5/4;

    int j = 3;
    while(logT < logTs[j]) j--;

    double alpha_beam = ((alphas[j+1] + (alphas[j+1] - alphas[j]) / (logTs[j+1] - logTs[j]) * (logT - logTs[j+1]))/4 );
    return alpha_beam;
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

    double beam = -2830. * alpha_beam * fac1 * fac2 * fac3 * fac4 * ppm;

    return beam;
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
                double omega0, double nu, double R1, double a, double mu, double tau){

    double alpha_11 = 15 * mu * (2 + tau) / (32 * (3 - mu));
    double alpha_21 = 3 * (15 + mu) * (1 + tau) / (20 * (3 - mu));
    double alpha_2b1 = 15 * (1 - mu) * (3 + tau) / (64 * (3 - mu));
    double alpha_01 = alpha_21 / 9;
    double alpha_0b1 = 3 * alpha_2b1 / 20;
    double alpha_31 = 5 * alpha_11 / 3;
    double alpha_41 = 7*alpha_2b1 / 4;

    double beta = (1 + e * cos(nu)) / (1 - SQR(e));
    double q = M2 / M1;
    // Rotational velocity = Angular velocity at the periapse
    double Prot = P * pow(1 - e, 3./2);

    // "Mean ellipsoidal terms"
    double AM1, AM2, AM3;
    // Sin terms 
    double S1, S3;
    // Cosine terms
    double C2_1, C2_2, C4;

    double ppm = 1.e-6;
    // Flag to include higher order terms in flux
    int enable_const = 1;
    double tot_flux = 0.;

    // First the lowest order Engel Terms
    AM1 = 13435 * 2 * alpha_01 * (2 - 3*SQR(sin(inc))) * (1 / M1) * (1 / SQR(Prot)) * CUBE(R1);
    AM2 = 13435 * 3 * alpha_01 * (2 - 3*SQR(sin(inc))) * (1 / M1) * q / (1 + q) * (1 / SQR(P)) * CUBE(beta*R1);
    C2_1 = 13435 * alpha_21 * SQR(sin(inc)) * (1 / M1) * q / (1 + q) * (1 / SQR(P)) * CUBE(beta * R1) * cos(2*(omega0 + nu));
    
    tot_flux += (AM1 + AM2 + C2_1) * ppm;

    if (enable_const == 1){
        AM3 = 759 * alpha_0b1 * (8 - 40*SQR(sin(inc)) + 35*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3)
        * pow(P, -10./3) * pow(beta * R1, 5);
        S1 = 3194 * alpha_11 * (4 * sin(inc) - 5 * CUBE(sin(inc))) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3)
        * QUAD(beta * R1) * sin(omega0 + nu);
        C2_2 =  759 * alpha_2b1 * (6*SQR(sin(inc)) - 7*QUAD(sin(inc))) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3)
                * (beta * R1) * QUAD(beta * R1) * cos(2*(omega0 + nu));
        S3 = 3194 * alpha_31 * CUBE(sin(inc)) * pow(M1, -4./3) * q / pow(1+q, 4./3) * pow(P, -8./3) * QUAD(beta * R1) 
                * sin(3*(omega0 + nu));
        C4 = 759 * alpha_41 * QUAD(sin(inc)) * pow(M1, -5./3) * q / pow(1+q, 5./3) * pow(P, -10./3) * (beta * R1) *
                QUAD(beta * R1) * cos(4*(omega0 + nu));

        tot_flux += (AM3 + S1 + C2_2 + S3 + C4) * ppm;
    }
    

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
double eclipse_area(double R1, double R2, double d){

    double h, r, dc, area, h_sq;
    
    // Function is call by value so values of R1 and R2 aren't swapped in main code
    if (R2 > R1) {
        double temp_ = R1;
        R1 = R2;
        R2 = temp_;
    }
    
    area = 0.;
    d = fabs(d)/RSUN;
    dc = sqrt(R1*R1-R2*R2);
    // Now find the observed overlapping area between the two stars
    if (d >= (R1+R2)) area = 0.;
    if (d < (R1-R2)) area = PI*R2*R2;

    if ((d > dc)&(d < (R1+R2))) { 
        h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d);
        h = sqrt(h_sq);

        double Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h);
        double Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h);
        area = Arh1 + Arh2;
        }

    if ((d <= dc)&(d >= (R1-R2))) { 
        h_sq = (4.*d*d*R1*R1- SQR(d*d-R2*R2+R1*R1))/(4.*d*d);
        h = sqrt(h_sq);
        double Arh1 = R1*R1*asin(h/R1)-h*sqrt(R1*R1-h*h);
        double Arh2 = R2*R2*asin(h/R2)-h*sqrt(R2*R2-h*h);
        area = PI*R2*R2-(-Arh1 + Arh2);
        }

    return area;
}

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
        inc:        Inclination (rad)
        Omega:      Long of ascending node (rad)
        omega0:     Angle or periastron (rad)
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
    double inc = pars[4];
    double omega0 = pars[5];
    double T0 = pars[6]*SEC_DAY;
    double rr1 = pars[7];
    double rr2 = pars[8];
    double alpha_Teff_1 = 0.;
    double alpha_Teff_2 = 0.;
    
    // Beaming coefficients
    int compute_alpha_beam = 1;
    double alpha_beam_1 = 1.;
    double alpha_beam_2 = 1.;
    double extra_alpha_beam_1 = 1.;
    double extra_alpha_beam_2 = 1.;
    double mu_1, mu_2, tau_1, tau_2, alpha_ref_1, alpha_ref_2;
    double blending = 0.;
    double flux_tune = 1.;
    
    if (ALPHA_FREE == 1)
    {
        // Limb and gravity darkening coefficients respectively
        mu_1 = pars[9];
        tau_1 = pars[10];
        mu_2 = pars[11];
        tau_2 = pars[12];
        // Reflection coefficients
        alpha_ref_1 = pars[13];
        alpha_ref_2 = pars[14];
	if (ALPHA_MORE ==1)
    {
	  //extra alphas
	  extra_alpha_beam_1 = exp(pars[15]);//pow(10., pars[16]);
	  extra_alpha_beam_2 = exp(pars[16]);//pow(10., pars[17]);
	  alpha_Teff_1 = pars[17];//pow(10., pars[18]);
	  alpha_Teff_2 = pars[18];//pow(10., pars[19]);

        if (BLENDING == 1){
        blending = pars[19];
        flux_tune = pars[20];
        }
	}


    }
    
    else
    {
        mu_1 = .16;
        tau_1 = .344;
        mu_2 = .16;
        tau_2 = .344;
        alpha_ref_1 = .1;
        alpha_ref_2 = .1;
    }

    double M1 = pow(10., logM1);
    double M2 = pow(10., logM2);

    
    // Parameters for the trajectory function
    double M1_cgs = M1 * MSUN;
    double M2_cgs = M2 * MSUN;
    double P_cgs = P;
    double T0_cgs = T0;

    double traj_pars[7] = {M1_cgs, M2_cgs, P_cgs, e, inc, omega0, T0_cgs};

    // Compute effective temperature and radius 
    double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

    calc_radii_and_Teffs(pars, &R1, &R2, &Teff1, &Teff2);

    // Flux normalization coefficients
    double Norm1, Norm2;
    Norm1 = SQR(R1) * QUAD(Teff1) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));
    Norm2 = SQR(R2) * QUAD(Teff2) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));

    // Set alpha_beam
    if (compute_alpha_beam == 1) 
    {
        alpha_beam_1 = get_alpha_beam(log10(Teff1));
        alpha_beam_2 = get_alpha_beam(log10(Teff2));
    }

    alpha_beam_1 *= extra_alpha_beam_1;
    alpha_beam_2 *= extra_alpha_beam_2;

    // Semi majot axis (cgs calculation)
    double Mtot = (M1+M2)*MSUN;
    double a = pow(G*Mtot*P*P/(4.0*PI*PI),1./3.);
    double ar = a / RSUN;

    // Fluxes from the stars stored here
    double Amag1[Nt];
    double Amag2[Nt];

    // Positon arrays for the stars (cylindrical separaion)
    double d_arr[Nt];
    double Z1_arr[Nt];
    double Z2_arr[Nt];

    // Radial separation
    double r_arr[Nt];
    // True anomaly
    double nu_arr[Nt];

    traj(times, traj_pars, d_arr, Z1_arr, Z2_arr, r_arr, nu_arr, Nt);


    // Calculate trajectory and store results in arrays + unitialize fluxes to 0
    for (int i = 0; i<Nt; i++)
    {

        Amag1[i] = 0.;
        Amag2[i] = 0.;

        double beam1 = beaming(Pdays, M1, M2, e, inc, omega0, nu_arr[i], alpha_beam_1);
        double ellip1 = ellipsoidal(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R1, ar, mu_1, tau_1);
        double ref1 = reflection(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R2, alpha_ref_1);

        Amag1[i] = Norm1 * (1 + beam1 + ellip1 + ref1);

        double beam2 = beaming(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], alpha_beam_2);
        double ellip2 = ellipsoidal(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R2, ar, mu_2, tau_2);
        double ref2 = reflection(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R1, alpha_ref_2);

        Amag2[i] = Norm2 * (1 + beam2 + ellip2 + ref2);

        // Eclipse contribution (delta F = F * (ecl area / tot area))
        double area = eclipse_area(R1, R2, d_arr[i]);
        if (Z2_arr[i] > Z1_arr[i]) Amag2[i] -= area * Norm2 / (PI * SQR(R2));
        else if (Z2_arr[i] < Z1_arr[i]) Amag1[i] -= area * Norm1 / (PI * SQR(R1));

        // Full lightcurve
        template[i] = (Amag1[i] + Amag2[i]);

    }


    // Normalize the lightcurve
    remove_median(template, 0, Nt);

    for (int i=0; i<Nt; i++)
    {
        template[i] += 1;
        template[i] = (1*blending + template[i]*(1 - blending)) * flux_tune;
    }
}

/*
Compute the radius and Teff for each star
Radii are in solar units
Teffs in K
*/
void calc_radii_and_Teffs(double params[],  double *R1, double *R2, double *Teff1, double* Teff2) {
    // Calculate R1, R2, T1, T2 from the parameters
    // Compute effective temperature and radius 
    double logM1 = params[0];
    double logM2 = params[1];
    double rr1 = params[7];
    double rr2 = params[8];
    double alpha_Teff_1 = 0.;
    double alpha_Teff_2 = 0.;
    if (ALPHA_MORE)
    {
        alpha_Teff_1 = params[17];
        alpha_Teff_2 = params[18];
    }

    *R1 = pow(10., _getR(logM1) + rr1*envelope_Radius(logM1)); 
    *R2 = pow(10., _getR(logM2) + rr2*envelope_Radius(logM2)); 

    *Teff1 = pow(10., _getT(logM1) + alpha_Teff_1*envelope_Temp(logM1));
    *Teff2 = pow(10., _getT(logM2) + alpha_Teff_2*envelope_Temp(logM2));
    //printf("logM1,alpha_Teff1: %f, %f\n", logM1,alpha_Teff_1);
    //printf("logM2,alpha_Teff2: %f, %f\n", logM2,alpha_Teff_2);
    //printf("getT,envTemp: %f, %f\n", _getT(logM2),envelope_Temp(logM2));
    //printf("Teff1,Teff2: %f, %f\n", *Teff1,*Teff2);
};

/*
Star color calculator
Calculate apparent magnitudes B, G, V, and T given the stellar radii
R1,R2 in cm, T1,T2 in K, and D in pc. 
Sets G, B-V, V-G and G-T mags
*/
void calc_mags(double params[],  double D, double *Gmg, double *BminusV, 
double *VminusG, double *GminusT)
{

    double logM1 = params[0];
    double logM2 = params[1];

    double rr1 = params[7];
    double rr2 = params[8];

    double alpha_Teff_1 = 0.;
    double alpha_Teff_2 = 0.;

    if (ALPHA_MORE)
    {
        alpha_Teff_1 = params[17];
        alpha_Teff_2 = params[18];
    }

    // Calculate R1, R2, T1, T2 from the parameters
    // Compute effective temperature and radius 
    double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

    R1 = pow(10., _getR(logM1) + rr1*envelope_Radius(logM1)); 
    R2 = pow(10., _getR(logM2) + rr2*envelope_Radius(logM2)); 

    Teff1 = pow(10., _getT(logM1) + alpha_Teff_1*envelope_Temp(logM1));
    Teff2 = pow(10., _getT(logM2) + alpha_Teff_2*envelope_Temp(logM2));

    // [Units are Rsun and K respectively]

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

  if (ALPHA_MORE && BLENDING)
  {
    blending = params[19];
  }
  for (j=0;j<4;j++) {
    nu[j]=C/(lam[j]*1e-7);
    f_nu[j]=PI*(R1*R1*(2.*h*CUBE(nu[j])/SQR(C)/(exp(h*nu[j]/(k*Teff1))-1.))+
        R2*R2*(2.*h*CUBE(nu[j])/SQR(C)/(exp(h*nu[j]/(k*Teff2))-1.)))
      /(SQR(D)*SQR(pc_cgs));
    // Correction from the blending
    f_nu[j] = f_nu[j] / (1 - blending);
  }
  double Bmag = -2.5*log10(f_nu[0])-48.6;
  double Vmag = -2.5*log10(f_nu[1])-48.6;
  double Gmag = -2.5*log10(f_nu[2])-48.6;
  double Tmag = -2.5*log10(f_nu[3])-48.6;

  // Return the differences
  *Gmg = Gmag;
  *BminusV = Bmag - Vmag;
  *VminusG = Vmag - Gmag;
  *GminusT = Gmag - Tmag;

  //printf("Gmag is %f\n", Gmag);
  //printf("Bmag is %f\n", Bmag);
  //printf("Vmag is %f\n", Vmag);
  //printf("Tmag is %f\n", Tmag);
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

double loglikelihood(double time[], double lightcurve[], double noise[],
		     long N, double params[], double mag_data[], double magerr[])
{
  double template[N];
  double residual;
  double chi2;
  long i;
	
  //compute template light curve
  calc_light_curve(time,N,params,template);

  //sum square of residual to get chi-squared
  chi2 = 0.;
  for (i=0;i<N;i++) 
    { // bound on the noise:
      if (noise[i] < 1.e-5) 
      {
          noise[i] = 1.e-5;
      }

      residual = (template[i]-lightcurve[i])/noise[i];
      chi2    += residual*residual;

    }

  if (USE_COLOR_INFO || USE_GMAG)
  {
    // Calculate magnitudes (Sets Gmag, B-V, V-G and G-T mags)
    double D = mag_data[0];
    double Gmg, BminusV, VminusG, GminusT;
    
    calc_mags(params, D, &Gmg, &BminusV, &VminusG, &GminusT);
    double computed_mags[4] = {Gmg, BminusV, VminusG, GminusT};

    // First entry in computed mags is the Gmag 
    if (USE_GMAG)
    {
        residual = (computed_mags[0] - mag_data[1])/magerr[0];
        chi2 += residual*residual;   
    }

    // The remaining entries in computed_mags are the colors
    if (USE_COLOR_INFO)
    {
        for (int i=1; i<4; i++)
        {
            residual = (computed_mags[i] - mag_data[i+1])/magerr[i];
            chi2 += residual*residual;
        }
    }

  }

  // Check for Roche Overflow
  int RocheOverFlowFlag = 0;
  RocheOverFlowFlag = RocheOverflow(params);

  if (RocheOverFlowFlag)
  {
    chi2 = BIG_NUM;
  }

  //return log likelihood
  return(-chi2/2.0);
}


/*
Function to write the lightcurve in a text file given the list of input parameters
for three periods
*/
void write_lc_to_file(double pars[], char fname[])
{
    // Construct the time array
    const int N = 10000;
    double period = pow(10., pars[2]);
    double Tmax = 30. + period;

    double dt = Tmax / (double) N;

    double times[N];
    times[0] = 0.; 
    for (int i=1; i<N; i++)
    {
        times[i] = times[i-1] + dt;
    }

    double *template;
    //printf("Allocating memory \n");
    if (SAVECOMP)
    {
        template = (double *)malloc(11*N*sizeof(double));
        
    }

    else
    {
        template = (double *)malloc(N*sizeof(double));
    }

    //printf("Size of template array is %d \n", sizeof(template));


    //printf("Allocated memory, calculating lc \n");

    calc_light_curve(times, N, pars,template);

    //printf("Calculated lc, moving on to saving \n");

    FILE *lcfile = fopen(fname, "w");

    for (int i=0; i<N; i++)
    {
        if (SAVECOMP)
        {
            fprintf(lcfile, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t", 
            times[i], template[i], template[N+i], template[2*N+i], template[3*N+i], template[4*N+i]);

            fprintf(lcfile, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n", 
            template[5*N+i], template[6*N+i], template[7*N+i], template[8*N+i], template[9*N+i],
            template[10*N+i]);

        }

        else
        {
            fprintf(lcfile, "%12.5e\t%12.5e\n", times[i], template[i]);
        }
    }

    fclose(lcfile);
    return;
}

// Function to monitor Roche Lobe overflow (Eggleton 1985)
// Returns Roche Lobe radius over the binary separation
double Eggleton_RL(double q)
{
    return 0.49 * pow(q, 2./3) / (0.6 * pow(q, 2./3) + log(1 + pow(q, 1./3)));
}

// Everything in cgs - I am comparing Roche Lobe to Binary Separation
// at periapse
// Returns 1 if Roche Lobe overflows...
int RocheOverflow(double *pars)
{   
    double M1 = pow(10., pars[0]) * MSUN;
    double M2 = pow(10., pars[1]) * MSUN;
    double q = M1 / M2;
    double period = pow(10., pars[2]) * SEC_DAY;
    double ecc = pars[3];
    double R1 = pow(10., _getR(pars[0]) + pars[7]*envelope_Radius(pars[0])) * RSUN; 
    double R2 = pow(10., _getR(pars[1]) + pars[8]*envelope_Radius(pars[1])) * RSUN; 
    double Bin_Sep = pow(G * (M1 + M2) * SQR(period) / (4.0*PI*PI), 1./3.);
    // The factor after comes from assuming eccentric orbits
    double RL1_over_sep = Eggleton_RL(q);
    double RL2_over_sep = Eggleton_RL(1/q);


    double R1_over_sep = R1 / (Bin_Sep * (1 - ecc));
    double R2_over_sep = R2 / (Bin_Sep * (1 - ecc));

    int flag = ((RL1_over_sep < R1_over_sep) || (RL2_over_sep < R2_over_sep)) ? 1 : 0;
    return flag;

}

// Structure for generating parameters
struct InputPars
{
    bounds limits;
    int flag;

};

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
  //limits on omega0, in rads
  limited[5].lo = 2;
  limits[5].lo = -PI;
  limited[5].hi = 2;
  limits[5].hi = PI;
  gauss_pars[5].flag = 0;
  //limits on T0, in MDJ-2450000
  limited[6].lo = 1;
  limits[6].lo = 0.;
  limited[6].hi = 1;
  limits[6].hi = LC_PERIOD;
  gauss_pars[6].flag = 0;
  //limits on log rr1, the scale factor for R1
  limited[7].lo = 1;
  limits[7].lo = -3.;
  limited[7].hi = 1;
  limits[7].hi = 3.;
  gauss_pars[7].flag = 1.;
  //limits on log rr2, the scale factor for R2
  limited[8].lo = 1;
  limits[8].lo = -3.;
  limited[8].hi = 1;
  limits[8].hi = 3.;
  gauss_pars[8].flag = 1.;
  if (ALPHA_FREE == 1){
    // Limits of the alpha_coefficients
    // limits on limb darkening coefficient for star 1
    limited[9].lo = 1;
    limits[9].lo = 0.12;
    limited[9].hi = 1;
    limits[9].hi = 0.20;
    gauss_pars[9].flag = 1;
    // limits on gravity darkening coefficient for star 1
    limited[10].lo = 1;
    limits[10].lo = 0.3;
    limited[10].hi = 1;
    limits[10].hi = 0.38;
    gauss_pars[10].flag = 1;
    // limits on limb darkening coefficient for star 2
    limited[11].lo = 1;
    limits[11].lo = 0.12;
    limited[11].hi = 1;
    limits[11].hi = 0.20;
    gauss_pars[11].flag = 1;
    // limits on gravity darkening coefficient for star 2
    limited[12].lo = 1;
    limits[12].lo = 0.3;
    limited[12].hi = 1;
    limits[12].hi = 0.38;
    gauss_pars[12].flag = 1;
    // limits on reflection coefficients on star 1
    limited[13].lo = 1;
    limits[13].lo = 0.8;
    limited[13].hi = 1;
    limits[13].hi = 1.2;
    gauss_pars[13].flag = 1;
    // limits on reflection coefficients on star 2
    limited[14].lo = 1;
    limits[14].lo = 0.8;
    limited[14].hi = 1;
    limits[14].hi = 1.2;
    gauss_pars[14].flag = 1;
    if (ALPHA_MORE == 1){
      // limits on extra (log) beaming coefficient for star 1
      limited[15].lo = 1;
      limits[15].lo = -0.1;
      limited[15].hi = 1;
      limits[15].hi = 0.1;
      gauss_pars[15].flag = 1;
      // limits on extra (log) beaming coefficient for star 2
      limited[16].lo = 1;
      limits[16].lo = -0.1;
      limited[16].hi = 1;
      limits[16].hi = 0.1;
      gauss_pars[16].flag = 1;
      // limits on (log) Teff coefficient for star 1
      limited[17].lo = 1;
      limits[17].lo = -3.;
      limited[17].hi = 1;
      limits[17].hi = 3.;
      gauss_pars[17].flag = 1.;
      // limits on (log) Teff coefficient for star 2
      limited[18].lo = 1;
      limits[18].lo = -3.;
      limited[18].hi = 1;
      limits[18].hi = 3.;
      gauss_pars[18].flag = 1.;
      if (BLENDING == 1){
        // Blending coefficient in the flux
        limited[19].lo = 1;
        limits[19].lo = 0.;
        limited[19].hi = 1;
        limits[19].hi = 1.;
        gauss_pars[19].flag = 0;
        // FLux tune coefficient
        limited[20].lo = 1;
        limits[20].lo = 0.99;
        limited[20].hi = 1;
        limits[20].hi = 1.01;
        gauss_pars[20].flag = 0;
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
	sigma[5]  = 1.0e-3;  //omega0 (rad)
	sigma[6]  = 1.0e-5;  //T0 (day)
	sigma[7]  = 1.0e-1;  //log rr1 normalization
	sigma[8]  = 1.0e-1;  //log rr2 normalization
  sigma[9] = 1.0e-2;   // mu 1
  sigma[10] = 1.0e-2;   // tau 1
  sigma[11] = 1.0e-2;   // mu 2
  sigma[12] = 1.0e-2;   // tau 2
  sigma[13] = 1.0e-2;   // ref 1
  sigma[14] = 1.0e-2;   // ref 2
  sigma[15] = 1.0e-2;   // extra alph 1
  sigma[16] = 1.0e-2;   // extra alph 2
  sigma[17] = 1.0e-1;   // temp 1
  sigma[18] = 1.0e-1;   // temp 2
  sigma[19] = 1.0e-3;   // blending
  sigma[20] = 1.0e-5;   // flux tune
  
  // Use bigger sigmas if not using color info
  if ((!USE_COLOR_INFO) || (!USE_GMAG))
  {
      for (int k=0; k<NPARS; k++)
      {
          sigma[0] = 1.e-1; // logM1
          sigma[1] = 1.e-1; // logM2
          sigma[4] = 1.e-2; // inc
          sigma[5] = 1.e-2; // omega0
          sigma[6] = 1.e-3; // T0
          sigma[9] = 1.e-1; // mu 1
          sigma[10] = 1.e-1; // tau 1
          sigma[11] = 1.e-1; // mu 2
          sigma[12] = 1.e-1; // tau 2
          sigma[13] = 1.e-1; // ref 1
          sigma[14] = 1.e-1; // ref 2
          sigma[15] = 1.e-1; // extra alpha 1
          sigma[16] = 1.e-1; // extra alpha 2
          sigma[17] = 1.e-1; // temp 1
          sigma[18] = 1.e-1; // temp 2

      }
  }

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

//Leave commented out unless for debugging purposes

/*
int main()
{
    
    FILE *par_file, *lc_file;
    //char *par_fname = "synthetic_lc_parameters.txt";
    //par_file = fopen(par_fname, "w");
    //fprintf(par_file, "logM1\tlogM2\tlogP\te\tinc\tOmega\tomega\tT0\trr1\trr2\tblending\n");
    double more_pars[NPARS];
    for (int lc_id=0; lc_id<2; lc_id ++)
    {
    
        more_pars[0] =          .3;//-0.2 + 2*((double)rand()/RAND_MAX - 0.5) * 1;  logM1 Uniform in log space from .06 to 6.3 Msun
        more_pars[1] =          .3;//-0.2 + 2*((double)rand()/RAND_MAX - 0.5) * 1;  logM2
        more_pars[2] =           .84;//0.5 + 2*((double)rand()/RAND_MAX - 0.5) * .5; logP Uniform in log space from to 1. to 10 days
        more_pars[3] =           .58;//0*0.35 + 0*((double)rand()/RAND_MAX - 0.5) * .35; e* Unifrom between 0 and 0.7
        more_pars[4] =           .74;//PI/2 + 2*((double)rand()/RAND_MAX - 0.5) * .35; inc (rad) Uniform between 70 and 110 degrees
        more_pars[5] =           0.   + 2*((double)rand()/RAND_MAX - 0.5) * PI;    //Omega (rad) Doesn't really matter!
        more_pars[6] =           1.18 + PI;//0.   + 2*((double)rand()/RAND_MAX - 0.5) * PI;     omega0 (rad)  Uniform between -PI and PI
        more_pars[7] =           0.   + 2*((double)rand()/RAND_MAX - 0.5) * 1000.;     //T0 (JD) Uniform between -1000 and 1000
        more_pars[8] =           0.;//0.15 + 2*((double)rand()/RAND_MAX - 0.5) * .15;    //rr1  Uniform between 0 and factor of 2 in radius (logspace)
        more_pars[9] =           0.;//0.15 + 2*((double)rand()/RAND_MAX - 0.5) * .15;    //rr2  Uniform between 0 and factor of 2 in radius (logspace)
        more_pars[10] =           0.15;   /*mu 1
        more_pars[11] =           0.35;   /*tau 1
        more_pars[12] =           0.15;   /*mu 2
        more_pars[13] =           0.35;   /*tau 2 
        more_pars[14] =           0.5;    /*ref 1
        more_pars[15] =           0.5;    /*ref 2
        more_pars[16] =           0.;     /*beam 1
        more_pars[17] =           0.;     /*beam ]
        more_pars[18] =           0.;   /*Teff 1
        more_pars[19] =           0.;  /*Teff 2
        more_pars[20] =           0.2 + 2*((double)rand()/RAND_MAX - 0.5) * .2; /*blending  uniform between 0 and 0.4 
        more_pars[21] =           1.;      /*Flux tune
        //for (int l=0;l<10;l++)
        //{
        //    fprintf(par_file, "%f\t", more_pars[l]);
        //}
        //fprintf(par_file, "%f\n", more_pars[20]);
        char lc_fname[100];
        sprintf(lc_fname, "thompson_lc.txt");//"/scratch/ssolanski/lc_samples/%06d.txt", lc_id);
        write_lc_to_file(more_pars, lc_fname);
    }
    //fclose(par_file);
    
    /*
    char *fname = "synthetic_lc.txt";
    write_lc_to_file(john_pars, fname2);
    printf("lcs written to files ");
    double tmp1, tmp2, tmp3, tmp4;
    calc_mags(more_pars, 2011., &tmp1, &tmp2, &tmp3, &tmp4);
    
    return 1;
    
}
*/
/*
//Leave commented out unless for debugging purposes
int main()
{
    double john_pars[NPARS] = {-0.872961250467939, -1.40912878790312, 0.395192, 0.380400534734164,
                               1.49943375380441, 1.98966931388, -0.264732952318595, 1818.44300372128,
                               0.467461878906405, 0.943673457666681, 0.159138696038522,
                               0.385436121432153, 0.156860600683433, 0.334831704749105,
                               0.832664900136882, 0.655167012105327, -0.0236897653984011,
                               0.148295531143382, -0.107481257137963, 0.1464336480184, 
                               0.603298776112969, 1.00037117832857};
    double siddhant_pars[NPARS] = {-1.12401733053, -0.798865215346, 0.395192, 0.412233604282,
                                -1.43629596554, 1.98966931388, 2.64545773095, 52.072939571,
                                 0.529544333685, 0.722604904589,  0.159138696038522,
                               0.385436121432153, 0.156860600683433, 0.334831704749105,
                               0.832664900136882, 0.655167012105327, -0.0236897653984011,
                               0.148295531143382, -0.107481257137963, 0.1464336480184, 
                               0.603298776112969, 1.00037117832857};
    char *fname1 = "generated_lightcurve.txt";
    char *fname2 = "generated_lightcurve_john.txt";
    write_lc_to_file(john_pars, fname2);
    write_lc_to_file(siddhant_pars, fname1);
    return 1;
    
}
*/