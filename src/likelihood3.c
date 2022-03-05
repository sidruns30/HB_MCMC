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
    double *sorted_arr;
    sorted_arr = (double *)malloc((end-begin)*sizeof(double));

    long Nt = end - begin;

    for (int i=0; i<Nt; i++) {sorted_arr[i] = arr[begin + i];}

    quickSort(sorted_arr, 0, Nt - 1); 

    int mid;
    if (Nt % 2 == 0) mid = (int) Nt/2;
    else mid = (int) Nt/2 + 1;
    
    double median = sorted_arr[mid];

    for (int i=0; i<Nt; i++) {arr[begin + i] -= median;}
    free(sorted_arr);
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
void traj(double t, double pars[], int i, double *X1, double *X2, double *Y1, double *Y2,
        double *Z1, double *Z2, double *rr, double *ff){

    double P, M1, M2, a, e, inc, Omega, omega0, T0;
    double Mtot, r, f;

    t = t*SEC_DAY;
    M1 = pow(10.,pars[0])*MSUN;
    M2 = pow(10.,pars[1])*MSUN;
    P = pow(10.,pars[2])*SEC_DAY;
    e = pars[3];
    inc = pars[4];
    Omega = pars[5];
    omega0 = pars[6];
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

    // Note that these are in cgs
    *X1 = XX*(M2/Mtot);
    *Y1 = YY*(M2/Mtot);
    *Z1 = ZZ*(M2/Mtot);
    *X2 = -XX*(M1/Mtot);
    *Y2 = -YY*(M1/Mtot);
    *Z2 = -ZZ*(M1/Mtot);
    *rr = r;
    *ff = f;
    return;
}

/*
Function to calculate alpha beam
Values taken from Fig 5 (Claret et. al 2020: Doppler beaming factors for white dwarfs
                        main sequence stars and giant stars)
Using values for g=5
Input: log Temperature in Kelvin [NOT NORMALIZED BY SOLAR TEMP]
*/
double get_alpha_beam(double logT){

    // Initializing the alpha and temperature values
    double alphas[4] = {6.5, 4.0, 2.5, 1.2};
    double logTs[4] = {3.5, 3.7, 3.9, 4.5};

    // Return endpoints if temperature is outside the domain
    if (logT > logTs[3]) return 1.2/4;
    if (logT < logTs[0]) return 6.5/4;

    int j = 4;
    while(logT < logTs[j]) j--;

    return (alphas[j+1] + (alphas[j+1] - alphas[j]) / (logTs[j+1] - logTs[j]) * (logT - logTs[j+1]))/4;
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
double eclipse_area(double R1, double R2, 
                double X1, double X2, double Y1, double Y2){

    double h, r, dc, area, h_sq;
    // Overlap function borrowed from likelihood2.c
    double d = sqrt(SQR(X2-X1) + SQR(Y2-Y1))/RSUN;
    
    // Function is call by value so values of R1 and R2 aren't swapped in main code
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
    double Omega = pars[5];
    double omega0 = pars[6];
    double T0 = pars[7]*SEC_DAY;
    double rr1 = pars[8];
    double rr2 = pars[9];
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
    
    if (ALPHA_FREE == 1){
        // Limb and gravity darkening coefficients respectively
        mu_1 = pars[10];
        tau_1 = pars[11];
        mu_2 = pars[12];
        tau_2 = pars[13];
        // Reflection coefficients
        alpha_ref_1 = pars[14];
        alpha_ref_2 = pars[15];
	if (ALPHA_MORE ==1){
	  //extra alphas
	  extra_alpha_beam_1 = exp(pars[16]);//pow(10., pars[16]);
	  extra_alpha_beam_1 = exp(pars[17]);//pow(10., pars[17]);
	  alpha_Teff_1 = pars[18];//pow(10., pars[18]);
	  alpha_Teff_2 = pars[19];//pow(10., pars[19]);

        if (BLENDING == 1){
        blending = pars[20];
        flux_tune = pars[21];
        }
	}


    }
    else{
        mu_1 = .16;
        tau_1 = .344;
        mu_2 = .16;
        tau_2 = .344;
        alpha_ref_1 = .1;
        alpha_ref_2 = .1;
    }
    
    double M1 = pow(10., logM1);
    double M2 = pow(10., logM2);

    // Compute effective temperature and radius 
    double R1 = 0., R2 = 0., Teff1 = 0., Teff2 = 0.;

    R1 = pow(10., _getR(logM1) + rr1*envelope_Radius(logM1)); 
    R2 = pow(10., _getR(logM2) + rr2*envelope_Radius(logM2)); 

    Teff1 = pow(10., _getT(logM1) + alpha_Teff_1*envelope_Temp(logM1));
    Teff2 = pow(10., _getT(logM2) + alpha_Teff_2*envelope_Temp(logM2));

    //printf("Radii and temperatures are %f \t %f \t %f \t %f \t %f \t %f \n", 
    //R1, R2, Teff1, Teff2, _getT(logM1), _getT(logM2));

    // Temperature and radii are now in Rsun and K respectively

    // Flux normalization coefficients
    double Norm1, Norm2;
    Norm1 = SQR(R1) * QUAD(Teff1) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));
    Norm2 = SQR(R2) * QUAD(Teff2) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));

    // Set alpha_beam
    if (compute_alpha_beam == 1) {
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

    // Write lightcurve to a file
    
    /*FILE *lc_file;
    char* lc_file_name = "output_lc.txt";
    lc_file = fopen(lc_file_name,"w");
    */

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

        double beam1 = beaming(Pdays, M1, M2, e, inc, omega0, nu_arr[i], alpha_beam_1);
        double ellip1 = ellipsoidal(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R1, ar, mu_1, tau_1);
        double ref1 = reflection(Pdays, M1, M2, e, inc, omega0, nu_arr[i], R2, alpha_ref_1);

        Amag1[i] = Norm1 * (1 + beam1 + ellip1 + ref1);

        double beam2 = beaming(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], alpha_beam_2);
        double ellip2 = ellipsoidal(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R2, ar, mu_2, tau_2);
        double ref2 = reflection(Pdays, M2, M1, e, inc, (omega0+PI), nu_arr[i], R1, alpha_ref_2);

        Amag2[i] = Norm2 * (1 + beam2 + ellip2 + ref2);

        // Eclipse contribution (delta F = F * (ecl area / tot area))
        double area = eclipse_area(R1, R2, X1_arr[i], X2_arr[i], Y1_arr[i], Y2_arr[i]);
        if (Z2_arr[i] > Z1_arr[i]) Amag2[i] -= area * Norm2 / (PI * SQR(R2));
        else if (Z2_arr[i] < Z1_arr[i]) Amag1[i] -= area * Norm1 / (PI * SQR(R1));

        // Full lightcurve
        template[i] = (Amag1[i] + Amag2[i]);

        // If individual components are turned on
        if (SAVECOMP)
        {
            if (sizeof(template) < sizeof(11*Nt*sizeof(double)))
            {
                printf("template array not allocated enough memory to save all lightcurve components \n");
                return;
            }

            template[Nt + i] = Norm1 * beam1 + Norm2 * beam2;
            template[2*Nt + i] = Norm1 * ellip1 + Norm2 * ellip2;
            template[3*Nt + i] = Norm1 * ref1 + Norm2 * ref2;
            template[4*Nt + i] = 0.;
            
            if (Z2_arr[i] > Z1_arr[i]) template[4*Nt + i] -= area * Norm2 / (PI * SQR(R2));
            else if (Z2_arr[i] < Z1_arr[i]) template[4*Nt + i] -= area * Norm1 / (PI * SQR(R1));

            template[5*Nt + i] = X1;
            template[6*Nt + i] = Y1;
            template[7*Nt + i] = Z1;
            template[8*Nt + i] = X2;
            template[9*Nt + i] = Y2;
            template[10*Nt + i] = Z2;

        }

        // For writing onto the file
      /*  
        double full_beam, full_ellip, full_ref;
        full_beam = Norm1*beam1 + Norm2*beam2;
        full_ellip = Norm1*ellip1 + Norm2*ellip2;
        full_ref = Norm1*ref1 + Norm2*ref2;
        fprintf(lc_file,"%lf\t%lf\t%lf\t%lf\t%lf\n",times[i], template[i], full_beam, full_ellip, full_ref);
        */
    }

    //printf("Added individual components, now proceeding to remove the median... \n"); 

    // Normalize the lightcurve
    remove_median(template, 0, Nt);

    for (int i=0; i<Nt; i++) {
        template[i] += 1;
        template[i] = (1*blending + template[i]*(1 - blending)) * flux_tune;
    }

    // If components are turned on
    if (SAVECOMP)
    {
        if (sizeof(template) < sizeof(11*Nt*sizeof(double)))
        {
            printf("template array not allocated enough memory to save all lightcurve components \n");
            return;
        }

        //printf("Removing median 1 \n");
        remove_median(template, Nt, 2*Nt);
        //printf("Removing median 2 \n");
        remove_median(template, 2*Nt, 3*Nt);
        //printf("Removing median 3 \n");
        remove_median(template, 3*Nt, 4*Nt);
        //printf("Removing median 4 \n");
        remove_median(template, 4*Nt, 5*Nt);

        //printf("Adding the lightcurve components \n");
        for (int k=1; k<5; k++)
        {
            for (int i=0; i<Nt; i++)
            {
                template[k*Nt + i] += 1;
                template[k*Nt + i] = (1*blending + template[k*Nt + i]*(1 - blending)) * flux_tune;
            }
        }
    }

    //fclose(lc_file);
}

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

    double rr1 = params[8];
    double rr2 = params[9];

    double alpha_Teff_1 = 0.;
    double alpha_Teff_2 = 0.;

    if (ALPHA_MORE)
    {
        alpha_Teff_1 = params[18];
        alpha_Teff_2 = params[19];
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
    blending = params[20];
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
		     long N, double params[], double mag_data[], double magerr[],
             double weight)
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
      if (noise[i] < 1.e-5) 
      {
          printf("Old noise: %f\n" ,noise[i]);
          noise[i] = 1.e-5;
      }

      residual = (template[i]-lightcurve[i])/noise[i];
      chi2    += residual*residual;
    }
  //printf("chain chi2 is: %.10e\n", chi2);
  //free memory
  free(template);
  
  // Calculate magnitudes (Sets G, B-V, V-G and G-T mags)
  double D = mag_data[0];
  double Gmg, BminusV, VminusG, GminusT;
  
  calc_mags(params, D, &Gmg, &BminusV, &VminusG, &GminusT);
  double computed_mags[4] = {Gmg, BminusV, VminusG, GminusT};
  //printf("D, G, B-V, V-G, G-T:%f %f %f %f %f\n", magerr[0], magerr[1],mag_data[2],mag_data[3],mag_data[4]);
  //printf("G, B-V, V-G, G-T:%f %f %f %f \n", Gmg, BminusV, VminusG, GminusT);
  
  for (i=0;i<4;i++)
  {
    residual = (computed_mags[i] - mag_data[i+1])/magerr[i];
    chi2 += weight*residual*residual;
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
}


//Leave commented out unless for debugging purposes

/*
Function to write the lightcurve in a text file given the list of input parameters
for three periods
*/
/*
void write_lc_to_file(double pars[], char fname[])
{
    // Construct the time array
    const int N = 1000;
    double period = pow(10., pars[2]);
    double times[N]; 
    for (int i=0; i<N; i++)
    {
        times[i] = (double)i * (3 * period) / (double) N;
    }

    double *template;
    printf("Allocating memory \n");
    if (SAVECOMP)
    {
        template = (double *)malloc(5*N*sizeof(double));
        
    }

    else
    {
        template = (double *)malloc(N*sizeof(double));
    }

    printf("Size of template array is %d \n", sizeof(template));


    printf("Allocated memory, calculating lc \n");

    calc_light_curve(times, N, pars,template);

    printf("Calculated lc, moving on to saving \n");

    FILE *lcfile = fopen(fname, "w");

    for (int i=0; i<N; i++)
    {
        if (SAVECOMP)
        {
            fprintf(lcfile, "%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\t%12.5e\n", 
            times[i], template[i], template[N+i], template[2*N+i], template[3*N+i], template[4*N+i]);
        }

        else
        {
            fprintf(lcfile, "%12.5e\t%12.5e\n", times[i], template[i]);
        }
    }

    fclose(lcfile);
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