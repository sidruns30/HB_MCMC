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

#define PI 3.14159265358979323846
#define G 6.6743e-8 //cgs
#define C 2.998e10
#define AU 1.496e13
#define MSUN 1.9885e33
#define RSUN 6.955e10
#define SEC_DAY 86400.0
#define ALPHA_FREE 1 // to set coefficitents as parameters in the model
#define ALPHA_MORE 0 // to add even more flexible coefficients
#if ALPHA_FREE == 1
  #if ALPHA_MORE == 1
    #define NPARS 20
  #else
    #define NPARS 16
  #endif
#else
    #define NPARS 10
#endif

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define QUAD(x) ((x)*(x)*(x)*(x))

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
  

void remove_median(double *arr, long Nt){
    // First sort the orignal array
    double *sorted_arr;
    sorted_arr = (double *)malloc(Nt*sizeof(double));

    for (int i=0; i<Nt; i++) {sorted_arr[i] = arr[i];}

    quickSort(sorted_arr, 0, Nt - 1); 

    int mid;
    if (Nt % 2 == 0) mid = (int) Nt/2;
    else mid = (int) Nt/2 + 1;
    
    double median = sorted_arr[mid];

    for (int i=0; i<Nt; i++) {arr[i] -= median;}
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
    if (logT > logTs[4]) return 1.2/4;
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
    double rr1 = pow(10., pars[8]);
    double rr2 = pow(10., pars[9]);
    double alpha_Teff_1 = 1.;
    double alpha_Teff_2 = 1.;
    
    // Beaming coefficients
    int compute_alpha_beam = 1;
    double alpha_beam_1 = 1.;
    double alpha_beam_2 = 1.;
    double extra_alpha_beam_1 = 1.;
    double extra_alpha_beam_2 = 1.;
    double mu_1, mu_2, tau_1, tau_2, alpha_ref_1, alpha_ref_2;
    
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
	  extra_alpha_beam_1 = pars[16];
	  extra_alpha_beam_1 = pars[17];
	  alpha_Teff_1 = pars[18];
	  alpha_Teff_2 = pars[19];
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
    Teff1 *= alpha_Teff_1;
    Teff2 *= alpha_Teff_2;

    // Flux normalization coefficients
    double Norm1, Norm2;
    Norm1 = SQR(R1) * QUAD(Teff1) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));
    Norm2 = SQR(R2) * QUAD(Teff2) / (SQR(R1) * QUAD(Teff1) + SQR(R2) * QUAD(Teff2));

    // Set alpha_beam
    if (compute_alpha_beam == 1) {
        alpha_beam_1 = get_alpha_beam(log10(Teff1 * 5580));
        alpha_beam_2 = get_alpha_beam(log10(Teff2 * 5580));
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

        // For writing onto the file
      /*  
        double full_beam, full_ellip, full_ref;
        full_beam = Norm1*beam1 + Norm2*beam2;
        full_ellip = Norm1*ellip1 + Norm2*ellip2;
        full_ref = Norm1*ref1 + Norm2*ref2;
        fprintf(lc_file,"%lf\t%lf\t%lf\t%lf\t%lf\n",times[i], template[i], full_beam, full_ellip, full_ref);
        */
    } 

    // Normalize the lightcurve
    remove_median(template, Nt);
    for (int i=0; i<Nt; i++) template[i] += 1;

    //fclose(lc_file);
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

typedef struct {
  double lo;
  double hi;
} bounds;

/* Set priors on parameters, and whether or not each parameter is bounded*/
/* Siddhant: Maybe just use an if condition/switch statment instead of limited*/
void set_limits(bounds limited[], bounds limits[])
{
  //limits on M1, in log10 MSUN
  limited[0].lo = 1; 
  limits[0].lo = -1.5;
  limited[0].hi = 1;
  limits[0].hi = 2.0;
  //limits on M2, in log10 MSUN
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
  //limits on inc, in rads
  limited[4].lo = 2;
  limits[4].lo = -PI/2;
  limited[4].hi = 2;
  limits[4].hi = PI/2;
  //limits on Omega, in rads
  limited[5].lo = 2;
  limits[5].lo = -PI;
  limited[5].hi = 2;
  limits[5].hi = PI;
  //limits on omega0, in rads
  limited[6].lo = 2;
  limits[6].lo = -PI;
  limited[6].hi = 2;
  limits[6].hi = PI;
  //limits on T0, in MDJ-2450000
  limited[7].lo = 1;
  limits[7].lo = -1000;
  limited[7].hi = 1;
  limits[7].hi = 1000.;
  //limits on log rr1, the scale factor for R1
  limited[8].lo = 1;
  limits[8].lo = 0;
  limited[8].hi = 1.;
  limits[8].hi = 1.0;
  //limits on log rr2, the scale factor for R2
  limited[9].lo = 1;
  limits[9].lo = 0.;
  limited[9].hi = 1.;
  limits[9].hi = 1.;
  if (ALPHA_FREE == 1){
    // Limits of the alpha_coefficients
    // limits on limb darkening coefficient for star 1
    limited[10].lo = 1;
    limits[10].lo = 0.;
    limited[10].hi = 1.;
    limits[10].hi = 1.;
    // limits on gravity darkening coefficient for star 1
    limited[11].lo = 1;
    limits[11].lo = 0.;
    limited[11].hi = 1.;
    limits[11].hi = 1.;
    // limits on limb darkening coefficient for star 2
    limited[12].lo = 1;
    limits[12].lo = 0.;
    limited[12].hi = 1.;
    limits[12].hi = 1.;
    // limits on gravity darkening coefficient for star 2
    limited[13].lo = 1;
    limits[13].lo = 0.;
    limited[13].hi = 1.;
    limits[13].hi = 1.;
    // limits on reflection coefficients on star 1
    limited[14].lo = 1;
    limits[14].lo = 0.;
    limited[14].hi = 1.;
    limits[14].hi = 1.;
    if (ALPHA_MORE == 1){
      // limits on extra beaming coefficient for star 1
      limited[16].lo = 1;
      limits[16].lo = 0.9;
      limited[16].hi = 1;
      limits[16].hi = 1.1;
      // limits on extra beaming coefficient for star 2
      limited[17].lo = 1;
      limits[17].lo = 0.9;
      limited[17].hi = 1;
      limits[17].hi = 1.1;
      // limits on Teff coefficient for star 1
      limited[18].lo = 1;
      limits[18].lo = 0.9;
      limited[18].hi = 1;
      limits[18].hi = 1.1;
      // limits on Teff coefficient for star 2
      limited[19].lo = 1;
      limits[19].lo = 0.9;
      limited[19].hi = 1;
      limits[19].hi = 1.1;
    }
  }
}
