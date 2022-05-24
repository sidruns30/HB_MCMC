#include "likelihood3.h"
#include "likelihood3_old.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>

/* Set priors on parameters, and whether or not each parameter is bounded*/
/* Siddhant: Maybe just use an if condition/switch statment instead of limited*/
void set_limits(bounds limited[], bounds limits[], gauss_bounds gauss_pars[], double LC_PERIOD)
{
  //limits on M1, in log10 MSUN
  limited[0].lo = 1; 
  limits[0].lo = -1.5;
  limited[0].hi = 1;
  limits[0].hi = 1.5;
  gauss_pars[0].flag = 0;
  //limits on M2, in log10 MSUN
  limited[1].lo = 1;
  limits[1].lo = -1.5;
  limited[1].hi = 1;
  limits[1].hi = 1.5;
  gauss_pars[1].flag = 0;
  //limits on P, in log10 days
  limited[2].lo = 1;
  limits[2].lo = 0.;
  limited[2].hi = 1;
  limits[2].hi = 0.5;
  gauss_pars[2].flag = 0;
  //limits on e
  limited[3].lo = 1;
  limits[3].lo = 0.0;
  limited[3].hi = 1;
  limits[3].hi = 0.8;
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
  limits[8].lo = -2.;
  limited[8].hi = 1;
  limits[8].hi = 2.;
  gauss_pars[8].flag = 1.;
  //limits on log rr2, the scale factor for R2
  limited[9].lo = 1;
  limits[9].lo = -2.;
  limited[9].hi = 1;
  limits[9].hi = 2.;
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
        limits[20].hi = 0.8;
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

int main(int argc, char* argv[])
{

    bounds limited[NPARS], limits[NPARS];
    gauss_bounds gauss_pars[NPARS];
    const double LC_DURATION = 10.;
    set_limits(limited, limits, gauss_pars, LC_DURATION);

    FILE *comparision;
    FILE *lc_generated;
    char *fname = "comparision.txt";


    comparision = fopen(fname, "w");

    double pars[NPARS];

    const int SIZE = 1000;
    const int NITER = 100;

    for (int p=0; p<NITER; p++)
    {
        char lc_fname[200];
        sprintf(lc_fname, "lc_comp_%02d.txt", p);
        lc_generated = fopen(lc_fname, "w");
        
        double lc1[SIZE];
        double lc2[SIZE];

        double times[SIZE];

        clock_t start1, start2, end1, end2;
        double cpu_time_used1, cpu_time_used2;

        for (int i=0; i<SIZE; i++)
        {
            times[i] = (double) i * LC_DURATION / (double) SIZE;
        }
        
        for (int j=0; j<NPARS; j++)
        {
            double rand_float = (float) rand() / (float) RAND_MAX;
            pars[j] =  (limits[j].lo + rand_float*(limits[j].hi - limits[j].lo));
        }

        double error = 0.;

        calc_light_curve(times, SIZE, pars, lc1);
        calc_light_curve_old(times, SIZE, pars, lc2);

        for (int k=0; k<SIZE; k++)
        {
            error += fabs(lc1[k] - lc2[k]);
        }
        
        error /= SIZE;

        if (error > 0.1)
        {
            for (int j=0; j<NPARS; j++)
            {
                printf("Par %d is %lf \n", j, pars[j]);
            }
            printf("\n**********\n");
        }

            fprintf(comparision, "%lf\n", error);

        for (int i=0; i<NPARS; i++)
        {
            fprintf(lc_generated, "%lf\t", pars[i]);
        }

        fprintf(lc_generated, "\n");

        for (int i=0; i<SIZE; i++)
        {
            fprintf(lc_generated, "%lf\t%lf\n", lc1[i], lc2[i]);
        }
        
        fclose(lc_generated);
    }

    /*
    for (int test_id=0; test_id<100; test_id++)
    {
        const int SIZE = 100 +  10 * test_id;
        const int NITER = 100000;

        printf("TESTING SIZE %ld \n", SIZE);

        double lc1[SIZE];
        double lc2[SIZE];

        double times[SIZE];

        clock_t start1, start2, end1, end2;
        double cpu_time_used1, cpu_time_used2;

        for (int i=0; i<SIZE; i++)
        {
            times[i] = (double) i * (2 * pow(10., LC_PERIOD)) / (double) SIZE;
        }

        for (int j=0; j<NPARS; j++)
        {
                double rand_float = (float) rand() / (float) RAND_MAX;
                pars[j] =  (limits[j].lo + rand_float*(limits[j].hi - limits[j].lo));
        }

        start1 = clock();
        for (int iter=0; iter<NITER; iter++)
        {
            calc_light_curve(times, SIZE, pars, lc1);
        }
        end1 = clock();
        cpu_time_used1 = ((double) (end1 - start1)) / CLOCKS_PER_SEC;


        start2 = clock();
        for (int iter=0; iter<NITER; iter++)
        {

            calc_light_curve_old(times, SIZE, pars, lc2);
        }
        end2 = clock();
        cpu_time_used2 = ((double) (end2 - start2)) / CLOCKS_PER_SEC;

        double dt1 = (double) (end1 - start1);
        double dt2 = (double) (end2 - start2);

        fprintf(comparision, "%ld\t%lf\t%lf\t%lf\t%lf\n",SIZE, dt1, dt2, cpu_time_used1, cpu_time_used2);

    }
    */
    fclose(comparision);
    return 0;
}

            /*
            double error = 0.;
            for (int k=0; k<SIZE; k++)
            {
                error += fabs(lc1[k] - lc2[k]);
            }
            error /= SIZE;

            if (error > 0.1)
            {
                for (int j=0; j<NPARS; j++)
                {
                    printf("Par %d is %lf \n", j, pars[j]);
                }
                printf("\n**********\n");
            }

            fprintf(comparision, "%lf\n", error);
        }

        for (int i=0; i<SIZE; i++)
        {
            fprintf(lc_generated, "%lf\t%lf\n", lc1[i], lc2[i]);
        }
        */