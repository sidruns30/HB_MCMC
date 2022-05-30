#include "likelihood3.h"
#include "likelihood3_old.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <string.h>


int main(int argc, char* argv[])
{

    bounds limited[NPARS], limits[NPARS];
    gauss_bounds gauss_pars[NPARS];

    bounds_old limited_old[NPARS], limits_old[NPARS];
    gauss_bounds_old gauss_pars_old[NPARS];

    const double LC_DURATION = 10.;
    set_limits(limited, limits, gauss_pars, LC_DURATION);
    set_limits_old(limited_old, limits_old, gauss_pars_old, LC_DURATION);

    FILE *comparision;
    FILE *lc_generated;
    char *fname = "comparision.txt";


    comparision = fopen(fname, "w");

    double pars[NPARS];
    double pars_old[NPARS + 1];

    const int SIZE = 1000;
    const int NITER = 100000;

    for (int p=0; p<NITER; p++)
    {
        char lc_fname[200];
        sprintf(lc_fname, "lc_comp_%02d.txt", p);
        // lc_generated = fopen(lc_fname, "w");
        
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

            if (j <= 4)
            {
                pars_old[j] = pars[j];
            }

            else
            {
                pars_old[j+1] = pars[j];
            }
            
        }

        // Make the period at least 1 day
        int RL_flag = RocheOverflow(pars);

        while (RL_flag != 0)
        {

            pars[0] =  (limits[0].lo + (float) rand() / (float) RAND_MAX*(limits[0].hi - limits[0].lo));
            pars[1] =  (limits[1].lo + (float) rand() / (float) RAND_MAX*(limits[1].hi - limits[1].lo));
            pars[2] =  (limits[2].lo + (float) rand() / (float) RAND_MAX*(limits[2].hi - limits[2].lo));
            pars[3] =  (limits[3].lo + (float) rand() / (float) RAND_MAX*(limits[3].hi - limits[3].lo));
            RL_flag = RocheOverflow(pars);
        }

        pars_old[0] = pars[0];
        pars_old[1] = pars[1];
        pars_old[2] = pars[2];
        pars_old[3] = pars[3];

        // Random Omega value
        pars_old[5] = 0.5;

        double error = 0.;


        calc_light_curve(times, SIZE, pars, lc1);
        calc_light_curve_old(times, SIZE, pars_old, lc2);

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

            fprintf(comparision, "%lf\t%d\n", error, RL_flag);

        /*
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
        */
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