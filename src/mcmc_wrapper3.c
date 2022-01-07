#include <math.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <unistd.h>
#include "likelihood3.c"
//#include "extra.c"

#define NCHAINS 50
#define NPAST 500 // Relevant for differential evolution
#define GAMMA (2.388 / sqrt(2. * NPARS))
#define NUM_ELEMENTS(x) (size_of(x) / size_of(x[0]))

/*************************Contain all the globabls****************************************/
struct globals
{
  // Input parameters
  int Niter;
  char *RUN_ID;
  double true_err, burn_in, guess_period;
  // Problem variables
  double SIGMAP, scale, logLmap;
  long subN, subi, Nt, acc;
  // Output file names
  char pfname[80], dfname[80], parname[80], subparname[80], outname[80], chainname[80], logLname[80];
  // File objects
  FILE *param_file, *data_file, *chain_file, *logL_file;
  // Data arrays
  int *index;
  double *a_data, *a_model, *e_data, *t_data, *q_data, *sigma, *light_curve, *light_curve_2, *rdata;
  double **P_, **x;
  double ***history;
  double xmap[NPARS];
  double heat[NCHAINS], temp[NCHAINS];
  double logLx[NCHAINS];
  // Random number stuff
  long seed;
  // Limits for parameters
  bounds limited[NPARS];
  bounds limits[NPARS];
  // For parallel tempering
  double dtemp;
};

/******************FUNCTION DECLARATIONS*************************************/
/* From likelihood.c */
void calc_light_curve(double t_data[], long Nt, double P_[], double light_curve[]);
double likelihood(double t_data[], double a_data[], double e_data[], long Nt, double P_[]);

/* From extra.c */
double ran2(long *idum);
double gasdev2(long *idum);
void free_1d(double *arr);
void free_2d(double **arr, int size);
void free_3d(double ***arr, int size1, int size2);
int file_exists(const char *fname);


void ptmcmc(int *index, double heat[], double logL[]);
void initialize_proposals(struct globals *glb);
void initialize_data(int argc, char *argv[], struct globals *glb);
void uniform_proposal(double *x, long *seed, bounds limits[], double *y);
void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y);
void hypercube_proposal(double *x, long *seed, double *y);
void differential_evolution_proposal(double *x, long *seed, double **history, double *y);
void run_chain(int chain_id, int iter, struct globals *glb);
void write_and_log(int iter, struct globals *glb);
/**************************************************************************/

int main(int argc, char *argv[])
{

  omp_set_num_threads(1);

  struct globals glb;
  initialize_data(argc, argv, &glb);

  // Main iteration loop
  for (int iter = 0; iter < glb.Niter; iter++)
  {
// Chain loop
//#pragma omp parallel for shared(glb)
    for (int chain_id = 0; chain_id < NCHAINS; chain_id++)
    {
      run_chain(chain_id, iter, &glb);
    }
    //#pragma omp barrier
    write_and_log(iter, &glb);
  }

  fclose(glb.chain_file);
  fclose(glb.logL_file);

  return 0;
}

void initialize_data(int argc, char *argv[], struct globals *glb)
{
  printf("Initialzing data...\n");

  glb->Niter = atoi(argv[1]);
  glb->RUN_ID = argv[2];
  glb->true_err = (double)atof(argv[3]);
  glb->burn_in = atoi(argv[4]);
  glb->guess_period = (double)atof(argv[5]);
  int run = atoi(argv[6]);

  // Set period jump
  if (glb->burn_in == 2)
  {
    glb->SIGMAP = 1.e-10;
  }
  else
  {
    glb->SIGMAP = 1.e-1;
  }

  printf("\t Setting filenames...\n");
  // Add in the logfile names
  //sprintf(glb->subparname, "%c", "");
  //sprintf(glb->parname, "%c", "");
  //sprintf(glb->chainname, "%c", "");
  //sprintf(glb->logLname, "%c", "");
  //sprintf(glb->outname, "%c", "");
  //sprintf(glb->dfname, "%c", "");
  
  strcat(glb->subparname, "../data/subpars/subpar.");
  strcat(glb->subparname, glb->RUN_ID);
  strcat(glb->subparname, ".dat");

  strcat(glb->parname, "../data/pars/par.");
  strcat(glb->parname, glb->RUN_ID);
  strcat(glb->parname, ".dat");

  strcat(glb->chainname, "../data/chains/chain.");
  strcat(glb->chainname, glb->RUN_ID);
  strcat(glb->chainname, ".dat");

  strcat(glb->logLname, "../data/logL/logL.");
  strcat(glb->logLname, glb->RUN_ID);
  strcat(glb->logLname, ".dat");

  strcat(glb->outname, "../data/lightcurves/mcmc_lightcurves/");
  strcat(glb->outname, glb->RUN_ID);
  strcat(glb->outname, ".out");

  strcat(glb->dfname, "../data/lightcurves/folded_lightcurves/");
  strcat(glb->dfname, glb->RUN_ID);
  strcat(glb->dfname, "_new.txt");

  char suffix[30];
  sprintf(suffix, "%d", run);
  strcat(glb->subparname, suffix);
  strcat(glb->parname, suffix);
  strcat(glb->chainname, suffix);
  strcat(glb->logLname, suffix);
  strcat(glb->outname, suffix);

  printf("\t Allocating memory...\n");
  glb->P_ = (double **)malloc(NCHAINS * sizeof(double));
  for (int i = 0; i < NCHAINS; i++)
    glb->P_[i] = (double *)malloc(NPARS * sizeof(double));

  glb->x = (double **)malloc(NCHAINS * sizeof(double));
  for (int i = 0; i < NCHAINS; i++)
    glb->x[i] = (double *)malloc(NPARS * sizeof(double));

  glb->history = (double ***)malloc(NCHAINS * sizeof(double));
  for (int i = 0; i < NCHAINS; i++)
  {
    glb->history[i] = (double **)malloc(NPAST * sizeof(double));
    for (int j = 0; j < NPAST; j++)
      glb->history[i][j] = (double *)malloc(NPARS * sizeof(double));
  }

  glb->sigma = (double *)malloc(NPARS * sizeof(double));
  glb->scale = 1.0 / ((double)(NPARS));

  printf("\t Setting parameter limits...\n");
  srand(glb->Niter);
  glb->seed = glb->Niter;
  set_limits(glb->limited, glb->limits);
  initialize_proposals(glb);

  // Set initial patameters
  printf("\t Assigning parameters...\n");
  int param_file_flag = file_exists(glb->parname);

  double tmp;

  if (param_file_flag)
  {
    printf("\t \t Opening parameter file... %s \n", glb->parname);
    glb->param_file = fopen(glb->parname, "r");
  }

  for (int i = 0; i < NPARS; i++)
  {
    if (param_file_flag)
    {
      fscanf(glb->param_file, "%lf\n", &tmp);
      printf("Par %d value: %.5e \n", i, tmp);
    }
    else
    {
      tmp = glb->limits[i].lo + ran2(&glb->seed) * (glb->limits[i].hi - glb->limits[i].lo);
      printf("\t \t Parameter file not found; Assinging random pars... \n");
    }

    glb->P_[0][i] = tmp;
    glb->x[0][i] = glb->P_[0][i];
    glb->xmap[i] = glb->P_[0][i];

    for (int j = 1; j < NCHAINS; j++)
    {
      glb->P_[j][i] = glb->P_[0][i];
      glb->x[j][i] = glb->P_[0][i];
    }

    for (int j = 0; j < NCHAINS; j++)
    {
      if (!param_file_flag)
      {
        if (glb->burn_in == 1)
        {
          glb->P_[j][i] = glb->limits[i].lo + ran2(&glb->seed) * (glb->limits[i].hi -
                                                                  glb->limits[i].lo);
          glb->x[j][i] = glb->P_[j][i];
        }

        if (glb->burn_in == 2)
        {
          glb->P_[j][i] = glb->limits[i].lo + ran2(&glb->seed) * (glb->limits[i].hi -
                                                                  glb->limits[i].lo);
          if (i == 2)
          {
            glb->P_[j][i] = glb->guess_period;
          } // Keep the period
          glb->x[j][i] = glb->P_[j][i];
        }
      }
    }
  }

  if (param_file_flag)
  {
    fclose(glb->param_file);
  }

  // Set an initial temperature
  glb->dtemp = 1.2;

  // Read data from the lightcurve file
  printf("\t Opening lightcurve file %s...\n", glb->dfname);
  int data_file_flag = file_exists(glb->dfname);

  if (!data_file_flag)
  {
    printf("\t Lightcurve file not found ...\n");
    exit;
  }

  if (glb->burn_in == 1)
  {
    double tmp1, tmp2, tmp3, tmp4;
    glb->data_file = fopen(glb->dfname, "r");
    tmp = 0;
    fscanf(glb->data_file, "%ld\n", &glb->Nt);
    glb->subN = 0;
    glb->rdata = (double *)malloc(4 * glb->Nt * sizeof(double));
    glb->index = (int *)malloc(NCHAINS * sizeof(int));

    for (int i = 0; i < glb->Nt; i++)
    {
      fscanf(glb->data_file, "%lf %lf %lf %lf\n", &tmp1, &tmp2, &tmp3, &tmp4);
      glb->rdata[i * 4] = tmp1;
      glb->rdata[i * 4 + 1] = tmp2;
      glb->rdata[i * 4 + 2] = tmp3;
      glb->rdata[i * 4 + 3] = tmp4;
      if (tmp4 == 0)
        glb->subN++;
    }
    printf("\t Closing lightcurve file... \n");
    fclose(glb->data_file);

    glb->a_data = (double *)malloc(glb->subN * sizeof(double));
    glb->a_model = (double *)malloc(glb->subN * sizeof(double));
    glb->t_data = (double *)malloc(glb->subN * sizeof(double));
    glb->e_data = (double *)malloc(glb->subN * sizeof(double));
    glb->q_data = (double *)malloc(glb->subN * sizeof(double));
    glb->light_curve = (double *)malloc(glb->subN * sizeof(double));
    glb->light_curve_2 = (double *)malloc(glb->subN * sizeof(double));

    glb->subi = 0;
    for (int i = 0; i < glb->Nt; i++)
    {
      tmp4 = glb->rdata[i * 4 + 3];
      if (tmp4 == 0)
      {
        glb->t_data[glb->subi] = glb->rdata[i * 4];
        glb->a_data[glb->subi] = glb->rdata[i * 4 + 1];
        glb->e_data[glb->subi] = sqrt(glb->rdata[i * 4 + 1]) * glb->true_err;
        glb->subi++;
      }
    }

    glb->dtemp = 1.3;
  }

  if (glb->burn_in == 2)
  {
    double tmp1, tmp2, tmp3;
    glb->data_file = fopen(glb->dfname, "r");
    tmp = 0;
    fscanf(glb->data_file, "%ld\n", &glb->Nt);
    glb->subN = 0;
    glb->rdata = (double *)malloc(3 * glb->Nt * sizeof(double));
    glb->index = (int *)malloc(NCHAINS * sizeof(int));

    for (int i = 0; i < glb->Nt; i++)
    {
      fscanf(glb->data_file, "%lf\t%lf\t%lf\n", &tmp1, &tmp2, &tmp3);
      glb->rdata[i * 3] = tmp1;
      glb->rdata[i * 3 + 1] = tmp2;
      glb->rdata[i * 3 + 2] = tmp3;
      glb->subN++;
    }
    printf("\t Closing data file... \n");
    fclose(glb->data_file);

    glb->a_data = (double *)malloc(glb->subN * sizeof(double));
    glb->a_model = (double *)malloc(glb->subN * sizeof(double));
    glb->t_data = (double *)malloc(glb->subN * sizeof(double));
    glb->e_data = (double *)malloc(glb->subN * sizeof(double));
    glb->q_data = (double *)malloc(glb->subN * sizeof(double));
    glb->light_curve = (double *)malloc(glb->subN * sizeof(double));
    glb->light_curve_2 = (double *)malloc(glb->subN * sizeof(double));

    glb->subi = 0;
    for (int i = 0; i < glb->Nt; i++)
    {
      glb->t_data[i] = glb->rdata[i * 3];
      glb->a_data[i] = glb->rdata[i * 3 + 1];
      glb->e_data[i] = glb->rdata[i * 3 + 2];
    }

    glb->dtemp = 1.3;
  }

  // Initialize parallel tempering
  glb->temp[0] = 1.0;
  glb->index[0] = 0;
  for (int i = 1; i < NCHAINS; i++)
  {
    glb->temp[i] = glb->temp[i - 1] * glb->dtemp;
    glb->index[i] = i;
  }

  // initialize likelihood
  glb->logLmap = loglikelihood(glb->t_data, glb->a_data, glb->e_data, glb->subN, glb->P_[0]);
  printf("\t Initial chi2 %g ... \n", -2 * glb->logLmap);

  printf("\t Opening chain and log files \n");
  glb->chain_file = fopen(glb->chainname, "w");
  glb->logL_file = fopen(glb->logLname, "w");

  /*
  calc_light_curve(glb->t_data, glb->subN, glb->xmap, glb->a_model);

  // Write the lightcurve to a file
  for (int k=0; k<NPARS; k++)
  {
    printf("Parameter %d value is %f \n", k, glb->P_[0][k]);
  }

  FILE *check_lc = fopen("generated_lc.txt", "w");

  for (int k=0; k<glb->subN; k++)
  {
      fprintf(check_lc, "%f \t %f \t %f  \t %f \t %f\n", glb->t_data[k], glb->a_data[k], 
      glb->a_model[k], glb->e_data[k], SQR((glb->a_model[k] - glb->a_data[k]) / glb->e_data[k]));
  }

  fclose(check_lc);
  */
}

void run_chain(int chain_id, int iter, struct globals *glb)
{

  double *y = (double *)malloc(NPARS * sizeof(double));
  double alpha = ran2(&glb->seed);
  double jscale = pow(10., -6. + 6. * alpha);
  long seed = glb->seed;
  int jump = 0;
  int DEtrial = 0;
  int DEacc = 0;
  int acc = 0;
  /* propose new solution */
  /* DE proposal; happens after 500 cycles */
  if (ran2(&seed) < 0.5 && iter > NPAST)
  {
    jump = 1;
  }
  // gaussian jumps along parameter directions
  if (jump == 0)
  {
    gaussian_proposal(glb->x[glb->index[chain_id]], &seed, glb->sigma,
                      jscale, glb->temp[chain_id], y);
  }
  // jump along correlations derived from chain history
  if (jump == 1)
  {
    if (glb->index[chain_id] == 0)
      DEtrial++;
    differential_evolution_proposal(glb->x[glb->index[chain_id]], &seed, glb->history[chain_id], y);
    double dx_mag = 0;
    for (int i = 0; i < NPARS; i++)
    {
      dx_mag += (glb->x[glb->index[chain_id]][i] - y[i]) * (glb->x[glb->index[chain_id]][i] - y[i]);
    }
    if (dx_mag < 1e-6)
      gaussian_proposal(glb->x[glb->index[chain_id]], &seed, glb->sigma, jscale, glb->temp[chain_id], y);
  }

  for (int i = 0; i < NPARS; i++)
  {
    // enforce priors (reflecting boundary conditions)
    if ((glb->limited[i].lo == 1) && (y[i] < glb->limits[i].lo))
      y[i] = 2.0 * glb->limits[i].lo - y[i];
    if ((glb->limited[i].hi == 1) && (y[i] > glb->limits[i].hi))
      y[i] = 2.0 * glb->limits[i].hi - y[i];
    // enforce priors (periodic boundary conditions)
    if ((glb->limited[i].lo == 2) && (y[i] < glb->limits[i].lo))
      y[i] = glb->limits[i].hi + (y[i] - glb->limits[i].lo);
    if ((glb->limited[i].hi == 2) && (y[i] > glb->limits[i].hi))
      y[i] = glb->limits[i].lo + (y[i] - glb->limits[i].hi);
  }

  // compute trial signal
  glb->logLx[glb->index[chain_id]] = loglikelihood(glb->t_data,
                                                   glb->a_data, glb->e_data, glb->subN, glb->x[glb->index[chain_id]]);
  double logLy = loglikelihood(glb->t_data, glb->a_data, glb->e_data, glb->subN, y);

  /* evaluate new solution */
  // Hasting's ratio
  double H = exp((logLy - glb->logLx[glb->index[chain_id]]) / glb->temp[chain_id]);
  // acceptance probability
  alpha = ran2(&seed);
  // printf("Log(Ly), Log(Lx): %12.5e %12.5e\n", logLy,logLx[index[j]]);
  // conditional acceptance of y


  if (alpha <= H)
  {
    //printf("Alpha is less than the heat \n");
    if (glb->index[chain_id] == 0)
      acc++;
    for (int i = 0; i < NPARS; i++)
    {
      glb->x[glb->index[chain_id]][i] = y[i];
    }
    glb->logLx[glb->index[chain_id]] = logLy;
    if (jump == 1 && glb->index[chain_id] == 0)
      DEacc++;
  }

  /* parallel tempering */
  ptmcmc(glb->index, glb->temp, glb->logLx);

  /*  fill history  */
  int k = iter - (iter / NPAST) * NPAST;
  for (int i = 0; i < NPARS; i++)
    glb->history[chain_id][k][i] = glb->x[glb->index[chain_id]][i];
  free(y);
  //#pragma omp barrier
}

void write_and_log(int iter, struct globals *glb)
{
  if (glb->logLx[glb->index[0]] > glb->logLmap)
  {
    for (int i = 0; i < NPARS; i++)
    {
      glb->xmap[i] = glb->x[glb->index[0]][i];
    }
    glb->logLmap = glb->logLx[glb->index[0]];
  }

  if (iter % 10 == 0)
  {
    // print parameter chains
    fprintf(glb->chain_file, "%ld %.12g ", iter / 10, glb->logLx[glb->index[0]]);
    for (int i = 0; i < NPARS; i++)
      fprintf(glb->chain_file, "%.12g ", glb->x[glb->index[0]][i]);
    fprintf(glb->chain_file, "\n");

    // print log likelihood chains
    fprintf(glb->logL_file, "%ld ", iter / 10);
    for (int i = 0; i < NCHAINS; i++)
      fprintf(glb->logL_file, "%.12g ", glb->logLx[glb->index[i]]);
    fprintf(glb->logL_file, "\n");
  }

  if (iter % 100 == 0)
  {
  	printf("%ld/%ld logL=%.10g \n",iter,glb->Niter,glb->logLx[glb->index[0]]);
    glb->data_file = fopen(glb->outname, "wb");


    if(!file_exists(glb->outname))
    {
      printf("Outfile %s does not exist \n", glb->outname);
      
    }

    if (glb->data_file == NULL)
    {
      printf("Cannot open the data file \n");
    }


    calc_light_curve(glb->t_data, glb->subN, glb->xmap, glb->a_model);
    fprintf(glb->data_file, "%d \n", glb->subN);

    for (int i = 0; i < glb->subN; i++)
    {
      fprintf(glb->data_file, "%12.5e %12.5e %12.5e\n", glb->t_data[i],
              glb->a_data[i], glb->a_model[i]);
    }
    fclose(glb->data_file);

    glb->param_file = fopen(glb->subparname, "w");
    fprintf(glb->param_file, "%12.5e %12.5e %12.5e %12.5e %12.5e\n",
            glb->x[glb->index[0]][0], glb->x[glb->index[0]][1], glb->x[glb->index[0]][2],
            glb->x[glb->index[0]][3], glb->x[glb->index[0]][4]);
    fprintf(glb->param_file, "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
            glb->x[glb->index[0]][5], glb->x[glb->index[0]][6], glb->x[glb->index[0]][7],
            glb->x[glb->index[0]][8], glb->x[glb->index[0]][9], glb->x[glb->index[0]][10]);
    if (ALPHA_FREE == 1)
    {
      fprintf(glb->param_file, "%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n",
              glb->x[glb->index[0]][11], glb->x[glb->index[0]][12], glb->x[glb->index[0]][13],
              glb->x[glb->index[0]][14], glb->x[glb->index[0]][15], glb->x[glb->index[0]][16],
              glb->x[glb->index[0]][17]);
    }
    fclose(glb->param_file);
  }
}

void initialize_proposals(struct globals *glb)
{
  int n, i, j;
  double junk;
  double x[NPARS];
  FILE *hfile;
  // jump sizes are set by hand based on what works
  glb->sigma[0] = 1.0e-1;      // log M1 (MSUN)
  glb->sigma[1] = 1.0e-1;      // log M2 (MSUN)
  glb->sigma[2] = glb->SIGMAP; // log P (days)
  glb->sigma[3] = 1.0e-1;      // e
  glb->sigma[4] = 1.0e+0;      // inc (deg)
  glb->sigma[5] = 1.0e+1;      // Omega (deg)
  glb->sigma[6] = 1.0e+1;      // omega0 (deg)
  glb->sigma[7] = 1.0e+0;      // T0 (day)
  glb->sigma[8] = 1.0e-2;      // log rr1 normalization
  glb->sigma[9] = 1.0e-2;      // log rr2 normalization
  if (ALPHA_FREE == 1)
  {
    glb->sigma[10] = 1.0e-2;     // limb darkening 1
    glb->sigma[11] = 1.0e-2;     // gravity darkening 1
    glb->sigma[12] = 1.0e-2;     // limb darkeing 2
    glb->sigma[13] = 1.0e-2;     // gravity darkening 2
    glb->sigma[14] = 1.0e-2;     // reflection 1
    glb->sigma[15] = 1.0e-2;     // reflection 2
    glb->sigma[16] = 1.0e-2;     // extra beaming 1
    glb->sigma[17] = 1.0e-2;     // extra beaming 2
    glb->sigma[18] = 1.0e-2;     // Teff 1
    glb->sigma[19] = 1.0e-2;     // Teff 2
    glb->sigma[20] = 1.0e-2;     // Blending
  }
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
  b = (int)((double)rand() / (RAND_MAX) * ((double)(NCHAINS - 1)));
  a = b + 1;

  olda = index[a];
  oldb = index[b];
  heat1 = heat[a];
  heat2 = heat[b];
  logL1 = logL[olda];
  logL2 = logL[oldb];
  dlogL = logL2 - logL1;
  H = (heat2 - heat1) / (heat2 * heat1);
  alpha = exp(dlogL * H);
  beta = (double)rand() / (RAND_MAX);
  if (alpha >= beta)
  {
    index[a] = oldb;
    index[b] = olda;
  }
}

void uniform_proposal(double *x, long *seed, bounds limits[], double *y)
{
  int n;
  double dx[NPARS];

  // compute size of jumps
  for (n = 0; n < NPARS; n++)
    dx[n] = ran2(seed) * (limits[n].hi - limits[n].lo);

  // uniform draw on prior range
  for (n = 0; n < NPARS; n++)
    y[n] = limits[n].lo + dx[n];
}

void gaussian_proposal(double *x, long *seed, double *sigma, double scale, double temp, double *y)
{
  int n;
  double gamma;
  double sqtemp;
  double dx[NPARS];

  // scale jumps by temperature
  sqtemp = sqrt(temp);

  // compute size of jumps
  for (n = 0; n < NPARS; n++)
    dx[n] = gasdev2(seed) * sigma[n] * sqtemp * scale;

    //printf("gasdev(seed) = %f \t sigma = %f \t sqtemp = %f \t scale = %f \n",
    //        gasdev2(seed), sigma[n], sqtemp, scale);

  // jump in parameter directions scaled by dx
  for (n = 0; n < NPARS; n++)
  {
    y[n] = x[n] + dx[n];
    //printf("dx_%d is %f \n", n, dx[n]);
    // printf("%12.5e %12.5e %12.5e %12.5e %12.5e\n",
    //	 scale,sqtemp,sigma[n],x[n],dx[n]);
  }
}

void differential_evolution_proposal(double *x, long *seed, double **history, double *y)
{
  int n;
  int a;
  int b;
  double dx[NPARS];

  // choose two samples from chain history
  a = ran2(seed) * NPAST;
  b = a;
  while (b == a)
    b = ran2(seed) * NPAST;

  // compute vector connecting two samples
  for (n = 0; n < NPARS; n++)
  {
    dx[n] = history[b][n] - history[a][n];
  }
  // Blocks?

  // 90% of jumps use Gaussian distribution for jump size
  if (ran2(seed) < 0.9)
    for (n = 0; n < NPARS; n++)
      dx[n] *= gasdev2(seed) * GAMMA;

  // jump along vector w/ appropriate scaling
  for (n = 0; n < NPARS; n++)
    y[n] = x[n] + dx[n];
}