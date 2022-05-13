compile with gcc:
gcc -O3 -o [bin_name] -lm -std=c99 -fopenmp mcmc_wrapper2.c likelihood3.c
compile with icc:
icc -O3 -o [bin_name] -std=c99 -qopenmp mcmc_wrapper2.c likelihood3.c

run:
./[bin_name] (int) N_iterations (int) TIC_ID (float) lc_period (int) run_id

comments:
- openmp flag is set by 'ENABLE_OPENMP' in mcmc_wrapper2.h
- lc data file must be located in 'data/lightcurves/folded_lightcurves/<TIC_ID>_new.txt'
- mag file must be located in 'data/magnitudes/<TIC_ID>.txt'
- the input format of lc file is expected to be:
line/col  1   2   3
1         Npts
2         t0  f0  e0
3         t1  f1  e1
        ...
Npts+1    tN  fN  eN

- the input format of the mag file is expected to be:
line/col  1   2   3
1        dist
2        Gmag Gmage
3        B-V  B-V_e
4        V-G  V-G_e
5        G-T  G-T_e
