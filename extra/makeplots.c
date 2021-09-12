/*  Plotting script for C
    Creates a temporary file 'gnudata.temp' with x-y (or x-y-z) value pairs
    Then creates a pipe 'gnuplotPipe' to send commands to the gnu interface
    Commands (see below) can be passed to the plotting function in an array 
    Borrowed the core ideas from:   https://stackoverflow.com/questions/3521209/making-c-code-plot-a-graph-automatically
                             and    http://www.gnuplot.info/docs_4.0/gpcard.pdf
    Make sure to have gnuplot installed and added to PATH

    Useful gnuplot commands:
    plot a data file (2d):      plot 'data'
    plot a data file (3d):      splot 'data'
    plot with errorbars  :      plot errorbars 'data'                
    plot specific columns:      plot 'data' {using <xcol>:<ycol>, <ycol>:<xdelta><ydelta> (if using errorplot), etc.}
    set x, y ranges      :      plot 'data' {[xmin:xmax], [ymin:ymax]}
    use different styles :      plot 'data' with {<style>}| choose from lines, points, linespoints, impulses, dots, steps, errorbars
    set title            :      plot 'data' title "<title>"

    FUTURE:     Allow plotting of functions
                Delete data files after program terminates (use tmpfile?)
                No need to input npoints and nargs (use a wrapper function/struct to calculate the args)
                If a data file is already being created then maybe launch a .py file with matplotlib
                Make cpp version I hate C
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

// Plot 2d data given x and y arrays
void plot_2d(double *xarr, double *yarr, char args[], const int npoints){
    // First write the data to a file
    FILE *file = fopen("gnudata.temp", "w");
    for (int i=0; i < npoints; i++){fprintf(file, "%lf \t %lf \n", xarr[i], yarr[i]);}
    // Create a pipe to gnuplot interface ("-persistent" keeps the program open)
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    // Now send commands to the interface
    fprintf(gnuplotPipe, "plot \"gnudata.temp\"");
    fprintf(gnuplotPipe, "\t %s \n", args);
    fclose(gnuplotPipe);
    fclose(file);
}

void plot_2d_file(char *fname, char args[], const int npoints){
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    // Plot file and send commands 
    fprintf(gnuplotPipe, "plot '%s'", fname);
    fprintf(gnuplotPipe, "%s \n", args);
    fclose(gnuplotPipe);
}

int main(){
    double x[] = {-1.,2.,3.,4., 5};
    double y[] = {-2., 4., 2., 5., 10};
    const int NPOINTS = 5;
    char args[] = "title \"Simple plot \" with lines";
    plot_2d_file("gnudata.temp", args, NPOINTS);
}


/*Future*/
/*
// Wrapper function for plotting - take in similar args as matplotlib; change them for 
void plot(){}

// Plot from just file names
void plot_3d_file(char *fname, char *args[], const int npoints, const int nargs){
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "splot '%s'", fname);
    for (int i=0; i < nargs; i++){fprintf(gnuplotPipe, "%s \n", args[i]);}
}

void plot_2d_file(char *fname, char *args[], const int npoints, const int nargs){
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "plot '%s'", fname);
    for (int i=0; i < nargs; i++){fprintf(gnuplotPipe, "%s \n", args[i]);}
}

/*Testing purposes only, bad code*/
/*********************/
// Test functions for plotting
/*
double f2(double x){return sin(x);}
double f3(double x, double y){return sin(x) + sin(y);}
// Make uniform grid for plotting
void makegrid(const int size, double xmin, double xmax, double ymin, double ymax, double *x, double *y){
    double dx = (xmax-xmin)/size, dy = (ymax-ymin)/size;
    for(int i=0;i<size;i++){x[i] = xmin + i*dx;    y[i] = ymin + i*dy;}
}
// Make a datafile with function values
void gen_datfile(){
    // Generate grid
    const int NPOINTS=50, NARGS=1;
    double xarr[NPOINTS];
    double yarr[NPOINTS];
    double xmin=-1., xmax=1., ymin=-1., ymax=1;
    makegrid(NPOINTS, xmin, xmax, ymin, ymax, xarr, yarr);
    FILE *file = fopen("gnudata.temp", "w");
    // Write to a file
    for(int i=0;i<NPOINTS;i++){
        for(int j=0;j<NPOINTS;j++){
            fprintf(file, "%lf \t %lf \t %lf \n", xarr[i], yarr[j], f3(xarr[i], yarr[j]));
        }
    }
    fclose(file);
}

//////////////////////////////////////
void plot_3d(double *xarr, double *yarr, double *zarr, char *args[], const int npoints, const int nargs);
void plot_3d_file(char *fname, char *args[], const int npoints, const int nargs);
void plot_2d_file(char *fname, char *args[], const int npoints, const int nargs);

// Structure to contain all the basic gnuplot commands
typedef struct{
    // Labels
    char title[100];
    char xlabel[100];
    char ylabel[100];
    // Limits
    double xmin;
    double xmax;
    double ymin;
    double ymax;
    // Line properties
    char linetype[100];
    char pointtype[100];
    char color[100];
    double linewidth;
}gnu_prop;

gnu_prop *plot_prop;

// A 'constructor' to initialize the struct values
void init_values(gnu_prop *plot_prop){
    strcpy(plot_prop->title, " ");
    printf("The title is: %s", (*plot_prop).title);
}



// No need to pass the plot/splot command in the args
void plot_3d(double *xarr, double *yarr, double *zarr, char *args[], const int npoints, const int nargs){
    // First write the data to a file
    FILE *file = fopen("gnudata.temp", "w");
    for (int i=0; i < npoints; i++){fprintf(file, "%lf \t %lf \t %lf \n", xarr[i], yarr[i], zarr[i]);}
    // Create a pipe to gnuplot interface ("-persistent" keeps the program open)
    FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
    // Now send commands to the interface; call splot for 3d
    fprintf(gnuplotPipe, "splot 'gnudata.temp'");
    for (int i=0; i < nargs; i++){fprintf(gnuplotPipe, "%s \n", args[i]);}
    fclose(file);
}

// reading data script
    // First read the data to a file
    FILE *file = fopen("gnudata.temp", "r");
    double **filedata;
    // Allocate memory
    filedata = (double **)malloc(nvars*sizeof(double*));
    for (int i=0;i<nvars;i++){filedata[i] = (double *)malloc(npoints*sizeof(double));}
    // Read the contents of the file
    for (int i=0;i<npoints;i++){
        for (int j=0;j<nvars;j++){
            fscanf(file, "%lf\t", filedata[j][i]);
        }
    }
    // Cleanup
    for (int i=0;i<nvars;i++){double *ptr = filedata[i]; free(ptr);}
    free(filedata);
    fclose(file);
/*********************/
