// test scripting file
#include <stdio.h>
#include <stdlib.h>

#define arr_size(x) (sizeof(x)/sizeof(x[0]))

char *str1 = "something ";
char str2[] = "something else";

static inline void swap(double *x, double *y){
    double temp = *x;
    *x = *y;
    *y = temp;
}

void wrapper(double x, double y){
    swap(&x, &y);
}

int main(){
    double x = 4, y=2;
    wrapper(x, y);
    printf("x and y are: %lf, %lf", x, y);
    //somefunc(x);
}
/*
        FILE *file = fopen("gnudata.temp", "r");
    double x[3] = {1,2,3};
    int size = arr_size(x);
    printf("arr size is: %d \n", sizeof(double*));
    double (*test)[2] = malloc(sizeof(double[2][3]));
    for (int i=0;i<2;i++){
        for (int j=0;j<3;j++){
            test[i][j] = i;
            printf("%d\n", test[i][j]);
            }
    }
        for (int i=0;i<2;i++){
        for (int j=0;j<3;j++){
            printf("%d\n", test[i][j]);
            }
    }
    free(test);
    FILE *file = fopen("../gnudata.temp", "r");
    if (file == NULL){printf("FAIL\n");}
    const int nrows=5, ncols=2;
    double **data;
    // We want to index data by nrows x ncols
    data = (double **)malloc(nrows*sizeof(double*));
    for (int i=0;i<nrows;i++){data[i] = (double*)malloc(ncols*(sizeof(double)));}
    // Read the data file
    char temp[100];
    for (int i=0;i<nrows*ncols;i++){
        fscanf(file, "%s\n",&temp);
        double temp_d = atof(temp);
        printf("value of temp is: %lf", temp_d);
    }

    //free the data
    for (int i=0;i<nrows;i++){double *ptr = data[i]; free(ptr);}
    free(data);
*/