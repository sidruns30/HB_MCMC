// test scripting file
#include <stdio.h>
#include <stdlib.h>

char *str1 = "something ";
char str2[] = "something else";

void somefunc(double xarr[3]){
    int len = sizeof(xarr);
    for (int i=0; i<len; i++){
        printf("%d \n", len);
        //printf("%f \n", xarr[i]);
    }
}

int main(){
    double temp = 10.;
    printf("temp is %f, %12.2e", temp, temp);
}