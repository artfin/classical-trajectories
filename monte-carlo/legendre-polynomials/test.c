#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "legendre_polynomial.h"

double* create_double_array(int length) {
    double *a = malloc(length * sizeof(double));

    // checking that memory allocating was successful
    if (a == NULL) {
        fprintf(stderr, "Error allocating memory in 'create_double_array'!");
        exit(1);
    }

    // returning a pointer to the allocated array
    return a;
}

double* populate_array(double* a, int SIZE, double start, double end) {
  
    double step = (end - start) / (SIZE - 1);

    for (int i = 0; i < SIZE; i++) {
        a[i] = start + step * i;
    }

    // returns the pointer to the populated array
    return a;
}

void print_double_array(double* a, int length) {
    for (int i = 0; i < length; i++) {
        printf("a[%d]: %lf\n", i, a[i]);
    }
}

double* sample_legendre_polynomial(double* x, int SIZE, int n) {
  
   int m = 0; 
   double _x[1];
   double *vals;

   double *y = create_double_array(SIZE);

   for ( int i = 0; i < SIZE; i++ ) {
       _x[0] = x[i];
       vals = pm_polynomial_value ( 1, n, m, _x ); 
       y[i] = vals[n];
       free(vals);
   }

   return y;
}

void save_values(double* x, double* y, int SIZE, char* filename) {
    FILE *tfile = fopen(filename, "w");

    for ( int i = 0; i < SIZE; i++ ) {
        fprintf(tfile, "%lf %lf\n", x[i], y[i]);
    } 
}

int main() {
    
    double start = -1.0;
    double end   =  1.0;
    int SIZE = 100;
      
    // here we have a pointer pointing to a newly dynamically allocated  array
    double *x_vec = create_double_array(SIZE);
    double *y_vec;

    int i = 5;

    x_vec = populate_array(x_vec, SIZE, start, end);
    y_vec = sample_legendre_polynomial(x_vec, SIZE, i);
    
    /*for ( int i = 0; i < SIZE; i++ ) {*/
    /*printf("x: %lf; y: %lf\n", x_vec[i], y_vec[i]);*/
    /*}*/
    
    char filename[15];
    sprintf(filename, "legendre_%d.dat", i);
    printf("filename: %s\n", filename);

    save_values(x_vec, y_vec, SIZE, filename);

    return 0;
}
