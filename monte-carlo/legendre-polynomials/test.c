#include <math.h>
#include <stdio.h>
#include <stdlib.h>
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

double* populate_array(double* a, int length, double start, double step) {
   // accepts a pointer, an array length, a start value and step 
   // and populates given array with values  

    for (int i = 0; i < length; i++) {
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

int main() {
    
    double start = 0.0;
    double end = 1.0;
    double step = 0.01;
    int x_size = (end - start) / step;
    
    // here we have a pointer pointing to a newly dynamically allocated  array
    double *x = create_double_array(x_size);

    x = populate_array(x, x_size, start, step);

    for (int i = 0; i < 10; i++) {
        // pm_polynomial():
        // int MM -- the number of evaluation points
        // int N -- the maximum first index of the Legendre function
        // int M -- the second index; 0 for simple Legednre functions
        // double x[MM] -- the points at which the function is to be evaluated

        value = pm_polynomial_value(1, i, 0, x);
        printf("N: %d \t\t X: %lf \t\t value: %lf", i, x[0], value);
    } 

    int n = 0;
    int m = 0;
    double x = 1.0;
    double x_vec[1];
    double *fx2_vec;

    // pm_polynomial_values ( &n_data, &n, &m, &x, &fx1 );
    
    x_vec[0] = x;
    fx2_vec = pm_polynomial_value (1, n, m, x_vec );

    printf("n: %d\n", n, m);
    printf("x: %lf\n", x);
    printf("fx2_vec: %lf\n", *fx2_vec); 

    return 0;
}
