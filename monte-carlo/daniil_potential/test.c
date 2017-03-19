#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "legendre_polynomial.h"

double potential(double *vars) {

    double V0, V2, V4, V6, V8, V10;

    double R = vars[0];
    double theta = vars[1];

    V0 = 22.4247 * exp(-0.716288 * R - 0.0869136 * R * R);

    if ( R >= 6.32925 ) {
       V0 -= -114.5 / pow(R, 6) + 2380. / pow(R, 8);
    } else {
        V0 += -0.599056 * exp(-0.650479 * R + 0.0320299 * R * R);
    }

    V2 = 63.5744 * exp(-0.811806 * R - 0.075313 * R * R);

    if ( R >= 6.90026 ) {
        V2 -= 26.6 / pow(R, 6) + 2080. / pow(R, 8);
    } else {
        V2 -= 0.207391 * exp(-0.620877 * R - 0.0310717 * R * R);
    }

    V4 = 99.1128 * exp(-1.17577 * R -0.0477138 * R * R);
    
    if ( R >= 6.96450) {
        V4  -= 410. / pow(R, 8);
    } else {
        V4 -= 0.0523497 * exp(-0.73534 * R - 0.0296750 * R * R);
    }

    V6 = 318.652 * exp(-1.88135 * R) - 0.0374994 * exp(-1.05547 * R - 0.0182219 * R * R);
    V8 = 332.826 * exp(-2.14596 * R) - 0.0137376 * exp(-1.03942 * R - 0.0494781 * R * R);
    V10 = 435.837 * exp(-2.44616 * R) - 0.108283 * exp(-2.04765 * R); 
   
    double cosine[1] = {cos(theta)};
   
    printf("given R: %lf\n", R);
    printf("given cosine value: %lf\n", cosine[0]);

    double *vals;
    int n;

    double legendre_values[6]; 
     
    for (int i = 0; i < 11; i++) {
        if ( i % 2 != 0 ) {
           continue;
        }
      
        vals = pm_polynomial_value (1, i, 0, cosine );

        legendre_values[i / 2] = vals[i];  
        free(vals);
    }

    for (int i = 0; i < 6; i++) {
        printf("legendre_values[%d]: %lf\n", i, legendre_values[i]);
    }
    
    double potential_value = V0 * legendre_values[0] + V2 * legendre_values[1] + V4 * legendre_values[2] + V6 * legendre_values[3] + V8 * legendre_values[4] + V10 * legendre_values[5];                                        
    
    printf("V0: %lf \n", V0);
    printf("V2: %lf \n", V2);
    printf("V4: %lf \n", V4);
    printf("V6: %lf \n", V6);
    printf("V8: %lf \n", V8);
    printf("V10: %lf \n", V10);
    return potential_value;
}

int main() {
    double vars[2] = {1., M_PI / 2};

    double val = potential(vars);
    printf("Potential value: %lf\n", val);
    return 0;
}
