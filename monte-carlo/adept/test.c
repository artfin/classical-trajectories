#include "adept.h"

using adept :: adouble;

adouble algorithm(const adouble x[2]) {
    adouble y = 4.0;
    adouble s = 2.0 * x[0] + 3.0 * x[1] * x[1];
    y *= sin(s);
    return y;
}

int main() {

    double *x[2] = {1.0, 2.0};
    double res = algorithm(x);

    return 0;
}
