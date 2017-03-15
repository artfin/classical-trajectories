//
//  main.cpp
//  monte-carlo
//
//  Created by Mac on 14.03.17.
//  Copyright (c) 2017 Mac. All rights reserved.
//

#include <iostream>
#include <stdlib.h>

using namespace std;

double integrate(double a, double b, int N);
double f(double x);

int main() {
    cout << integrate(0, 1, 100000000) << endl;
}

double integrate(double a, double b, int iterations) {
    if (a > b) return integrate(b, a, iterations);
    
    double sum = 0;
    double x;
    
    for (int i = 0; i < iterations; i++) {
        x = (double) rand() / RAND_MAX;
        sum += f(a + (b - a) * x);
    }
    
    sum *= (b - a) / iterations;
    return sum;
}

double f(double x) {
    return x * x;
}
