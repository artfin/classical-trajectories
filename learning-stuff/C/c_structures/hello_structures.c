#include <stdio.h>

// declaring a structure in the global scope
struct pixel {
    unsigned char r;
    unsigned char g;
    unsigned char b;
};

void print_pixel(struct pixel p) {
    printf("%u %u %u\n", p.r, p.g, p.b);
}

int main() {

    // could do like that 
    // struct pixel p = {0, 255, 0};
    
    // or using dot-notation
    struct pixel p;
    p.r = 0;
    p.g = 255;
    p.b = 0;

    print_pixel(p);

    return 0;
}
