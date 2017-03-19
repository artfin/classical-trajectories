#include <stdio.h>

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} pixel; // providing a structure name

void print_pixel(pixel p) {
    printf("%u %u %u\n", p.r, p.g, p.b);
}

int main() {

    pixel p;
    p.r = 0;
    p.g = 255;
    p.b = 0;

    printf("sizeof(pixel): %lu\n", sizeof(pixel));
    // sizeof(pixel) gives 3 bytes

    print_pixel(p);

    return 0;
}

