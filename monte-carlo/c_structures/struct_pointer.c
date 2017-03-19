#include <stdio.h>

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} pixel;

void print_pixel(pixel *p) {
    printf("%u %u %u\n", p->r, p->g, p->b);
}

void make_redder(pixel *p) {
    p->r = 255;
}

int main() {
    pixel p = {0, 255, 0};
    
    print_pixel(&p);
    make_redder(&p);
    print_pixel(&p);

    return 0;
}
