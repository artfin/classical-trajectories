#include <stdio.h>
#include <stdlib.h>

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} pixel;

void make_screen_redder(pixel *screen, int len) {
    for ( int i = 0; i < len; i++ ) {
        screen[i].r = 255;
    }
}

void print_screen(pixel *screen, int len) {
    for ( int i = 0; i < len; i++ ) {
        printf("pixel [%d]: %u %u %u\n", i, screen[i].r, screen[i].g, screen[i].b);
    }
} 

int main() {
    pixel *screen = malloc(10 * sizeof(pixel));
       
    for ( int i = 0; i < 10; i++ ) { 
        screen[i].r = 0;
        screen[i].g = 0;
        screen[i].b = 0;
    }
    
    printf("Just initialized screen:\n");
    print_screen(screen, 10);
    
    make_screen_redder(screen, 10);
    
    printf("Reddening screen:\n");
    print_screen(screen, 10);

    return 0;
}
