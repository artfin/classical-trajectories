#include <stdio.h>
#include <stdlib.h>

typedef struct {
    unsigned char r;
    unsigned char g;
    unsigned char b;
} pixel;

typedef struct {
    int width;
    int height;
    pixel *data;
} screen;

// typically we need a function that will construct
screen* screen_new ( int width, int height ) {

    screen *s = malloc(sizeof(screen));

    s->width = width;
    s->height = height;
    s->data = malloc(width*height*sizeof(pixel));
    
    for ( int i = 0; i < s->height; i++ ) {
        for ( int j = 0; j < s->width; j++ ) {
            int index = i*s->width + j;
            s->data[index].r = i;
            s->data[index].g = j;
            s->data[index].b = 0;
        }
    }
    return s;
}

void screen_free(screen *s) {
    free(s->data);
    free(s);
}

void screen_print(screen *s) {

    for ( int i = 0; i < s->height; i++ ) {
        for ( int j = 0; j < s->width; j++) {
            int index = i * s->width + j;
            printf("(%u, %u, %u)\n", s->data[index].r, s->data[index].g, s->data[index].b);
        }
    }
}

int main() {
    screen *s = screen_new(5, 5);
    screen_print(s);
    screen_free(s);

    return 0;
}


