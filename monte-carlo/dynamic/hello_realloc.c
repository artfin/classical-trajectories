#include <stdlib.h>
#include <stdio.h>

int main() {
    int a_size = 10;
    int new_a_size = 0;
    int *a = malloc(a_size * sizeof(int));

    if (a == NULL) {
        fprintf(stderr, "Error allocating memory!\n");
        exit(1);
    }

    for( int i = 0; i < a_size; i++) {
        a[i] = i * 10;
        printf("a[%d]: %d\n", i, a[i]);
    }
    
    printf("Reallocating a\n");

    // creating new pointer
    int *tmp = NULL;
    new_a_size = a_size * 2;

    tmp = realloc(a, new_a_size * sizeof(int));

    if (tmp == NULL) {
        fprintf(stderr, "Error reallocating a!\n");
        exit(1);
    }

    a = tmp;

    for ( int i = 0; i < new_a_size; i++) {
        printf("a[%d]: %d\n", i, a[i]);
    }

    printf("Done!\n");

    free(a);

    return 0;
}
