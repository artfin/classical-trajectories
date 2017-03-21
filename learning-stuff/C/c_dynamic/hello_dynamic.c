#include <stdlib.h>
// malloc function is in here

#include <stdio.h>

int main() {
    int a_size = 10;
    int *a = NULL;
    // declaring a pointer that points nowhere
    
    a = malloc(a_size * sizeof(int));
    
    // malloc returns NULL if it can't allocate memory
    // it's very rare but could happen
    if (a == NULL) {
        fprintf(stderr, "Error allocating memory.\n");
        exit(1);
    }

    for (int i = 0; i < a_size; i++) {
        a[i] = i * 10;
        printf("a[%d]: %d\n", i, a[i]);
    }
    
    // freeing memory
    free(a);

    return 0;
}
