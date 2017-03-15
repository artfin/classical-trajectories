#include <stdlib.h>
#include <stdio.h>

void print_int_array(int *b, int b_size) {
   for (int i = 0; i < b_size; i++) {
      printf("b[%d]: %d\n", i, b[i]);
   }
}

int* create_int_array(int length) {
    int *a = malloc(length * sizeof(int));
    if (a == NULL) {
        return a;
    }
    for (int i = 0; i < length; i++) {
        a[i] = i * 1000;
    }
    return a;
}

void double_int_array(int *a, int length) {
    for (int i = 0; i < length; i++) {
        a[i] *= 2;
    }
}

int main() {
    // creating an array of a given size and returning the pointer
    int a_size = 10;
    int *a = create_int_array(a_size);    

    print_int_array(a, a_size);
    printf("Doubling an array\n");
    double_int_array(a, a_size);
    print_int_array(a, a_size);

    free(a);
    return 0;
}
