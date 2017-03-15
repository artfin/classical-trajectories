#include <stdlib.h>
#include <stdio.h>

void change_int_array(int *a, int a_size) {
    for (int i = 0; i < a_size; i++) {
        a[i] = i * 100;
    }
}

void print_int_array(int *a, int a_size) {
    for (int i = 0; i < a_size; i++) {
        printf("a[%d]: %d\n", i, a[i]);
    }
}

int main() {
    
    // declaring a static array
    int a_size = 10;
    int a[a_size];
    int *ap = (int*) &a; // specifically casting it to an integer pointer
    
    print_int_array(ap, a_size);
    printf("Changing an array!\n");
    change_int_array(ap, a_size);
    print_int_array(ap, a_size); 

    return 0;
}
