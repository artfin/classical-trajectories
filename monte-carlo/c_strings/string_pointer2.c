#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main() {
    // allocating space in the static memory
    char *a = malloc(6 * sizeof(char));
    printf("%p\n", a); // printing the address of allocated memory
    strncpy(a, "hello", 5); // copying the string to the allocated memory (5 characters at most)
    printf("%p\n", a); // making sure that the copy is placed exactly in the allocated space
    a[5] = '\0'; // explicitly adding null-terminating symbol to the end of the string

    puts(a);
    
    return 0;
}


