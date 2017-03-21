#include <stdio.h>

int main() {
    char a[] = "hello";
    char b[6] = "hello";
    char c[5] = "hello";
    
    char d[] = {'h', 'e', 'l', 'l', 'o', '\0'}; 
    char e[] = {'h', 'e', 'l', 'l', 'o'};

    // just a quicker way to print a string 
    // puts(a);
    
    printf("a: \"%s\"\n", a);
    printf("b: \"%s\"\n", b);
    printf("c: \"%s\"\n", c);
    printf("d: \"%s\"\n", d);
    printf("e: \"%s\"\n", e);
    
    return 0;
}
