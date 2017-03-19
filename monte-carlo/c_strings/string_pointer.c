#include <stdio.h>
int main() { 
    // read-only string
    char *a = "hello";
    char b[] = "hello";

    printf("a: \"%s\"\n", a);
    printf("b: \"%s\"\n", b);
    
    b[0] = 'b';
    printf("b: \"%s\"\n", b);
    // a[0] = 'c';
    // printf("a: \"%s\"\n", a);
    // will produce a Segmentation Fault

    return 0;
}
