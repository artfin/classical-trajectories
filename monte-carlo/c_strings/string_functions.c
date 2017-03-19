#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// brief overview of functions in string.h
int main() {
    char a[] = "hello world";
    
    // creating a pointer to the first element of character array
    char *ap = (char*) &a[0]; 
    char *bp = (char*) &a;
    char *cp = (char*) &a[6];

    printf("*ap: %p\n", ap);
    printf("data at ap: \"%s\"\n\n", ap);

    printf("*bp: %p\n", bp);
    printf("data at bp: \"%s\"\n\n", bp);
    
    printf("*cp: %p\n", cp);
    printf("data at cp: \"%s\"\n\n", cp);

    // explicitly putting null-terminating character in the middle of the string
    a[5] = '\0';
    
    printf("*ap: %p\n", ap);
    printf("data at ap: \"%s\"\n\n", ap);

    printf("*bp: %p\n", bp);
    printf("data at bp: \"%s\"\n\n", bp);
    
    printf("*cp: %p\n", cp);
    printf("data at cp: \"%s\"\n\n", cp);
    
    return 0;
} 
