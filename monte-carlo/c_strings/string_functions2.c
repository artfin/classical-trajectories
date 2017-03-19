#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    char a[] = "Humpy dumpty sat on a wall";
    char search_char = 'u';
    char replace_char = 'o';
    
    char *ap = (char*) &a[0]; // placing a pointer at the very beginning of a character array, explicitly casting it to be a character pointer
    
    printf("Before: %s\n", a);

    // construction in the braces will return NULL as soon as it can't find any more occurences of the search_char in the given string, which we're going to truncate moving pointer along 
    while ( (ap = strchr(ap, search_char)) != NULL ) {
        printf("In the middle: %s\n", ap);
        ap[0] = replace_char;
    }

    printf("After: %s\n", a);

    return 0;
}
