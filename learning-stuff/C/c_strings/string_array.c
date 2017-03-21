#include <stdio.h>
#include <string.h>

int main() {
    // this is a pointer to the string array
    char *string_array[5] = {
        "hello", 
        "my",
        "name",
        "is",
        "Artem\n"
    };

    for ( int i = 0; i < 5; i++ ) {
        printf("%d: \"%s\"\n", i, string_array[i]);
    }

    return 0;
}
