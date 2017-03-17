#include <stdlib.h>
#include <stdio.h>

#define NUM_POINTS 5
#define NUM_COMMANDS 2

int main() {
    char *commandsForGnuplot[] = {"set title \"TITLEEEEE\"",
        "set xr[0.0:10.0]"};
    double xvals[NUM_POINTS] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double yvals[NUM_POINTS] = {5.0, 3.0, 1.0, 3.0, 5.0};

    FILE * temp = fopen("data.temp", "w");
    FILE *gnuplotPipe = popen( "gnuplot -persistent", "w" );

     /*
     Opens an interface that one can use to send commands as if they were typing into the gnuplot command line.  "The -persistent" keeps the plot open even after your C program terminates.
     */


    for (int i = 0; i < NUM_POINTS; i++) {
       fprintf(temp, "%lf %lf\n", xvals[i], yvals[i]);
    }

    printf("Command to gnuplot\n");
  
    for (int i = 0; i < NUM_COMMANDS; i++) {
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
    }

    fprintf(gnuplotPipe, "%s\n", "plot 'data.temp' with lp"); 
    fprintf(gnuplotPipe, "e");

    return 0;
} 
