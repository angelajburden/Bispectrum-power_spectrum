#include <stdlib.h>
#include <stdio.h>
#define NUM_POINTS 5
#define NUM_COMMANDS 6
//set term pngcairo;                                                                                               
//set output "filename.png";

void xy_plot(double *k, double *pk, int log_flag, int NB){
  char * commandsForGnuplot[] = {"set title \"power spectrum\"", "plot 'data.temp'", " set logscale xy", "set xrange [0.005:0.5]", "set xlabel {\"k [h/Mpc]\"}", "set ylabel {\"P(k)\"}"};

  FILE * temp = fopen("data.temp", "w");
  FILE * gnuplotPipe = popen ("gnuplot -persist", "w");
  int i;
  for (i=0; i < NB; i++) fprintf(temp, "%lf  %lf \n", k[i], pk[i]); //Write the data to a temporary file
  for (i=0; i < NUM_COMMANDS; i++) fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]); 
  //  set term pngcairo;
  //  set output "filename.png";
  //replot;
  //set output;
  //  return 0;
}
