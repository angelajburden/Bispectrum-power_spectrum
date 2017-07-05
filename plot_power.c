
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <fstream>
#include <iomanip>
#include<complex>
#include<omp.h>
#define NUM_COMMANDS 6
#define pi2 (pi*pi)
#define eps 1.0e-3// numerical accuracy for integrations                                                                     
#define NEVAL 10000
char *power_file = (char*)malloc( sizeof(char) * ( 100 ) );
int NBIN_MAX = 2100;
int NBIN =0;
struct power_spectrum {
  double k;
  double pk;
};
struct power_spectrum *pk1;

int main(int argc, char **argv) {

  if(!(pk1 = (struct power_spectrum*)malloc(NBIN_MAX*sizeof(struct power_spectrum))-1))
    printf("memory allocation problem for pk1\n");
  FILE *fp_open;
  if((fp_open=fopen("/home/fas/padmanabhan/ab2733/Cosmic_emu/CosmicEmu_v1.0/power_emu.txt","r"))==NULL) printf("emu file not opened\n");
  const int bsz=300; char buf[bsz];
  while((fgets(buf, bsz, fp_open))!=NULL) {
    double k, p1;
    sscanf(buf,"%lf %*lf %*lf %*lf %*lf %*lf %lf\n",&k, &p1);
    if(++(NBIN) > NBIN_MAX) { NBIN--; break; }
    pk1[NBIN].k = k;
    pk1[NBIN].pk = p1;
  }
  printf("p1.k=%lf, p1=%lf\n", pk1[10].k, pk1[10].pk);   
  char * commandsForGnuplot[] = {"set title \"power spectrum\"", "plot 'data.temp'", " set logscale xy", "set xrange [0.005:0.2]"};

    FILE * temp = fopen("data.temp", "w");
    FILE * gnuplotPipe = popen ("gnuplot -persist", "w");
    int i;
    for (i=1; i < NBIN; i++) fprintf(temp, "%lf  %lf \n", pk1[i].k, pk1[i].pk); //Write the data to a temporary file                           
    printf("data in temp file\n");
    for (i=0; i < 4; i++) fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
  
  

}
