#include "header.h"
#include<stdio.h>
#include<rfftw.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <fstream>
#include <iomanip>
#include<complex>
#include<omp.h>
#define pi2 (pi*pi)
#define eps 1.0e-3// numerical accuracy for integrations
#define NEVAL 10000

//struct basic_gal *gal; 

char *data_file = cvector(0,200);
char *power_file = cvector(0,200);

int main(int argc, char **argv) {
  fprintf(stderr,"starting program\n");
  time_t t_codestart,t_codeend;                                                                
  t_codestart=time(NULL);  
  int NFCHAR = 1000;
  int snapshot;
  int NB, NG;//NB is number of k bins for each dimension, NG is number of 1d grid points
  double L = 0.0;
  long int NGAL_MAX = 0;
  int input_flag; 
  char snapshot_base[NFCHAR], inpath[NFCHAR], outfname[NFCHAR], power_file[NFCHAR], delta_file[NFCHAR];
  if (argc !=12) {
    fprintf(stderr,"\n\n wrong argument(s).  Specify: `%s' \n\n",argv[0]);
    fprintf(stderr,"<inpath>      (input path without trailing slash)\n");
    fprintf(stderr,"<basename>    (basename of snapshot files)\n");
    fprintf(stderr,"<snapshot>    (number of snapshot)\n");
    fprintf(stderr,"<outpath>     (output path)\n");
    fprintf(stderr,"<power_file>  (power_spectrum_file)\n");
    fprintf(stderr,"<NB>  (no. k bins)\n");
    fprintf(stderr,"<NG>  (no. grid points (1D))\n");
    fprintf(stderr,"<L>  (box length (1D))\n");
    fprintf(stderr,"<NGAL_MAX>  (no.particle ^(1/3) (1D))\n");
    fprintf(stderr,"<input_flag>  (readin in binary=0, read in delta=1\n");
    fprintf(stderr,"<delta_file>  (delta_file\n");
    fprintf(stderr,"\n exiting ..\n");
    exit(EXIT_FAILURE);
  }
  else{
    printf("copying inputs\n");
    strcpy(inpath,argv[1]);
    strcpy(snapshot_base,argv[2]);
    snapshot=atoi(argv[3]);
    strcpy(outfname,argv[4]);
    strcpy(power_file, argv[5]);
    NB=atoi(argv[6]);
    NG=atoi(argv[7]);
    L=atoi(argv[8]);
    NGAL_MAX=atoi(argv[9]);
    input_flag=atoi(argv[10]);
    strcpy(delta_file,argv[11]);
  }
  NGAL_MAX = NGAL_MAX*NGAL_MAX*NGAL_MAX;
  fprintf(stderr," for grid length %d, %d k bins will be computed, box length =%lf, no. gals =%d \n", NG, NB, L, NGAL_MAX);

  long int NGAL =0;
  long int *ipNGAL =&NGAL;

  //parameters for box
  const int NGK = NG*NG*2*(NG/2+1);
  //read in binary data
  if (input_flag==0)fprintf(stderr,"read in binary file\n");
  if (input_flag==1)fprintf(stderr,"read in delta file\n");
  fprintf(stderr,"allocate mem for delta\n");
  double *delta = (double*)malloc(sizeof(double)*NGK);
  fprintf(stderr,"number of gals = %d\n", NGAL);
  for (int i=0; i<NGK; i++) delta[i]=0.0;
  fprintf(stderr,"delta initiallised\n");
  double rho=1;
  if (input_flag==0) {
    read_in_binary_gadget(delta, NG, L, inpath, snapshot_base, snapshot, outfname, ipNGAL);
    fprintf(stderr,"number of gals = %d\n", NGAL);
    rho = (double)NGAL/(double)(NG*NG*NG);
  }
  if (input_flag==1)rho = read_in_delta_field(delta, NG, delta_file);
  
  printf("rho = %lf\n", rho);


  for (int i=0; i<NG; i++)
    for (int j=0; j<NG; j++)
      for (int k=0; k<NG; k++) {
        int ind = k + 2*(NG/2 + 1)*(j +NG*i);
        delta[ind] = (delta[ind]-rho)/ rho;
  }
  fun_calc_bispectrum(delta, NG, NGAL, NB, L, power_file); 
  free(delta);
  fprintf(stderr," for a grid size %d with %d k bins the bispectrum code takes \n", NG, NB);
  t_codeend = time(NULL);
  print_time(t_codestart,t_codeend,"this long"); 
}

/**********************/
 

      
 
