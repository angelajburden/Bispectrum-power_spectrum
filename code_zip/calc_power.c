//g++ calc_power.c power_functions.c util.c integration.c io_gadget_slim.c -L/gpfs/apps/hpc/Libs/FFTW/2.1.5/lib -I/gpfs/apps/hpc/Libs/FFTW/2.1.5/include -lrfftw -lfftw -lm -o power_sp
//./power_sp  /gpfs/scratch60/fas/padmanabhan/np274/grace0/small/output/small128/sims/delta0.0/r1 snapshot_000 outputfile
//module load Libs/FFTW/2.1.5
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

//const int NG = 512;

//struct basic_gal *gal; 

//const double L = 1500.0;//3500.0;
//const double L = 300.0;
//const int NGAL_MAX = 1024*1024*1024;//130*130*130 + 20;

//parameters for box
//const int NGMAX = NG*NG*NG;
//const double RG= 15.0;	      
//const int NGK = NG*NG*2*(NG/2+1);
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
  char snapshot_base[NFCHAR], inpath[NFCHAR], outfname[NFCHAR], power_file[NFCHAR];
  if (argc !=10) {
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
  }
  NGAL_MAX = NGAL_MAX*NGAL_MAX*NGAL_MAX;
  fprintf(stderr," for grid length %d, %d k bins will be computed, box length =%lf, no. gals =%d \n", NG, NB, L, NGAL_MAX);

  long int NGAL =0;
  long int *ipNGAL =&NGAL;
  //parameters for box                                                                                                      
  const int NGK = NG*NG*2*(NG/2+1);

  //read in binary data                                                                                       

  fprintf(stderr,"read in binary file\n");
  fprintf(stderr,"allocate mem for delta\n");

  double *delta = (double*)malloc(sizeof(double)*NGK);
  //  fprintf(stderr,"number of gals = %d\n", NGAL);
  for (int i=0; i<NGK; i++) delta[i]=0.0;
  fprintf(stderr,"delta initiallised\n");
  read_in_binary_gadget(delta, NG, L, inpath, snapshot_base, snapshot, outfname, ipNGAL);

  printf("number of gals = %d\n", NGAL);
  double rho = (double)NGAL/(double)(NG*NG*NG);
  printf("rho = %lf\n", rho);

  FILE *fp_del_out;
  if((fp_del_out = fopen("del_out.txt","w"))==0)
    printf("cannot open pow output file\n");
 
  for (int i=0; i<NG; i++)
    for (int j=0; j<NG; j++)
      for (int k=0; k<NG; k++) {
        int ind = k + 2*(NG/2 + 1)*(j +NG*i);
        delta[ind] = (delta[ind]-rho)/ rho;
	fprintf(fp_del_out,"%g\n", delta[ind]);
   }
  printf("del computed\n");
  fclose(fp_del_out); 
  //  fun_calc_power(delta, NG, NGAL, L, power_file); 
   //   free(delta);
}

/**********************/
 

      
 
