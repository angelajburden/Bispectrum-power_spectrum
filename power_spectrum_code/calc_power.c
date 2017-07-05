//g++ calc_power.c power_functions.c util.c integration.c io_gadget_slim.c -lfftw3 -lm -o power_sp
//./power_sp  /gpfs/scratch60/fas/padmanabhan/np274/grace0/small/output/small128/sims/delta0.0/r1 snapshot_000 outputfile
//Libs/FFTW/3.3-openmpi-gcc
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

const int NG = 512;

struct basic_gal *gal; 

//const double L = 1150.0;//3500.0;
const double L = 1500.0;
const int NGAL_MAX = 1024*1024*1024;//130*130*130 + 20;

//parameters for box
const int NGMAX = NG*NG*NG;
const double RG= 15.0;	      
const int NGK = NG*NG*2*(NG/2+1);
char *data_file = cvector(0,200);
char *power_file = cvector(0,200);

int main(int argc, char **argv) {
  printf("starting program\n");
  int NFCHAR = 1000;
  int snapshot;
  char snapshot_base[NFCHAR], inpath[NFCHAR], outfname[NFCHAR], power_file[NFCHAR]; 
  if (argc !=6) {
    fprintf(stderr,"\n\n wrong argument(s).  Specify: `%s' \n\n",argv[0]);
    fprintf(stderr,"<inpath>      (input path without trailing slash)\n");
    fprintf(stderr,"<basename>    (basename of snapshot files)\n");
    fprintf(stderr,"<snapshot>    (number of snapshot)\n");
    fprintf(stderr,"<outpath>     (output path)\n");
    fprintf(stderr,"<power_file>  (power_spectrum_file)\n");
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
  }
  printf("allocating galaxy structure memory\n");
  if(!(gal = (struct basic_gal*)malloc(NGAL_MAX*sizeof(struct basic_gal))-1))
    printf("memory allocation problem for galaxies\n");
  printf("memory allocated\n");
  int NGAL =0;
  int *ipNGAL =&NGAL;

  //read in binary data
  printf("read in binary file\n");
  read_in_binary_gadget(gal, inpath, snapshot_base, snapshot, outfname, ipNGAL);
  double *delta = (double*)malloc(sizeof(double)*NGK);
  printf("number of gals = %d\n", NGAL);
  for (int i=0; i<NGK; i++) delta[i]=0.0;
  printf("delta initiallised\n"); 
  for (int i=0; i<10; i++) printf("x=%lf, y=%lf, z=%lf\n",gal[i].cp[0], gal[i].cp[1], gal[i].cp[2]);
  for (int i=NGAL-1; i>NGAL-10; i--) printf("x=%lf, y=%lf, z=%lf\n",gal[i].cp[0], gal[i].cp[1], gal[i].cp[2]);
  double rho = (double)NGAL/(double)(L*L*L);
  printf("rho = %lf\n", rho);
  for (int i=0; i<NGAL; i++) {
    int xp = (int)(((gal[i].cp[0])/(double)L) *(double) NG);
    int yp = (int)(((gal[i].cp[1])/(double)L) * (double)NG);
    int zp = (int)(((gal[i].cp[2])/(double)L) *(double)NG);
    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 
    if (i <10 || i >(NGAL-10)) printf("xp =%d, yp=%d, zp=%d\n", xp, yp, zp);
    if (xp>NG || xp<0) printf("gal xp[%d] = %d, galx = %lf\n", i, xp, gal[i].cp[0]);
    if (yp>NG || yp<0)printf("gal yp[%d] = %d, galy = %lf\n", i, yp, gal[i].cp[1]);
    if (zp>NG || zp<0) printf("gal zp[%d] = %d\n", i, zp);
    delta[ind] = delta[ind]+1.0; 
  }
  for (int i=0; i<NG; i++)
    for (int j=0; j<NG; j++)
      for (int k=0; k<NG; k++) {
        int ind = k + 2*(NG/2 + 1)*(j +NG*i);
        delta[ind] = (delta[ind]/rho) -1.0;
      }
   calc_power(delta, NG, NGAL, L, power_file); 
   free(delta);
}

/**********************/
 

      
 
