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

char *power_file = cvector(0,200);

int main(int argc, char **argv) {
  fprintf(stderr,"starting program\n");
  time_t t_codestart,t_codeend;
  t_codestart=time(NULL);


  //compute redshift space kernels
  double *Z1 = (double*)malloc(sizeof(double)*NKBIN);
  double *Z2 = (double*)malloc(sizeof(double)*(NKBIN*NKBIN));
  double k_mag_sq=0.0;
  double k_mag=0.0;
  double kx1, kx2, ky1, ky2, kz1, kz2;
  for (int i=0; i<Nk; i++)  {
    printf("i=%d\n",i);
    for (int j=0; j<Nk; j++)  {
      for (int l=0; l<Nk; l++) {
	int ink = l + Nk*(j + i*Nk);
       
        kx1 = (double) i*0.012 + 0.012;
        ky1 = (double) j*0.012 + 0.012;
        kz1 = (double) l*0.012 + 0.012;
        k1_mag_sq = (kx1*kx1 + ky1*ky1 + kz1*kz1);
        k1_mag = sqrt(k1_mag_sq);
        if (k1_mag_sq ==0) continue;
        k1_out_vals[ink] = k1_mag;
        double theta = acos(kz1/k1_mag);
        double phi = atan2(ky1,kx1);
        double mu1 = 1.0/(k1_mag*sqrt(3))*(kx1 + ky1 + kz1);
        for (int i2=0; i2<Nk; i2++) {
          for (int j2=0; j2<Nk; j2++) {
            for (int l2=0; l2<Nk; l2++) {
              int ink2 = l2 + Nk*(j2 + Nk*i2);
              kx2 = (double) i2*0.012 + 0.012;
	      ky2 = (double) j2*0.012 + 0.012;
	      kz2 = (double) l2*0.012 + 0.012;
              k2_mag_sq = (kx2*kx2 + ky2*ky2 + kz2*kz2);
              k2_mag = sqrt(k2_mag_sq);
	      double mu2 = 1.0/(k2_mag*sqrt(3))*(kx2 + ky2 + kz2);
              double alpha =(kx1*kx2 + ky1*ky2 + kz1*kz2)/(k2_mag*k1_mag);

	      //double damp = exp(-k_mag_sq*(1 + f*(f+2)*mu*mu)*A);
              //double da = exp(-k_mag_sq*A);
              double PLk=0.0;
	      //splint(mode, P_lin, y2, Nk_read, k_mag, &PLk);
              double Z1a = (b1 + f*mu1*mu1);
              double Z1b = (b1 + f*mu2*mu2);
              double F12 = 5./7. + (1./2.)*cos(alpha)*(k1_mag/k2_mag + k2_mag/k1_mag) + (2./7.)*cos(alpha)*cos(alpha);
              double G12 = 3./7. + (1./2.)*cos(alpha)*(k1_mag/k2_mag + k2_mag/k1_mag) + (4./7.)*cos(alpha)*cos(alpha);
              double mu = (mu1*k1_mag + mu2*k2_mag)/(kx1+kx2 ) 
              double Spi = sqrt(pi)/2;

        double av_lin_a = Spi*mu_0[ink]*da*PLk;
        double av_lin_b = -Spi*mu_1[ink]*da*PLk*2*f;
        double av_lin_c = Spi*mu_2[ink]*da*PLk*f*f;

        double factq = (k_mag*k_mag*k_mag/(4*pi*pi));

        //loop over all x values and r values                                                                
        for (int ir=0; ir<Nr; ir++) {
          double r = r_val[ir];

          for (int ix=0; ix<Nx; ix++) {
            double x = x_val[ix];
    
