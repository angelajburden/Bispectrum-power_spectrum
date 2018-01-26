//Calculates the power spectrum in redshift space with LRT based on https://arxiv.org/pdf/0711.2521.pdf
//Resumming Cosmological Perturbations via the Lagrangian Picture

//NB this is the old program that atempts to calculate the integrals.
//They do not converge so use the compare_LRT.c program with the integrals pre-calculated in mathematica.
#include "header.h"
#include<stdio.h>
#include<srfftw.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include <fstream>
#include <iomanip>
#include<complex>
#include<omp.h>
#define pi2 (pi*pi)
#define eps 1.0e-3 // numerical accuracy for integrations
#define JMAX 20
#define NEVAL 10000
#define CHUNKSIZE 1; 

int Nk =10; //number of k bins
int Nr =100;
int Nx =100;
int Nk_read = 450;
int NKT = Nk*Nk*Nk;
//declare functions
double function_lin (double);
double function_q1 (double);
double function_q2 (double);
double function_q3 (double);
double function_q4 (double);
double function_qr1a (double);
double function_qr2a (double);
double function_qr3a (double);
double function_qr4a (double);
double function_qr1b (double);
double function_qr2b (double);
double function_qr3b (double);
double function_qr4b (double);
double function_qr1c (double);
double function_qr2c (double);
double function_qr3c (double);
double function_qr4c (double);
double function_qr1d (double);
double function_qr2d (double);
double function_qr3d (double);
double function_qr4d (double);
double function_qr1e (double);
double function_qr2e (double);
double function_qr3e (double);
double function_qr4e (double);
double function_r1 (double);
double function_r2 (double);
double mu1 (double, int);
double *P_red = (double*)malloc(sizeof(double)*NKT);
double *k_out_vals = (double*)malloc(sizeof(double)*NKT);
double *P_average = (double*)malloc(sizeof(double)*NKT);
double *mode = (double*)malloc(sizeof(double)*Nk_read);
double *r_val = (double*) malloc(sizeof(double)*Nr);
double *x_val = (double*) malloc(sizeof(double)*Nx);

double *P_lin = (double*)malloc(sizeof(double)*Nk_read);
double *y2 = (double*)malloc(sizeof(double)*Nk_read);

double *func_x1 = (double*)malloc(sizeof(double)*Nx);
double *func_x2 = (double*)malloc(sizeof(double)*Nx);
double *func_x3 = (double*)malloc(sizeof(double)*Nx);
double *func_x4 = (double*)malloc(sizeof(double)*Nx);
double *y2x1 = (double*)malloc(sizeof(double)*Nx);
double *y2x2 = (double*)malloc(sizeof(double)*Nx);
double *y2x3 = (double*)malloc(sizeof(double)*Nx);
double *y2x4 = (double*)malloc(sizeof(double)*Nx);

double *func_rq1a = (double*)malloc(sizeof(double)*Nr);
double *func_rq1b = (double*)malloc(sizeof(double)*Nr);
double *func_rq1c = (double*)malloc(sizeof(double)*Nr);
double *func_rq1d = (double*)malloc(sizeof(double)*Nr);
double *func_rq1e = (double*)malloc(sizeof(double)*Nr);

double *func_rq2a = (double*)malloc(sizeof(double)*Nr);
double *func_rq2b = (double*)malloc(sizeof(double)*Nr);
double *func_rq2c = (double*)malloc(sizeof(double)*Nr);
double *func_rq2d = (double*)malloc(sizeof(double)*Nr);
double *func_rq2e = (double*)malloc(sizeof(double)*Nr);

double *func_rq3a = (double*)malloc(sizeof(double)*Nr);
double *func_rq3b = (double*)malloc(sizeof(double)*Nr);
double *func_rq3c = (double*)malloc(sizeof(double)*Nr);
double *func_rq3d = (double*)malloc(sizeof(double)*Nr);
double *func_rq3e = (double*)malloc(sizeof(double)*Nr);

double *func_rq4a = (double*)malloc(sizeof(double)*Nr);
double *func_rq4b = (double*)malloc(sizeof(double)*Nr);
double *func_rq4c = (double*)malloc(sizeof(double)*Nr);
double *func_rq4d = (double*)malloc(sizeof(double)*Nr);
double *func_rq4e = (double*)malloc(sizeof(double)*Nr);

double *y2q1a = (double*)malloc(sizeof(double)*Nx);
double *y2q2a = (double*)malloc(sizeof(double)*Nx);
double *y2q3a = (double*)malloc(sizeof(double)*Nx);
double *y2q4a = (double*)malloc(sizeof(double)*Nx);

double *y2q1b = (double*)malloc(sizeof(double)*Nx);
double *y2q2b = (double*)malloc(sizeof(double)*Nx);
double *y2q3b = (double*)malloc(sizeof(double)*Nx);
double *y2q4b = (double*)malloc(sizeof(double)*Nx);

double *y2q1c = (double*)malloc(sizeof(double)*Nx);
double *y2q2c = (double*)malloc(sizeof(double)*Nx);
double *y2q3c = (double*)malloc(sizeof(double)*Nx);
double *y2q4c = (double*)malloc(sizeof(double)*Nx);

double *y2q1d = (double*)malloc(sizeof(double)*Nx);
double *y2q2d = (double*)malloc(sizeof(double)*Nx);
double *y2q3d = (double*)malloc(sizeof(double)*Nx);
double *y2q4d = (double*)malloc(sizeof(double)*Nx);

double *y2q1e = (double*)malloc(sizeof(double)*Nx);
double *y2q2e = (double*)malloc(sizeof(double)*Nx);
double *y2q3e = (double*)malloc(sizeof(double)*Nx);
double *y2q4e = (double*)malloc(sizeof(double)*Nx);


double *y2r1 = (double*)malloc(sizeof(double)*Nr);
double *y2r2 = (double*)malloc(sizeof(double)*Nr);
double *func_rR1 = (double*)malloc(sizeof(double)*Nr);
double *func_rR2 = (double*)malloc(sizeof(double)*Nr);

char *linear_power = cvector(0,200);
char *file_out = cvector(0,200);
char *mu_n_zero = cvector(0,200);
char *mu_n_one = cvector(0,200);
char *mu_n_two = cvector(0,200);
char *mu_n_three = cvector(0,200);
char *mu_n_four = cvector(0,200);

double *E00 = (double*)malloc(sizeof(double)*NKT);
double *E11 = (double*)malloc(sizeof(double)*NKT);
double *E12 = (double*)malloc(sizeof(double)*NKT);
double *E22 = (double*)malloc(sizeof(double)*NKT);
double *E23 = (double*)malloc(sizeof(double)*NKT);
double *E24 = (double*)malloc(sizeof(double)*NKT);
double *E33 = (double*)malloc(sizeof(double)*NKT);
double *E34 = (double*)malloc(sizeof(double)*NKT);
double *E44 = (double*)malloc(sizeof(double)*NKT);

double *mu_0 = (double*)malloc(sizeof(double)*NKT);
double *mu_1 = (double*)malloc(sizeof(double)*NKT);
double *mu_2 = (double*)malloc(sizeof(double)*NKT);
double *mu_3 = (double*)malloc(sizeof(double)*NKT);
double *mu_4 = (double*)malloc(sizeof(double)*NKT);

int main(int argc, char *argv[]) {

  double redshift;	
  sscanf(argv[1],"%lf",&redshift);

  //read in linear power spectrum
  for (int i=0; i<Nk_read; i++ ) P_lin[i] = mode[i] = y2[i] =0.0;
  sprintf(linear_power,"/users/angela/lustre/angela/DR10_power/CMASS_DR10_matterpower_z0.dat");
  int ind=0;
  FILE *fp_in;
  if((fp_in=fopen(linear_power,"r"))==NULL) printf("CAMB file not opened\n");
  const int bsz=200; char buf[bsz];  
  fgets(buf, bsz, fp_in);   
  while((fgets(buf, bsz, fp_in))!=NULL) { 
    if (buf[0] =='#'|| buf[0] ==0) continue; 
    double km, Pk;
    sscanf(buf,"%lf %lf\n", &km, &Pk);
    if(++(ind) > Nk_read) { ind--; break; }
    mode[ind] = km;
    P_lin[ind] = Pk;
  }
  fclose(fp_in);

  //define modes, x vals and r vals
  for (int i=0; i<NKT; i++) {
    E00[i]=E11[i]=E12[i]=0.0;
    E22[i]=E23[i]=E24[i]=0.0;
    E33[i]=E34[i]=E44[i]=0.0;
    P_red[i] = P_average[i] =k_out_vals[i] = 0.0;
    mu_0[i] = mu_1[i] = mu_2[i] = mu_3[i] = mu_4[i] =0.0;
  } 
  for (int i=0; i<Nx; i++) x_val[i] = i*0.02 -1.0;
  for (int i=0; i<Nr; i++) r_val[i] = i*0.1 + 0.00001;   //not sure about these values**TEST THIS***

  //read in files for spherical average integration
  sprintf(mu_n_zero,"/users/angela/lustre/angela/f0s.dat");
  int ind0=0;
  FILE *fp_mu_0;
  if((fp_mu_0=fopen(mu_n_zero,"r"))==NULL) printf("f0s file not opened\n");
  const int bsz_m0=200; char buf_m0[bsz_m0];  
  fgets(buf_m0, bsz_m0, fp_mu_0);   
  while((fgets(buf_m0, bsz_m0, fp_mu_0))!=NULL) { 
    if (buf_m0[0] =='#'|| buf_m0[0] ==0) continue; 
    double mu0;
    sscanf(buf_m0,"%lf\n", &mu0);
    if(++(ind0) > NKT) { ind0--; break; }
    mu_0[ind0] = mu0;
  }
  fclose(fp_mu_0); 

  sprintf(mu_n_one,"/users/angela/lustre/angela/f1s.dat");
  int ind1=0;
  FILE *fp_mu_1;
  if((fp_mu_1=fopen(mu_n_one,"r"))==NULL) printf("f1s file not opened\n");
  const int bsz_m1=200; char buf_m1[bsz_m1];  
  fgets(buf_m1, bsz_m1, fp_mu_1);   
  while((fgets(buf_m1, bsz_m1, fp_mu_1))!=NULL) { 
    if (buf_m1[0] =='#'|| buf_m1[0] ==0) continue; 
    double mu1;
    sscanf(buf_m1,"%lf\n", &mu1);
    if(++(ind1) > NKT) { ind1--; break; }
    mu_1[ind1] = mu1;
  }
  fclose(fp_mu_1);   

  sprintf(mu_n_two,"/users/angela/lustre/angela/f2s.dat");
  int ind2=0;
  FILE *fp_mu_2;
  if((fp_mu_2=fopen(mu_n_two,"r"))==NULL) printf("f2s file not opened\n");
  const int bsz_m2=200; char buf_m2[bsz_m2];  
  //fgets(buf_m2, bsz_m2, fp_mu_2);   
  while((fgets(buf_m2, bsz_m2, fp_mu_2))!=NULL) { 
    if (buf_m2[0] =='#'|| buf_m2[0] ==0) continue; 
    double mu2;
    sscanf(buf_m2,"%lf\n", &mu2);
    if(++(ind2) > NKT) { ind2--; break; }
    mu_2[ind2] = mu2;
  }
  fclose(fp_mu_2);  

  sprintf(mu_n_three,"/users/angela/lustre/angela/f3s.dat");
  int ind3=0;
  FILE *fp_mu_3;
  if((fp_mu_3=fopen(mu_n_three,"r"))==NULL) printf("f3s file not opened\n");
  const int bsz_m3=200; char buf_m3[bsz_m3];  
  //fgets(buf_m3, bsz_m3, fp_mu_3);   
  while((fgets(buf_m3, bsz_m3, fp_mu_3))!=NULL) { 
    if (buf_m3[0] =='#'|| buf_m3[0] ==0) continue; 
    double mu3;
    sscanf(buf_m3,"%lf\n", &mu3);
    if(++(ind3) > NKT) { ind3--; break; }
    mu_3[ind3] = mu3;
  }
  fclose(fp_mu_3);  

  sprintf(mu_n_four,"/users/angela/lustre/angela/f4s.dat");
  int ind4=0;
  FILE *fp_mu_4;
  if((fp_mu_4=fopen(mu_n_four,"r"))==NULL) printf("f4s file not opened\n");
  const int bsz_m4=200; char buf_m4[bsz_m4];  
  //fgets(buf_m4, bsz_m4, fp_mu_4);   
  while((fgets(buf_m4, bsz_m4, fp_mu_4))!=NULL) { 
    if (buf_m4[0] =='#'|| buf_m4[0] ==0) continue; 
    double mu4;
    sscanf(buf_m4,"%lf\n", &mu4);
    if(++(ind4) > NKT) { ind4--; break; }
    mu_4[ind4] = mu4;
  }
  fclose(fp_mu_4);  

  //spline the power spectrum.
  spline(mode, P_lin, Nk_read, 1.0e30, 1.0e30, y2); 

  //intergate the linear power spectrum
  double kmin=0.0000;
  double kmax =0.75;
  double PL_int = qsimp(function_lin, kmin,kmax);
  double f = 0.7441;
  double A = 1.0/(6.0*pi*pi)*PL_int;
  printf("A = %lf\n",A);

  double kx, ky, kz, k_mag, k_mag_sq;
  double kc =0.4;

  double Q1 =0.0;
  double Q2 =0.0;
  double Q3 =0.0;
  double Q4 =0.0;

  double epsilon =0.0001;//small value for integration

  for (int i=0; i<Nk; i++)  {
    printf("i=%d\n",i);
    for (int j=0; j<Nk; j++)  {
      for (int l=0; l<Nk; l++)  {

        int ink = l + Nk*( j + i*Nk);
        kx = (double) i*0.012 + 0.012;
        ky = (double) j*0.012 + 0.012;
        kz = (double) l*0.012 + 0.012;
        
        k_mag_sq = (kx*kx + ky*ky + kz*kz);
        k_mag = sqrt(k_mag_sq);
        if (k_mag_sq ==0) continue;
        k_out_vals[ink] = k_mag;
        double theta = acos(kz/k_mag);
        double phi = atan2(ky,kx);        
        double mu = 1.0/(k_mag*sqrt(3))*(kx + ky + kz);
        double damp = exp(-k_mag_sq*(1 + f*(f+2)*mu*mu)*A);
        double da = exp(-k_mag_sq*A);
        double PLk=0.0; 
        splint(mode, P_lin, y2, Nk_read, k_mag, &PLk);
        double lin_term = (1 + f*mu*mu)*(1 + f*mu*mu) *PLk;
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
       
            double q1 = r*r*(1-x*x)*(1-x*x);
            double q2 = (1 - x*x)*r*x*(1-r*x);
            double q3 = x*x*(1-r*x)*(1-r*x);          
            double q4 = 1 -x*x;

            double y = (1 + r*r -2*r*x);
            double arg = k_mag*sqrt(y);
            double PLy =0.0;
            splint(mode, P_lin, y2, Nk_read, arg, &PLy);

            func_x1[ix] = PLy*q1/(y*y); 
            func_x2[ix] = PLy*q2/(y*y); 
            func_x3[ix] = PLy*q3/(y*y);
            func_x4[ix] = PLy*q4/(y*y);
          }
          //have values of the integrand func_x, integrate over them
          double rs= r*r;
          spline(x_val, func_x1, Nx, 1.0e30, 1.0e30, y2x1);
          spline(x_val, func_x2, Nx, 1.0e30, 1.0e30, y2x2);
          spline(x_val, func_x3, Nx, 1.0e30, 1.0e30, y2x3);
          spline(x_val, func_x4, Nx, 1.0e30, 1.0e30, y2x4);

          double max_b = (1+ rs - epsilon*epsilon)/(2*r);
          double max_c = max_b;
          double min_e = (1 + rs -(kc/k_mag)*(kc/k_mag))/(2*r);

          double x1a_int = qsimp(function_q1, -1., 1.); 
          double x1b_int = qsimp(function_q1, -1., max_b);
          double x1c_int = qsimp(function_q1, -1., max_c);
          double x1d_int = x1a_int;
          double x1e_int = qsimp(function_q1, min_e, 1.);

          double x2a_int = qsimp(function_q2, -1., 1.);
          double x2b_int = qsimp(function_q2, -1., max_b);
          double x2c_int = qsimp(function_q2, -1., max_c);
          double x2d_int = x2a_int;
          double x2e_int = qsimp(function_q2, min_e, 1.);

          double x3a_int = qsimp(function_q3, -1., 1.);
          double x3b_int = qsimp(function_q3, -1., max_b);
          double x3c_int = qsimp(function_q3, -1., max_c);
          double x3d_int = x3a_int;
          double x3e_int = qsimp(function_q3, min_e, 1.);

          double x4a_int = qsimp(function_q4, -1., 1.);
          double x4b_int = qsimp(function_q4, -1., max_b);
          double x4c_int = qsimp(function_q4, -1., max_c);
          double x4d_int = x4a_int;
          double x4e_int = qsimp(function_q4, min_e, 1.);


          double Rr1 = -(2.0/(rs)) * (1.0 + rs) * (3.0 - 14.0*rs + 3.0*rs*rs) 
                        + (3.0/(r*rs))*(rs -1.0)*(rs -1.0)*(rs -1.0)*(rs -1.0) *log(fabs((1.0+r)/(1.0-r)));

          double Rr2 =  (2.0/(rs)) * (1.0 - rs) * (3.0 - 2.0*rs + 3.0*rs*rs) 
                        + (3.0/(r*rs))*(rs -1.0)*(rs -1.0)*(rs -1.0)*(rs +1.0) *log(fabs((1.0+r)/(1.0-r)));

          //then calculate the function(r) to integrate over
          double kr = k_mag*r;
          double PLr =0.0;
          splint(mode, P_lin, y2, Nk_read, kr, &PLr);

          func_rq1a[ir] = PLr * x1a_int;
          func_rq1b[ir] = PLr * x1b_int;
          func_rq1c[ir] = PLr * x1c_int;
          func_rq1d[ir] = PLr * x1d_int;
          func_rq1e[ir] = PLr * x1e_int;

          func_rq2a[ir] = PLr * x2a_int; 
          func_rq2b[ir] = PLr * x2b_int;
          func_rq2c[ir] = PLr * x2c_int;
          func_rq2d[ir] = PLr * x2d_int;
          func_rq2e[ir] = PLr * x2e_int;
             
          func_rq3a[ir] = PLr * x3a_int;
          func_rq3b[ir] = PLr * x3b_int;
          func_rq3c[ir] = PLr * x3c_int;
          func_rq3d[ir] = PLr * x3d_int;
          func_rq3e[ir] = PLr * x3e_int;

          func_rq4a[ir] = PLr * x4a_int;
          func_rq4b[ir] = PLr * x4b_int;
          func_rq4c[ir] = PLr * x4c_int;
          func_rq4d[ir] = PLr * x4d_int;
          func_rq4e[ir] = PLr * x4e_int;
          
          func_rR1[ir] = PLr * Rr1;
          func_rR2[ir] = PLr * Rr2;
        }

        double rmin = r_val[0];
        double rmax = r_val[Nr];
        //printf("r and x loops done\n");

        spline(r_val, func_rq1a, Nr, 1.0e30, 1.0e30, y2q1a);
        spline(r_val, func_rq2a, Nr, 1.0e30, 1.0e30, y2q2a);
        spline(r_val, func_rq3a, Nr, 1.0e30, 1.0e30, y2q3a);
        spline(r_val, func_rq4a, Nr, 1.0e30, 1.0e30, y2q4a);

        spline(r_val, func_rq1b, Nr, 1.0e30, 1.0e30, y2q1b);
        spline(r_val, func_rq2b, Nr, 1.0e30, 1.0e30, y2q2b);
        spline(r_val, func_rq3b, Nr, 1.0e30, 1.0e30, y2q3b);
        spline(r_val, func_rq4b, Nr, 1.0e30, 1.0e30, y2q4b);

        spline(r_val, func_rq1c, Nr, 1.0e30, 1.0e30, y2q1c);
        spline(r_val, func_rq2c, Nr, 1.0e30, 1.0e30, y2q2c);
        spline(r_val, func_rq3c, Nr, 1.0e30, 1.0e30, y2q3c);
        spline(r_val, func_rq4c, Nr, 1.0e30, 1.0e30, y2q4c);

        spline(r_val, func_rq1d, Nr, 1.0e30, 1.0e30, y2q1d);
        spline(r_val, func_rq2d, Nr, 1.0e30, 1.0e30, y2q2d);
        spline(r_val, func_rq3d, Nr, 1.0e30, 1.0e30, y2q3d);
        spline(r_val, func_rq4d, Nr, 1.0e30, 1.0e30, y2q4d);

        spline(r_val, func_rq1e, Nr, 1.0e30, 1.0e30, y2q1e);
        spline(r_val, func_rq2e, Nr, 1.0e30, 1.0e30, y2q2e);
        spline(r_val, func_rq3e, Nr, 1.0e30, 1.0e30, y2q3e);
        spline(r_val, func_rq4e, Nr, 1.0e30, 1.0e30, y2q4e);

        //printf("splines to all q functions done\n");
        double rmax_a = 1.- epsilon;
        double rmin_a = epsilon;

        double rmax_b = 1.;
        double rmin_b = rmax_a;

        double rmax_c = 1.+ epsilon;
        double rmin_c = 1.;

        double rmax_d = kc/k_mag -1;
        double rmin_d = rmax_c;

        double rmax_e = kc/k_mag;
        double rmin_e = rmax_d;
        
        double Q1a = factq * qsimp(function_qr1a, rmin_a, rmax_a);
        double Q1b = factq * qsimp(function_qr1b, rmin_b, rmax_b);
        double Q1c = factq * qsimp(function_qr1c, rmin_c, rmax_c);
        double Q1d = factq * qsimp(function_qr1d, rmin_d, rmax_d);
        double Q1e = factq * qsimp(function_qr1e, rmin_e, rmax_e);

        double Q2a = factq * qsimp(function_qr2a, rmin_a, rmax_a);
        double Q2b = factq * qsimp(function_qr2b, rmin_b, rmax_b);
        double Q2c = factq * qsimp(function_qr2c, rmin_c, rmax_c);
        double Q2d = factq * qsimp(function_qr2d, rmin_d, rmax_d);
        double Q2e = factq * qsimp(function_qr2e, rmin_e, rmax_e);

        double Q3a = factq * qsimp(function_qr3a, rmin_a, rmax_a);
        double Q3b = factq * qsimp(function_qr3b, rmin_b, rmax_b);
        double Q3c = factq * qsimp(function_qr3c, rmin_c, rmax_c);
        double Q3d = factq * qsimp(function_qr3d, rmin_d, rmax_d);
        double Q3e = factq * qsimp(function_qr3e, rmin_e, rmax_e);

        double Q4a = factq * qsimp(function_qr4a, rmin_a, rmax_a);
        double Q4b = factq * qsimp(function_qr4b, rmin_b, rmax_b);
        double Q4c = factq * qsimp(function_qr4c, rmin_c, rmax_c);
        double Q4d = factq * qsimp(function_qr4d, rmin_d, rmax_d);
        double Q4e = factq * qsimp(function_qr4e, rmin_e, rmax_e);

        double Q1 = Q1a + Q1b + Q1c + Q1d + Q1e;
        double Q2 = Q2a + Q2b + Q2c + Q2d + Q2e;
        double Q3 = Q3a + Q3b + Q3c + Q3d + Q3e;
        double Q4 = Q4a + Q4b + Q4c + Q4d + Q4e;

        //printf("Q1 = %lf\n",Q1);
        //printf("Q2 = %lf\n",Q2);
        //printf("Q3 = %lf\n",Q3);
        //printf("Q4 = %lf\n",Q4);

        //printf("integrations for all q functions done\n");

        spline(r_val, func_rR1, Nr, 1.0e30, 1.0e30, y2r1);
        spline(r_val, func_rR2, Nr, 1.0e30, 1.0e30, y2r2);

        //printf("splines to all r functions done\n");

        double rmax1 = kc/k_mag;
        double R1a = (1./48.) * PLk * factq * qsimp(function_r1, rmin_a, rmax_a);
        double R1b = (1./48.) * PLk * factq * qsimp(function_r1, rmin_b, rmax_b);
        double R1c = (1./48.) * PLk * factq * qsimp(function_r1, rmin_c, rmax_c);
        double R1d = (1./48.) * PLk * factq * qsimp(function_r1, rmin_d, rmax_d);
        double R1e = (1./48.) * PLk * factq * qsimp(function_r1, rmin_e, rmax_e);
        double R1 = R1a;// + R1b + R1c + R1d + R1e;
        
        //printf("R1 done\n");

        double R2a = (1./48.) * PLk * factq * qsimp(function_r2, rmin_a, rmax_a);
        double R2b = (1./48.) * PLk * factq * qsimp(function_r2, rmin_b, rmax_b);
        double R2c = (1./48.) * PLk * factq * qsimp(function_r2, rmin_c, rmax_c);
        double R2d = (1./48.) * PLk * factq * qsimp(function_r2, rmin_d, rmax_d);
        double R2e = (1./48.) * PLk * factq * qsimp(function_r2, rmin_e, rmax_e);
        double R2 = R2a + R2b + R2c + R2d + R2e;
        //printf("R2 done\n");

        //printf("R1 = %lf\n",R1);
        //printf("R2 = %lf\n",R2);      
       // printf("now sum up all the contributions\n");
        E00[ink] = (9./98.)*Q1 + (3./7.)  *Q2 + (1./2.)*Q3 +              (10./21.)*R1 + (6./7.) *R2;
        E11[ink] = 4.0*E00[ink];
        E12[ink] = -(3./14)*Q1 - (3./2.)  *Q2 +              (1./4.)*Q4 - (6./7.)  *R1;
        E22[ink] = (57./98)*Q1 + (51./14.)*Q2 + (3.0)  *Q3 - (1./4.)*Q4 + (16./7.) *R1 + (30./7.)*R2;
        E23[ink] = -(3./7.)*Q1 - (3.0)    *Q2 +              (1./2.)*Q4 - (6./7.)  *R1;
        E24[ink] = (3./16) *Q1;
        E33[ink] = (3./7.) *Q1 + (27./7.) *Q2 + (2.0)  *Q3 - (1./2.)*Q4 + (6./7.)  *R1 + (12./7.)*R2;
        E34[ink] = -(3./8.)*Q1 - (3./2.)  *Q2 +              (1./4.)*Q4;
        E44[ink] = (3./16.)*Q1 + (3./2.)  *Q2 + (1./2.)*Q3 - (1./4.)*Q4;

        double mus = mu*mu;
        double other_terms = 
        E00[ink] 
        + (mus * f * E11[ink])
        + (mus * f*f * E12[ink])
        + (mus*mus *f*f * E22[ink])
        + (mus*mus *f*f*f * E23[ink])
        + (mus*mus *f*f*f*f * E24[ink])
        + (mus*mus*mus *f*f*f * E33[ink])
        + (mus*mus*mus *f*f*f*f * E34[ink])
        + (mus*mus*mus*mus *f*f*f*f *E44[ink]);
        double ratio_lin = lin_term/other_terms;

        double av_E00 = Spi*da*mu_0[ink]*E00[ink];
        double av_E11 = Spi*da*mu_1[ink]*E11[ink];
        double av_E12 = Spi*da*mu_1[ink]*E12[ink];
        double av_E22 = Spi*da*mu_2[ink]*E22[ink];
        double av_E23 = Spi*da*mu_2[ink]*E23[ink];
        double av_E24 = Spi*da*mu_2[ink]*E24[ink];
        double av_E33 = Spi*da*mu_3[ink]*E33[ink];
        double av_E34 = Spi*da*mu_3[ink]*E34[ink];
        double av_E44 = Spi*da*mu_4[ink]*E44[ink];

        double av_other_terms = av_E00 + av_E11 + av_E12 + av_E22 + av_E23 + av_E24
                              + av_E33 + av_E34 + av_E44;
        //if (ratio_lin<1) printf("kmag = %lf, ratio = %lf\n",k_mag,ratio_lin);
        //printf("damp = %lf\n",damp);
        //printf("kmag = %lf, linear terms = %lf\n",k_mag, lin_term);
        //printf("other terms = %lf\n",other_terms);
        P_red[ink] = damp*(lin_term + other_terms);
        //spherically average using Legendre polynomials
        //double Legendre_2 = 0.5*(3*mu*mu -1);
        //P_average[ink] = (5.0)*(P_red[ink]) *Legendre_2;
        P_average[ink] = av_lin_a + av_lin_b + av_lin_c + av_other_terms;
      }
    }
  }

  sprintf(file_out,"/users/angela/lustre/angela/LRT_redshift_power_av6.txt",redshift);
  FILE *fp_out;
  if((fp_out= fopen(file_out,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=0 ;i<Nk; i++) 
    for (int j=0; j<Nk; j++)
      for (int l=0; l<Nk; l++) {
        int ind = l + Nk*( j + Nk*i);
        //kx = (double) i*0.012 + 0.012;
        //ky = (double) j*0.012 + 0.012;
        //kz = (double) l*0.012 + 0.012;       
        k_mag_sq = (kx*kx + ky*ky + kz*kz);
        k_mag = k_out_vals[ind];//sqrt(k_mag_sq);
        double muT =  1.0/(k_mag*sqrt(3))*(kx + ky + kz);
        fprintf(fp_out,"%lf %lf\n",k_mag, P_average[ind]);
        //fprintf(fp_out,"%lf\n",k_mag);
  }

}
//end of main

//----------------------------------------------------------
double function_lin (double x) {
  double value=0;
  splint(mode, P_lin, y2, Nk_read, x, &value);
}
double function_q1 (double x) {
  double value=0;
  splint(x_val, func_x1, y2x1, Nx, x, &value);
}
double function_q2 (double x) {
  double value=0;
  splint(x_val, func_x2, y2x2, Nx, x, &value);
}
double function_q3 (double x) {
  double value=0;
  splint(x_val, func_x3, y2x3, Nx, x, &value);
}
double function_q4 (double x) {
  double value=0;
  splint(x_val, func_x4, y2x4, Nx, x, &value);
}
//---------------------------------------------------------
double function_qr1a (double x) {
  double value=0;
  splint(r_val, func_rq1a, y2q1a, Nr, x, &value);
}
double function_qr2a (double x) {
  double value=0;
  splint(r_val, func_rq2a, y2q2a, Nr, x, &value);
}
double function_qr3a (double x) {
  double value=0;
  splint(r_val, func_rq3a, y2q3a, Nr, x, &value);
}
double function_qr4a (double x) {
  double value=0;
  splint(r_val, func_rq4a, y2q4a, Nr, x, &value);
}
//--
double function_qr1b (double x) {
  double value=0;
  splint(r_val, func_rq1b, y2q1b, Nr, x, &value);
}
double function_qr2b (double x) {
  double value=0;
  splint(r_val, func_rq2b, y2q2b, Nr, x, &value);
}
double function_qr3b (double x) {
  double value=0;
  splint(r_val, func_rq3b, y2q3b, Nr, x, &value);
}
double function_qr4b (double x) {
  double value=0;
  splint(r_val, func_rq4b, y2q4b, Nr, x, &value);
}
//----
double function_qr1c (double x) {
  double value=0;
  splint(r_val, func_rq1c, y2q1c, Nr, x, &value);
}
double function_qr2c (double x) {
  double value=0;
  splint(r_val, func_rq2c, y2q2c, Nr, x, &value);
}
double function_qr3c (double x) {
  double value=0;
  splint(r_val, func_rq3c, y2q3c, Nr, x, &value);
}
double function_qr4c (double x) {
  double value=0;
  splint(r_val, func_rq4c, y2q4c, Nr, x, &value);
}
//------
double function_qr1d (double x) {
  double value=0;
  splint(r_val, func_rq1d, y2q1d, Nr, x, &value);
}
double function_qr2d (double x) {
  double value=0;
  splint(r_val, func_rq2d, y2q2d, Nr, x, &value);
}
double function_qr3d (double x) {
  double value=0;
  splint(r_val, func_rq3d, y2q3d, Nr, x, &value);
}
double function_qr4d (double x) {
  double value=0;
  splint(r_val, func_rq4d, y2q4d, Nr, x, &value);
}
//----
double function_qr1e (double x) {
  double value=0;
  splint(r_val, func_rq1e, y2q1e, Nr, x, &value);
}
double function_qr2e (double x) {
  double value=0;
  splint(r_val, func_rq2e, y2q2e, Nr, x, &value);
}
double function_qr3e (double x) {
  double value=0;
  splint(r_val, func_rq3e, y2q3e, Nr, x, &value);
}
double function_qr4e (double x) {
  double value=0;
  splint(r_val, func_rq4e, y2q4e, Nr, x, &value);
}
//----------------------------------------------------------

double function_r1 (double x) {
  double value=0;
  splint(r_val, func_rR1, y2r1, Nr, x, &value);
}
double function_r2 (double x) {
  double value=0;
  splint(r_val, func_rR2, y2r2, Nr, x, &value);
}

double mu1 (double x, int flag) {
  double val =0.0;
  double exx= exp(-x);
  double stx = sqrt(x);
  double spi = sqrt(pi);
  if (flag==0) {
  }
  if (flag ==1) {
    val = -(1/2)*spi*(exx/(spi*stx) - erf(stx)/(2*stx*stx*stx));
  }
  if (flag ==2) {
    val = (1/2)*spi*( -exx/(spi*x*x) + ((-exx/(2*spi*stx*stx*stx)) - (exx/(spi*stx)))/stx + 3*erf(stx)/(4*stx*stx*stx*stx*stx));
  }
  if (flag ==3) {
    val = -(1/2)*spi*((9*exx)/(4*spi*x*x*x) - 3*((-exx/(2*spi*stx*stx*stx)) - (exx/(spi*stx)))/(2*stx*stx*stx) 
        + ((3*exx/(4*spi*stx*stx*stx*stx*stx)) + (exx/(spi*stx*stx*stx)) + (exx/(spi*stx)))/(stx) - 15*erf(stx)/(8*x*x*x*stx));
  }
  if (flag ==4) {
     val = (1/2)* spi* (-((15*exx)/(2*spi* x*x*x*x)) + (
   9* (-(exx/(2 *spi*stx*stx*stx)) - exx/(spi*stx)))/(
   2*stx*stx*stx*stx*stx) - (
   2 *((3*exx)/(4 *spi*stx*stx*stx*stx*stx) + exx/(spi*stx*stx*stx) + 
      exx/(spi*stx)))/(stx*stx*stx) + (-((15*exx)/(8*spi*x*x*x*stx)) - (9*exx)/(
    4*spi* x*x*stx) - (3*exx)/(2*spi*x*stx) - exx/(
    spi*stx))/stx + (105*erf(stx))/(16* x*x*x*x*stx)) ;
  }
    
}

/**********************/
 

      
 

