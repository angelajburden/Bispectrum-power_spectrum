//Calculates the power spectrum in realshift space with LRT .
//uses the r1,r2,q1,q2..integrals from a program run in Mathematica. These need to be done seperately for
//different redshifts or a different selection of k vals.
//g++ LRT_real_space.c util.c integration.c -lsrfftw -lsfftw -lm -o LRT_real

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
#define eps 1.0e-5 // numerical accuracy for integrations
#define JMAX 20
#define NEVAL 10000
#define CHUNKSIZE 1; 

//const double redshift_no = 0.0;//change to 0.57
 
int Nk =10; //number of k bins
int Nr =100;
int Nx =100;
int Nk_read =450;// 564 for camb_pk_linear_057.dat;//450 matterpower_z0_k08.dat;
//int Nk_read =564;
int NKT = 60;

//declare functions
double function_lin (double);
double mu1 (double, int);
void read_in_integrals(char*, double*, int);
void read_in_integralsR(char*, double*, int);
double *P_real = (double*)malloc(sizeof(double)*NKT);
double *P1 = (double*)malloc(sizeof(double)*NKT);
double *P2 = (double*)malloc(sizeof(double)*NKT);
double *P3 = (double*)malloc(sizeof(double)*NKT);
double *P4 = (double*)malloc(sizeof(double)*NKT);
double *P5 = (double*)malloc(sizeof(double)*NKT);
double *lin_only = (double*)malloc(sizeof(double)*NKT);
double *k_out_vals = (double*)malloc(sizeof(double)*NKT);
double *mode = (double*)malloc(sizeof(double)*Nk_read);

double *P_lin = (double*)malloc(sizeof(double)*Nk_read);
double *y2 = (double*)malloc(sizeof(double)*Nk_read);

char *linear_power = cvector(0,200);
char *file_out = cvector(0,200);
char *Q1_file = cvector(0,200);
char *Q2_file = cvector(0,200);
char *Q3_file = cvector(0,200);
char *Q4_file = cvector(0,200);
char *R1_file = cvector(0,200);
char *R2_file = cvector(0,200);

double *Q1 = (double*)malloc(sizeof(double)*NKT);
double *Q2 = (double*)malloc(sizeof(double)*NKT);
double *Q3 = (double*)malloc(sizeof(double)*NKT);
double *Q4 = (double*)malloc(sizeof(double)*NKT);
double *R1 = (double*)malloc(sizeof(double)*NKT);
double *R2 = (double*)malloc(sizeof(double)*NKT);

const double Om = 0.274;//today
const double OL = 0.726;
const double a =1.;
const double H = 70.;
double D_fact = 0;

int main(int argc, char *argv[]) {

  double redshift;	
  sscanf(argv[1],"%lf",&redshift);

  //define modes, x vals and r vals
  for (int i=0; i<NKT; i++) {
    Q1[i]=Q2[i]=Q3[i]=Q4[i]=0.0;
    R1[i]=R2[i]=0.0;
    P_real[i] =k_out_vals[i] = 0.0;
    P1[i] = P2[i] = P3[i] = P4[i] = P5[i] = lin_only[i] = 0.0;
  } 
  double growth1 = Growth(redshift)/0.994755;
  printf("growth1=%lf\n",growth1*growth1);
  D_fact = growth1*growth1;

  //read in linear power spectrum
  for (int i=0; i<Nk_read; i++ ) P_lin[i] = mode[i] = y2[i] =0.0;
  sprintf(linear_power,"/users/angela/lustre/angela/LRT_stuff/matterpower_z0_k08.dat");
  int ind=0;
  FILE *fp_in;
  if((fp_in=fopen(linear_power,"r"))==NULL) printf("CAMB file not opened\n");
  const int bsz=200; char buf[bsz];  
  //fgets(buf, bsz, fp_in);   
  while((fgets(buf, bsz, fp_in))!=NULL) { 
    if (buf[0] =='#'|| buf[0] ==0) continue; 
    double km, Pk;
    sscanf(buf,"%lf %lf\n", &km, &Pk);
    if(++(ind) > Nk_read) { ind--; break; }
    mode[ind] = km;
    P_lin[ind] = Pk*D_fact;
  }
  fclose(fp_in);

  char *q1a = "q1_all";
  char *q2a = "q2_all";
  char *q3a = "q3_all";
  char *q4a = "q4_all";
  char *R1a = "R1_all";
  char *R2a = "R2_all";

  read_in_integrals(q1a,Q1,NKT);
  read_in_integrals(q2a,Q2,NKT);
  read_in_integrals(q3a,Q3,NKT);
  read_in_integrals(q4a,Q4,NKT);

  read_in_integralsR(R1a,R1,NKT);
  read_in_integralsR(R2a,R2,NKT);

  //spline the power spectrum.
  spline(mode, P_lin, Nk_read, 1.0e30, 1.0e30, y2); 

  //integrate the linear power spectrum
  double kmin=0.0;
  double kmax =0.80;
  double PL_int = qsimp(function_lin, kmin,kmax);
  double A = 1.0/(6.0*pi*pi)*PL_int;
  printf("A = %lf\n",A);

  double kx, ky, kz, k_mag, k_mag_sq;

  int NK=200;

  for (int i=0;i<NKT; i++) {
        int ink =i;
        printf("%d\n",i);
        double k_mag = i*0.005 + 0.005;
        if (k_mag ==0) continue;      
        double da = exp(-k_mag*k_mag*A);
        double PLk=0.0; 
        splint(mode, P_lin, y2, Nk_read, k_mag, &PLk);
        double lin_part = da*PLk;
        double factq = (k_mag*k_mag*k_mag/(4.*pi*pi))*D_fact*D_fact;
        double factr = (1/48)*(k_mag*k_mag*k_mag/(4.*pi*pi))*D_fact*D_fact;
        P1[ink] = da*factq*(9./98.)*Q1[i];
        P2[ink] = da*factq*(3./7.)  *Q2[i];
        P3[ink] = da*factq*(1./2.)*Q3[i];
        P4[ink] = da*factq*(10./21.)*R1[i];
        P5[ink] = da*factq*(6./7.)*R2[i];  
        P_real[ink] = lin_part + P1[ink] +P2[ink] +P3[ink] +P4[ink] +P5[ink] ;
        lin_only[ink] = lin_part;
  }
  int RED_file =(int) (redshift*10);
  if (redshift ==0.0) RED_file =0;
  printf("REDFILE = %d\n",RED_file);
  sprintf(file_out,"/users/angela/lustre/angela/LRT_stuff/power_spectra/LRT_real_power_z%d.txt",RED_file);
  FILE *fp_out;
  if((fp_out= fopen(file_out,"w"))==0) 
    printf("cannot open pow output file\n"); 
  for (int i=0;i<NKT; i++) {
        double k_mag = i*0.005 + 0.005;
        //printf("kmag = %lf\n",k_mag);
        if (k_mag ==0) continue; 
        fprintf(fp_out,"%lf %lf\n",k_mag,P_real[i]);
        //fprintf(fp_out,"%lf\n",k_mag);
  }

}
//end of main

//----------------------------------------------------------
double function_lin (double x) {
  double value=0;
  splint(mode, P_lin, y2, Nk_read, x, &value);
}

//----------------------------------------------------------

double mu1 (double x, int flag) {//this function calculates the spherically averaged terms
  double val =0.0;
  double val1 =0.0;
  double val2 =0.0;
  double val3 =0.0;
  double val4 =0.0;
  double exx= exp(-x);
  double stx = sqrt(x);
  double spi = sqrt(pi);
  if (flag==0) {
    val = (spi*erf(stx))/(2.*stx);
    return val;
  }
  if (flag ==1) {
    val1 = -(1./2.)*spi*(exx/(spi*x) - erf(stx)/(2.*stx*stx*stx));
    //printf("VAL=%lf, flag =%d\n",val1,flag);
    return val1;
  }
  if (flag ==2) {
    val2 = (1./2.)*spi*( -exx/(spi*x*x) + ((-exx/(2.*spi*stx*stx*stx)) - (exx/(spi*stx)))/stx + 3.*erf(stx)/(4.*x*x*stx));
    return val2;
  }
  if (flag ==3) {
    val3 = -(1./2.)*spi*((9.*exx)/(4.*spi*x*x*x) - 3.*((-exx/(2.*spi*x*stx)) - (exx/(spi*stx)))/(2.*x*stx) 
        + ((3.*exx/(4.*spi*x*x*stx)) + (exx/(spi*x*stx)) + (exx/(spi*stx)))/(stx) - 15.*erf(stx)/(8.*x*x*x*stx));
    return val3;
  }
  if (flag ==4) {
     val4 = (1./2.)* spi* (-((15.*exx)/(2.*spi* x*x*x*x)) + (
   9.* (-(exx/(2. *spi*stx*stx*stx)) - exx/(spi*stx)))/(
   2.*stx*stx*stx*stx*stx) - (
   2. *((3.*exx)/(4. *spi*stx*stx*stx*stx*stx) + exx/(spi*stx*stx*stx) + 
      exx/(spi*stx)))/(stx*stx*stx) + (-((15.*exx)/(8.*spi*x*x*x*stx)) - (9.*exx)/(
    4.*spi* x*x*stx) - (3.*exx)/(2.*spi*x*stx) - exx/(
    spi*stx))/stx + (105.*erf(stx))/(16.* x*x*x*x*stx)) ;
    return val4;
  }
    
}

/**********************/
 
void read_in_integrals(char *Q_val, double *Qx, int NKT){
  char *file =cvector(0,200); 
  sprintf(file,"/users/angela/lustre/angela/LRT_stuff/Q_and_R_z0/%sTz0.dat",Q_val);
  printf("file = %s\n",file);
  int ind=0;
  FILE *fp_q;
  if((fp_q=fopen(file,"r"))==NULL) printf("%s not opened\n",file);
  const int bsz_q=200; char buf_q[bsz_q];    
  while((fgets(buf_q, bsz_q, fp_q))!=NULL) { 
    if (buf_q[0] =='#'|| buf_q[0] ==0) continue; 
    double qx;
    sscanf(buf_q,"%lf\n", &qx);
    if(++(ind) > NKT) { ind--; break; }
    Qx[ind] = qx;
  }
  fclose(fp_q);
} 
      
 void read_in_integralsR(char *Q_val, double *Qx, int NKT){
  char *file =cvector(0,200); 
  sprintf(file,"/users/angela/lustre/angela/LRT_stuff/Q_and_R_z0/%sTz0.dat",Q_val);
  printf("file = %s\n",file);
  int ind=0;
  FILE *fp_q;
  if((fp_q=fopen(file,"r"))==NULL) printf("%s not opened\n",file);
  const int bsz_q=200; char buf_q[bsz_q];    
  while((fgets(buf_q, bsz_q, fp_q))!=NULL) { 
    if (buf_q[0] =='#'|| buf_q[0] ==0) continue; 
    double qx;
    sscanf(buf_q,"%lf\n", &qx);
    if(++(ind) > NKT) { ind--; break; }
    Qx[ind] =(1./48.)*qx;
  }
  fclose(fp_q);
}

double Growth(double red) {
  double Ok = 1. - Om_0 - OL_0;
  double H2 = Om_0*(1.+red)*(1.+red)*(1.+red) + Ok*(1.+red)*(1.+red) + OL_0;
  double Om = Om_0*(1.+red)*(1.+red)*(1.+red)/H2;
  double OL = OL_0/H2;
  double pow_val = exp((0.57142857)*log(Om));//4./7.
  double D_z = (5.*Om)/((2.*(pow_val)-OL + (1. + 0.5*Om)*(1. + (OL/70.)))*(1.+red));
  double f_val = exp((0.57142857)*log(Om)) + (1 + Om/2)*OL/70;
  return D_z;//f_val;//
}
