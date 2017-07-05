#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <string>
//#include <cpgplot.h>

const int   NR_END     = 1;  

const double pi    = 3.1415926536; // Pi
const double G_cst = 6.67e-11;     // gravitational cst G
const double Mpc   = 3.09e22;      // Mpc in meters
const double M_sun = 1.99e30;      // Solar mass in kg
const double MSTAR = -20.44;

//cosmological parameters
const double Om_0 = 0.27;
const double OL_0 = 1.0-Om_0;
const double a_0 =1;
const double H0 = 70;
const double bias = 1.85;
const double f_growth = 0.7441;//exp(0.6*log(Om_0)) + (OL_0/70.)*( 1 + 0.5*Om_0);

// define universe we're working in - read from file
extern double sig8;         // value of sigma_8
extern double rho;          // current density of universe (solar mass/Mpc^3)
extern double G0;           // present day density cst
extern double Gl;           // Cosmological Constant densty
//extern double H0;           // Hubble cst km/s.MPc
extern double pow_n;        // Scalar spectral index
extern double bfrac;        // Baryon fraction of total density
extern double pow_norm;     // power spectrum normalisation z=0 only
extern double a_st;         // Sheth & Tormen constant a
extern double p_st;         // Sheth & Tormen constant p
extern double N_st;         // Sheth & Tormen constant N
extern double gom_m;        // time varying Omega_m
extern double gom_v;        // time varying Omega_v



extern unsigned long ija[];  //for spare matrix alogrithm and 
extern double sa[];          //precon. conjugate grad. method

// set up for complex structures
#ifndef _FCOMPLEX_DECLARE_T_
  typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif



// house keeping functions: util.c
void   err_handler(const char*);
double *dvector(long,long);
float  poidev(float, long*);
void   free_dvector(double*,long,long);
float  *vector(long,long);
void   free_vector(float*,long,long);
int    *ivector(long,long);
void   free_ivector(int*,long,long);
char   *cvector(long,long);
void   free_cvector(char*,long,long);
float  **matrix(long,long,long,long);
void   free_matrix(float**,long,long,long,long);
int    **imatrix(long,long,long,long); 
void   free_imatrix(int**,long,long,long,long);
double **dmatrix(long,long,long,long);
void   free_dmatrix(double**,long,long,long,long);
float  ***f3tensor(long,long,long,long,long,long);
void   free_f3tensor(float***,long,long,long,long,long,long);
float  ****f4tensor(long,long,long,long,long,long,long,long);
void   free_f4tensor(float****,long,long,long,long,long,long,long,long);
int    **vimatrix(long,long,long,int*);
void   free_vimatrix(int**,long,long,long);
float  ran2(long*);
float  gasdev(long*);

// chisq optimisation routines: fit_param.c
double golden(double,double,double,double(*f)(double),double,double*);
void powell(double[],double**,int,double,int*,double*,double (*func)(double []));

// integration functions: integration.c
double qsimp(double (*func)(double),double,double);
double qsimp2(double (*func)(double),double,double);
double qmidinf(double (*func)(double,double),double,double,double);
double qmidinf2(double (*func)(double,double,double),double,double,double,double);
double qsimpmid(double (*func)(double,double,double),double,double,double,double);
double qsimpmid2(double (*func)(double,double),double,double,double);

// root finding functions: root.c
double dfridr(double (*func)(double,double),double,double,double*,double);
double zbrent(double (*func)(double,double),double,double,double,double);
void   zroots(float a[], int, fcomplex roots[], int);
void   laguer(fcomplex a[], int, fcomplex*, int*);

// spline fitting routines: spline.c
void   spline(double[],double[],int,double,double,double[]);
void   splint(double[],double[],double[],int,double,double*);

//  matrix inversion: mat_inv.c
void mat_inv(double**,double**,int,double*);
void dsvdcmp(double**,int,int,double[],double**);

// PS mass functions: ps_mass_func.c
void   nm_power();
double Pk(double); 
double Sigth(double,double);
double dSigthbydR(double,double);

//preconditioned conjugate gradient method for multiplying matrices: linbcg.c
void linbcg(unsigned long,double,double,int,double,int,int,double,double,unsigned long);
void atimes(unsigned long,double,double,int,double,unsigned long,double);
void asolve(unsigned long,double,double,int);
void dsprsax(double,unsigned long,double,double,unsigned long);
void dsprstx(double,unsigned long,double,double,unsigned long);
//void sprsin(float**,int,float,unsigned long,float,unsigned long);
double snrm(unsigned long,double,int);

//power spectrum calculation: power.c
//void calc_power(float*,int,int,double, char*,double,double,double);
void fun_calc_power(double*,int,int,double, char*);
void fun_calc_bispectrum(double*,int,long int,int,double,char*);
//binning of data
double calc_dp(double);//calculate distance from redshift
double Hubble(double);//calculate Hubble parameter
double Growth(double);//calculate growth function.

//functions io_gadget_slim.c                                                                               
struct io_header get_gadget_header(const char *fname);
int get_Numpart(struct io_header *header);
//void read_in_binary_gadget(struct basic_gal*, char*, char*, int, char*, long int*);
void read_in_binary_gadget(double del[], int, double, char*, char*, int, char*, long int*);
struct basic_gal { 
    float cp[3]; 
  //  double vel[3];
  //  double id; 
}; 

void calculate_fkp(double, double, double, int, struct basic_gal*, float nbar_slice[],int, float*, float*, double*, double*, double*, int);
void bin_gals_NGP(struct basic_gal*, double, double, double, int, float dr[], float dr_power[], int, double, int, double);
void smooth_field(double, float dr[], float dr_uni[], int, double, int);
void calc_redshift(struct basic_gal*, int, double*, double*);
void read_in_binary_mock(const char*, struct basic_gal*, long int*, long int);
void xy_plot(double*, double* , int, int);
void print_time(time_t, time_t, const char*);
