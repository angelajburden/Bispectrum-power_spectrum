//module for bispectrum functions
#include "header.h"
#include <stdio.h>
#include <rfftw.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fstream>
#include <iomanip>
#include<complex>

#define ASEED 12346653
// function to calculate bispectrum of cubic grid and output to file

void fun_calc_bispectrum(double *del, int NG, long int NGAL, int NB, double BOX, char *name) {
  const int NGMAX = NG*NG*2*(NG/2+1);
  rfftwnd_plan dp1;
  dp1 = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)-1,
			    FFTW_MEASURE | FFTW_IN_PLACE);
  printf("fourier transform\n");
  rfftwnd_one_real_to_complex(dp1,(fftw_real*)del,NULL);
  
  //free memory                                                                                
  rfftwnd_destroy_plan(dp1);
  printf("memory for r2c plan freed\n");

  double k_ny = pi*(double)NG/2.*BOX;
  double k_fund = 2.*pi/BOX;
   
  // Set up bins for P(k) data 

  const double mink = 0.00;                // minimum k to be binned
  //  const double maxk = k_ny;//2.*pi*0.25;        // maximum k to be binned
  double binwidth=k_fund;
  int num_bins = (NB*NB*NB);
  double *Bg;
  
  Bg = (double *)malloc(sizeof(double) * num_bins);
  for(int i=0;i<num_bins; i++){
     Bg[i]=0.0;
   }
  fprintf(stderr,"allocate 2d array memory\n");
  double **arr = (double **)malloc(NB*sizeof(double*));
  int *Nkcount = (int *) malloc (NB*sizeof(int*));
  fprintf(stderr,"loop through to allocate memory for 2d\n");
  for (int i=0; i<NB; ++i){
    arr[i] =(double *) malloc(NGMAX*sizeof(double));
    Nkcount[i]=0;
    if (arr[i]==NULL) printf("memory allocation problem\n");  /* Error checking */
  }
  for (int i=0; i<NB; i++) {
    for (int j=0; j<NGMAX; j++) {
      arr[i][j]=0.0;
    }
  }
  // Note that arr[i][j] is same as *(*(arr+i)+j)
  double fktot_prev =0.0;
  // Work out |k| for each frequency component
  double fx, fy, fz;
  fprintf(stderr,"loop though grid\n");
  //now loop through grid checking for these three modes
  for(int i=0;i<NG;i++) {
  // frequency in x
    if(i<=NG/2) fx = (double)i/BOX; 
    else        fx = ((double)i-(double)NG)/BOX;

    for(int j=0;j<NG;j++) {
      // frequency in y
      if(j<=NG/2) fy = (double)j/BOX; 
      else        fy = ((double)j-(double)NG)/BOX;
	  
      for(int k=0;k<=NG/2;k++) {
        // frequency in z
	fz=(double)k/BOX;
	    
	// length of k vector and bin
	double fktot=2.*pi*sqrt(fx*fx+fy*fy+fz*fz);
	int bin = 1+ (int)( (double)(fktot-mink)/binwidth);
	    
        if( fktot < k_ny && fktot > 0.00 && bin < NB) {

	  double deltak = fabs(fktot - fktot_prev);
          fktot_prev=fktot;
          arr[bin][(2*k  )+2*(NG/2+1)*(j+NG*i)]+=del[(2*k  )+2*(NG/2+1)*(j+NG*i)]*binwidth;
	  arr[bin][(2*k+1)+2*(NG/2+1)*(j+NG*i)]+=del[(2*k+1)+2*(NG/2+1)*(j+NG*i)]*binwidth; 
   	  Nkcount[bin] ++;        
	}
      }
    }
  } 
  fprintf(stderr,"reverse fftw\n");
  rfftwnd_plan dp_c2r1;
  dp_c2r1 = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)+1,
                                FFTW_MEASURE | FFTW_IN_PLACE);
  for (int i=0; i<NB; i++) {
    rfftwnd_one_complex_to_real(dp_c2r1,(fftw_complex*)arr[i],NULL);
      for (int j=0; j<(NGMAX); j++){
        if( arr[i][j] <0.1) arr[i][j]=0.0;
        arr[i][j] /= (double)(NG*NG*NG);
      }
  }
  fprintf(stderr, "destroy RFFTW plan\n");
  rfftwnd_destroy_plan(dp_c2r1);
  fprintf(stderr,"bin the bispectrum triplets\n");
  double delR = (BOX/(double)NG)*(BOX/(double)NG)*(BOX/(double)NG);
  for (int ibin=1; ibin<NB; ibin++){
    for (int jbin=(ibin+1); jbin<NB; jbin++) {
      for (int kbin=(jbin+1); kbin<NB; kbin++) {
	
        int iii = ibin + NB*(ibin + NB*ibin);
	int ijj = ibin + NB*(jbin + NB*jbin);
        int jij = jbin + NB*(ibin + NB*jbin);
        int jji = jbin + NB*(jbin + NB*ibin);
        int iij = ibin + NB*(ibin + NB*jbin);
        int iji = ibin + NB*(jbin + NB*ibin);
        int jii = jbin + NB*(ibin + NB*ibin);
	 	 
        int kjj = kbin + NB*(jbin + NB*jbin);
        int jkj = jbin + NB*(kbin + NB*jbin);
        int jjk = jbin + NB*(jbin + NB*kbin);
        int kkj = kbin + NB*(kbin + NB*jbin);
        int kjk = kbin + NB*(jbin + NB*kbin);
        int jkk = jbin + NB*(kbin + NB*kbin);

        int kii = kbin + NB*(ibin + NB*ibin);
        int iki = ibin + NB*(kbin + NB*ibin);
        int iik = ibin + NB*(ibin + NB*kbin);
        int kki = kbin + NB*(kbin + NB*ibin);
        int kik = kbin + NB*(ibin + NB*kbin);
        int ikk = ibin + NB*(kbin + NB*kbin);

        int ijk = ibin + NB*(jbin + NB*kbin);
        int jki = jbin + NB*(kbin + NB*ibin);
        int kij = kbin + NB*(ibin + NB*jbin);
        int ikj = ibin + NB*(kbin + NB*jbin);
        int kji = kbin + NB*(jbin + NB*ibin);
        int jik = jbin + NB*(ibin + NB*kbin);

        for (int gridi =0; gridi<NG; gridi++)
	  for (int gridj =0; gridj<NG; gridj++)
	    for (int gridk =0; gridk<NG; gridk++){
	      int grid = gridk + 2*(NG/2 +1)*(gridj + NG*gridi);
	      Bg[iii]+=arr[ibin][grid]*arr[ibin][grid]*arr[ibin][grid]*delR; 
	      Bg[ijj]+=arr[ibin][grid]*arr[jbin][grid]*arr[jbin][grid]*delR;
	      Bg[jij]=Bg[jji]=Bg[ijj];
	      Bg[iij]+=arr[ibin][grid]*arr[ibin][grid]*arr[jbin][grid]*delR;
              Bg[iji]=Bg[jii]= Bg[iij];
	      Bg[kjj]+=arr[kbin][grid]*arr[jbin][grid]*arr[jbin][grid]*delR;
	      Bg[jkj]=Bg[jjk]=Bg[kjj];
              Bg[kii]+=arr[kbin][grid]*arr[ibin][grid]*arr[ibin][grid]*delR;
              Bg[iki]=Bg[iik]=Bg[kii];
	      Bg[kkj]+=arr[kbin][grid]*arr[kbin][grid]*arr[jbin][grid]*delR;
	      Bg[kjk]=Bg[jkk]= Bg[kkj];
              Bg[kki]+=arr[kbin][grid]*arr[kbin][grid]*arr[ibin][grid]*delR;
              Bg[kik]=Bg[ikk]= Bg[kki]; 
              Bg[ijk]+=arr[ibin][grid]*arr[jbin][grid]*arr[kbin][grid]*delR;
              Bg[jki]=Bg[kij]=Bg[ikj]=Bg[kji]=Bg[jik]=Bg[ijk];
        }
      }
    }  
  }  
  for (int i=0; i<NB; i++) free(arr[i]);
   free(arr);
   FILE *fout;
  if((fout = fopen(name,"w"))==0)  printf("cannot open pow output file\n");
  fprintf(stderr,"# Nyquist: %g\n",pi*(double)NG/BOX);
  double k1, k2, k3;
  for(int i=0;i<NB;i++) 
    for (int j=0;j<NB;j++)
      for (int k=0; k<NB; k++){
        int ind = k + NB*(j + i*NB);
        k1 = mink + binwidth * (i +0.5);
	k2 = mink + binwidth * (j +0.5);
        k3 = mink + binwidth * (k +0.5);
        fprintf(fout,"%d %g %g %g %g\n",i,k1,k2,k3,(Bg[ind]*pi/(k1*k2*k3*binwidth*binwidth*binwidth)));
      }
  fclose(fout);
  free(Bg); 
  fprintf(stderr,"bispectrum calculated for %s\n",name);

}

//function to find distance given redshift.
double calc_dp(double red) {
  //printf("red = %lf\n",red);
  double qsimp(double (*func)(double), double, double);
  double dpbit(double); 
  return qsimp(dpbit,0.,red);
}

double dpbit(double z) { 
  return 2997.92458/sqrt((1.+z)*(1.+z)*(1.+z)*Om_0 + (1-Om_0)); //3000.
} 
double Hubble(double z) {
  return H0*sqrt((1.+z)*(1.+z)*(1.+z)*Om_0 + (1-Om_0)); 
}

double Growth(double red) {
  double Ok = 1. - Om_0 - OL_0;
  double H2 = Om_0*(1.+red)*(1.+red)*(1.+red) + Ok*(1.+red)*(1.+red) + OL_0;
  double Om = Om_0*(1.+red)*(1.+red)*(1.+red)/H2;
  printf("Om = %lf\n",Om);
  double OL = OL_0/H2;
  printf("OL = %lf\n",OL);
  double pow_val = exp((0.550140)*log(Om));
  double D_z = (2.5*Om)/(((pow_val)-OL + (1. + 0.5*Om)*(1. + (OL/70.)))*(1.+red));
  return D_z;
}

double read_in_delta_field(double *delta, int NG, char *delta_file) {
  FILE *fp;
  int NGRID =0;
  int NVAL =0;
  double counter=0;
  double rho;
  if((fp=fopen(delta_file,"r"))==NULL) printf("delta file not opened\n");
  const int bsz=300; char buf[bsz];
  while((fgets(buf, bsz, fp))!=NULL) {
    double del;
    int i = floor(NGRID/(double)(NG*NG));
    int j = floor((NGRID-(i*NG*NG))/(double)NG);
    int k = NGRID - j*NG - i*NG*NG;
    sscanf(buf,"%lf\n",&del);
    NVAL = k + 2*(NG/2 + 1)*(j +NG*i);  
    delta[NVAL]=del;
    NGRID ++;
    counter = counter+del;
  }
  rho = counter/(double)NGRID;
  return rho;
}

