//module for power spectrum functions
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
// function to calculate power spectrum of cubic grid and output to file

void fun_calc_power(double *dd, int NG, int NGAL, double BOX, char *name) {

  rfftwnd_plan dp1;
  dp1 = rfftw3d_create_plan(NG,NG,NG,FFTW_REAL_TO_COMPLEX,FFTW_MEASURE | FFTW_IN_PLACE);
  //  rfftwnd_one_real_to_complex(dp1,(fftw_real*)dd,NULL);
  // rfftwnd_destroy_plan(dp1);
  /*
  printf("fourier transformed\n");
  // Set up bins for P(k) data 
  const int    NB   = 400;                // Number of bins
  const double mink = 0.0;                // minimum k to be binned
  const double maxk = 0.7;//2.*pi*0.25;        // maximum k to be binned
  double binwidth=(maxk-mink)/(double)NB; // bin width
  double *Pg=dvector(1,NB);               // binned power spectrum
  double *Nk=dvector(1,NB);  
  int    *Pgfill=ivector(1,NB);            
  for(int i=0;i<NB;i++) Pg[i]=Nk[i]=0.0;
  for(int i=0;i<NB; i++) Pgfill[i]=0;

  // Fourier transform density field
  printf("Now Fourier transforming overdensity field\n");

  double Nyquist = pi*(double)NG/2*BOX;    
  // Work out |k| for each frequency component
  double fx, fy, fz;
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
	int bin = (int)( (double)(fktot-mink)/binwidth + 0.5);
	    
	if(bin>=1 && bin<NB  && Pgfill[bin]==0) {
	      
	  double sinc_x = 1.0, sinc_y=1.0, sinc_z=1.0;
	  double ax = pi*fx*(BOX/(double)NG);
	  double ay = pi*fy*(BOX/(double)NG);
	  double az = pi*fz*(BOX/(double)NG);
	  if(fx!=0.0) sinc_x = sin(ax)/ax;
	  if(fy!=0.0) sinc_y = sin(ay)/ay;
	  if(fz!=0.0) sinc_z = sin(az)/az;
	  double grid_cor = 1.0/(sinc_x*sinc_y*sinc_z);
	  
	  double dkr = dd[(2*k  )+2*(NG/2+1)*(j+NG*i)];
	  double dki = dd[(2*k+1)+2*(NG/2+1)*(j+NG*i)]; 
	  //for data use first instance, for periodic box, second.
	  //Pg[bin] += ((dkr*dkr+dki*dki - (tran_weightsq + tgal_weightsq))*grid_cor*grid_cor);
          Pg[bin] += ((dkr*dkr+dki*dki)- 1./(NGAL));//*grid_cor*grid_cor);
          if (bin ==0) Pg[bin]= 0.0;
	  Nk[bin]++;	      
	}	
      }
    }
  }
  printf("normalise\n");
  for(int i=0;i<NB;i++) {
    if(Nk[i]>0.0 && Pgfill[i]==0 ) {
      Pg[i] /= (Nk[i]);
      Pg[i] /= (double)(BOX*BOX*BOX);
      Pgfill[i] =1; 
    } 
  }

  // output power spectrum values*/
  const int NGMAX = NG*NG*2*(NG/2+1);
  double *dr = (double*)malloc(sizeof(double)*NGMAX);
  double nbar_grid = (double)NGAL/(double)(NG*NG*NG);
  for(int igg=0;igg<NGMAX;igg++) dr[igg]=dd[igg]*nbar_grid;
  printf("nbar = %lf\n",nbar_grid);
  // Set up bins for P(k) data 
  const int    NB   = 400;                // Number of bins
  const double mink = 0.0;                // minimum k to be binned
  const double maxk = 2.*pi*0.256;        // maximum k to be binned
  double binwidth=(maxk-mink)/(double)NB; // bin width
  double *Pg=dvector(1,NB);               // binned power spectrum
  double *Nk=dvector(1,NB);               // number of modes in each bin
  for(int i=1;i<=NB;i++) Pg[i]=Nk[i]=0.0;
  double *GaussFT = (double*)malloc(sizeof(double)*NGMAX);
  // Fourier transform density field
  printf("Now Fourier transforming overdensity field\n");
      
  // fftw_execute(std_p_r2c);
  rfftwnd_one_real_to_complex(dp1,(fftw_real*)dr,NULL);
      
  // Work out |k| for each frequency component
  double fx, fy, fz;
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
	//int q = k + 2*(NG/2 +1)*(j +(NG)*i);   
	// length of k vector and bin
	double fktot=2.*pi*sqrt(fx*fx+fy*fy+fz*fz);
	int bin = (int)( (double)(fktot-mink)/binwidth );
	//GaussFT[q] = exp(-0.5*fktot*fktot*RG*RG);   
	if(bin>=1 && bin<=NB) {
	        
	  // set up correction for gridding - in effect we're
	  // convolving the density field with a top-hat function in
	  // each direction, so we're multiplying each Fourier mode by
	  // a sinc function. To correct this, we therefore divide by
	  // the sinc functions.
	  double sinc_x = 1.0, sinc_y=1.0, sinc_z=1.0;
	  double ax = pi*fx*(BOX/(double)NG);
	  double ay = pi*fy*(BOX/(double)NG);
	  double az = pi*fz*(BOX/(double)NG);
	  if(fx!=0.0) sinc_x = sin(ax)/ax;
	  if(fy!=0.0) sinc_y = sin(ay)/ay;
	  if(fz!=0.0) sinc_z = sin(az)/az;
	  double grid_cor = 1.0/(sinc_x*sinc_y*sinc_z);
	    
	  double dkr = dr[(2*k  )+2*(NG/2+1)*(j+NG*i)];
	  double dki = dr[(2*k+1)+2*(NG/2+1)*(j+NG*i)]; 
	        
	  double fac = ((k==0) || (k==NG/2)) ? 1 : 2;
          Pg[bin] += fac*(dkr*dkr+dki*dki);
          Nk[bin]+= fac;       
	}
      }
    }
  }
  //  printf ("NGAL_MAX = %d\n", NGAL_MAX);
  printf("NGAL = %d\n",NGAL);
  for(int i=1;i<=NB;i++) {
    if(Nk[i]>0.0) Pg[i] /= (Nk[i]); else Pg[i]=0.0;
    Pg[i] *= (BOX*BOX*BOX)/((double)NGAL*(double)NGAL);
  }

  // output power spectrum values
  FILE *fout;
  if((fout = fopen(name,"w"))==0) 
    printf("cannot open pow output file\n");
  //    printf("# Nyquist: %g\n",pi*(double)NG/BOX);
  double *kp=dvector(1,NB);
  for(int i=0;i<NB;i++) {
    kp[i] = mink+((double)i+0.5)*binwidth;
    fprintf(fout,"%d %g %g\n",i,kp[i],Pg[i]);
  }
  fclose(fout);
  //plot power-spectrum
  //  xy_plot(kp, Pg, 1, NB);
  free_dvector(Pg,1,NB); 
  free_dvector(Nk,1,NB);
  printf("power spectrum calculated for %s\n",name);

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


/*void calculate_fkp(double min_red, double max_red, double SURVEY_SIZE, int NGAL, struct basic_gal *gal, float nbar_slice [], int NRB, float *redslice, float *red_vol, double *tgal_weights, double *tgal_weightsq, double *tgal_weight_nbar, int flag) {
 

  double FULL_AREA = (4*pi);
  *tgal_weights =0.0;
  *tgal_weightsq =0.0;
  *tgal_weight_nbar =0.0;

    //calculate density of redshift bin.
    for (int i =1; i<NGAL; i++) {
      int bin = (int)(NRB*(gal[i].z-min_red)/(max_red -min_red));
      if (bin<0) printf("bin = %d, gal[%d].z = %lf\n",bin,i,gal[i].z);
      redslice[bin] +=gal[i].weight;
    }

    for (int i=0; i<NRB; i++) {
      double red_shell_min = ((double)i*(max_red-min_red)/NRB) + min_red;
      double red_shell_max = ((double)(i+1)*(max_red-min_red)/NRB) + min_red;
      double d_min = calc_dp(red_shell_min);
      double d_max = calc_dp(red_shell_max);
      red_vol[i] = ((4.0/3.0)*(pi)*((d_max*d_max*d_max)-(d_min*d_min*d_min)))*(SURVEY_SIZE/FULL_AREA);
      nbar_slice[i] =(double)redslice[i]/red_vol[i];
    }

  
  if (flag ==2 || flag ==3 ||flag ==4) {
    for (int i=1; i<NGAL; i++) {

      int bin = (int)(NRB*(gal[i].z-min_red)/(max_red -min_red));
      if (bin<0) printf("bin = %d, gal[%d].z = %lf\n",bin,i,gal[i].z);

      gal[i].nbar = nbar_slice[bin];
      gal[i].fkp = 1.0/((1.0+20000.0*nbar_slice[bin]));
  }
}
}

 
void bin_gals_NGP(struct basic_gal *gal, double  min_x, double min_y, double min_z, int  NG, float dr[], float dr_power[], int NGAL, double L, int flag, double alpha) {
  //flag =0 for data/mocks
  //flag =1 for all rands.

  int NGK = NG*NG*2*(NG/2+1);
  printf("flag = %d\n",flag);

  for (int i=1; i<NGAL; i++) {
    int xp = (int)(((gal[i].cp[0]-min_x)/L) * NG);
    int yp = (int)(((gal[i].cp[1]-min_y)/L) * NG);
    int zp = (int)(((gal[i].cp[2]-min_z)/L) * NG);
   
    int ind = zp + 2*(NG/2 +1)*(yp + NG*xp); 


    if (ind >NGK || ind<0) {
      printf("xp= %d, yp =%d, zp=%d\n",xp,yp,zp);
      printf("i =%d, galx = %lf, dist =%lf, z= %lf\n",i,gal[i].cp[0],gal[i].dist,gal[i].z);
    }

    if (flag ==0) {
      dr_power[ind] += (gal[i].fkp*gal[i].weight); 
      dr[ind] +=gal[i].weight;
    }

    if (flag ==1) {
      dr_power[ind] -= (gal[i].fkp) *alpha; 
      dr[ind] +=gal[i].fkp;
    }
    if (flag ==2) {
      dr_power[ind] -= (gal[i].weight *alpha); 
      dr[ind] +=gal[i].weight;
    }
  }  
}

void smooth_field(double RG, float dr[],float dr_uni[], int NG, double L, int flag) {

  double GaussFT=0.0;
  for (int i=0;i<NG;i++) 
    for (int l=0;l<NG;l++)
      for (int m=0;m<=(NG/2);m++) {
        double fx, fy, fz;
        if (i<=NG/2 && i>-NG/2) fx =(double)i*2*pi/L;
        if (i>NG/2) fx = (double)(i-NG)*2*pi/L;
       
        if (l<=NG/2 && l>-NG/2) fy =(double)l*2*pi/L;
        if (l>NG/2) fy = (double)(l-NG)*2*pi/L;

        if (m<=NG/2 && m>-NG/2) fz =(double)m*2*pi/L;

        double ks = (fx*fx + fy*fy + fz*fz);  
        GaussFT = exp(-0.5*ks*RG*RG);

        int RE = (2*m) +2*(NG/2+1)*(l+NG*i);
        int IM = (2*m+1)+2*(NG/2+1)*(l+NG*i);

        if (flag ==0) {      
          dr[RE] *=GaussFT;
          dr[IM] *=GaussFT;
	  dr_uni[RE] *=GaussFT;
          dr_uni[IM] *=GaussFT;
        }

  }
}


void read_in_binary_mock(const char *data, struct basic_gal *gal, long int *NGAL, long int NGAL_MAX)  {

   typedef struct
  {
    unsigned long long int id;
    double M;
    double x[3], v[3], ID;
    int N;
  } cat;
  cat catdata;
  unsigned long long int id,n;
  int nproc,iproc,ngood,igood,i,dummy;
  char runname[100],filename[100],redshift[100];
  FILE *fd;


  sprintf(filename,"%s",data);
  printf("Opening file %s\n",filename);

  long int j=0;
  fd=fopen(filename,"r");
  if (fd==0x0)
    {
      printf("Error: catalog file %s not found\n",filename);
    }

  n=0;
 
  fread(&dummy,sizeof(int),1,fd);
  fread(&nproc,sizeof(int),1,fd);
  fread(&dummy,sizeof(int),1,fd);
  printf("The file has been written by %d tasks\n",nproc);
  for (iproc=0; iproc<nproc; iproc++)
    {
      //      fread(&dummy,sizeof(int),1,fd);
      fread(&ngood,sizeof(int),1,fd);
      // fread(&dummy,sizeof(int),1,fd);
      *NGAL+=ngood;
      n+=ngood;

      for (igood=0; igood<ngood; igood++)
	{
	  fread(&dummy,sizeof(int),1,fd);
	  fread(&catdata,dummy,1,fd);
	  fread(&dummy,sizeof(int),1,fd);
	  gal[j].cp[0]=catdata.x[0];
	  gal[j].cp[1]=catdata.x[1];
	  gal[j].cp[2]=catdata.x[2];
	  gal[j].vel[0]=catdata.v[0];
	  gal[j].vel[1]=catdata.v[1];
	  gal[j].vel[2]=catdata.v[2];
	  j++;
	}
    }
  fclose(fd);

  printf("I found %Ld halos in the file, i=%ld\n",n, j);
  */
//}
