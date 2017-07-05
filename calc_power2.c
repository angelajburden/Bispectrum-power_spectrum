#include "header.h"
#include <rfftw.h>

/*struct basic_gal {
  double cp[3];
};
struct basic_gal *gal;*/
const char *data = ("new_positionsM.dat");///home/angela/fof0p156_Oriana_ir1001gauss_z0p000.sorted");//

void calc_power(double*,int,int,double);
const int NGAL_MAX = 1024.*1024.*1024;
const int NG=512;
const double BOX = 1500.0;
int NGAL = 0;
// parameters for data
const int NGMAX = NG*NG*2*(NG/2+1);
const double delx = BOX/(double)NG; //width of each grid Mpc/h.
//const double RG= 10.0/delx;
int main() {

  // parameters for HV simulation and this analysis
 
  double *dr = (double*)malloc(sizeof(double)*NGMAX);
  /*allocate memory*/
  // if(!(gal = (struct basic_gal*)malloc(NGAL_MAX*sizeof(struct basic_gal))-1)) 
  // printf("memory allocation problem for galaxies");

  for(int igg=0;igg<NGMAX;igg++) dr[igg]=0.0;

  /*read in data*/
  /* FILE *fp;
  if((fp=fopen(data,"r"))==NULL) printf("Data file not opened");
  const int bsz=400; char buf[bsz];

  while(fgets(buf,bsz,fp)) {
    double x,y,z; 
    sscanf(buf,"%lf,%lf,%lf",&x,&y,&z);//%*d %*d %*f %lf %lf %lf",&x,&y,&z);//
    if(++NGAL>NGAL_MAX){ NGAL--; break;}
    gal[NGAL].cp[0]= x;
    gal[NGAL].cp[1]= y;
    gal[NGAL].cp[2]= z;
  }
  fclose(fp);*/


  for (int i=1; i<=NGAL_MAX; i++) {

    int xp = (int)((gal[i].cp[0]/BOX) * NG);
    int yp = (int)((gal[i].cp[1]/BOX) * NG);
    int zp = (int)((gal[i].cp[2]/BOX) * NG);

    //offset in pixel
    double dx = (double)NG*(gal[i].cp[0]/BOX) - (double)xp;
    double dy = (double)NG*(gal[i].cp[1]/BOX) - (double)yp;
    double dz = (double)NG*(gal[i].cp[2]/BOX) - (double)zp;

    while(xp>=NG) xp-=NG; // fold back onto NGRID
    while(yp>=NG) yp-=NG; // fold back onto NGRID
    while(zp>=NG) zp-=NG; // fold back onto NGRID; 

    while(xp<0) xp+=NG; // fold back onto NGRID
    while(yp<0) yp+=NG; // fold back onto NGRID
    while(zp<0) zp+=NG; // fold back onto NGRID;

    double tx=0.0, ty=0.0, tz=0.0;
    for (int i=0; i<=1; i++)
      for (int j=0; j<=1; j++) 
	for (int k=0; k<=1; k++) {
	  int a=xp+i, b=yp+j, c=zp+k;
	  while (a>=NG) a-=NG;
	  while (b>=NG) b-=NG;
	  while (c>=NG) c-=NG;
	  int ind = c+2*(NG/2+1)*(b+NG*a);
	  if (ind <0 || ind > NGMAX) printf("the grid data is out of range\n");
	  if (i==0) tx=1.0-dx; else tx=dx;
	  if (j==0) ty=1.0-dy; else ty=dy;
	  if (k==0) tz=1.0-dz; else tz=dz;

	  dr[ind] += tx*ty*tz;        
	}
  }


  // convert from galaxy distribution to overdensity
  double nbar_grid = (double)NGAL/(double)(NG*NG*NG);
  for(int igg=0;igg<NGMAX;igg++) dr[igg] = (dr[igg]-nbar_grid) / nbar_grid;

  // estimate power

  calc_power(dr,NG,NGAL,BOX);
   
  return 1;
}


// function to calculate power spectrum of cubic grid and output to file
void calc_power(double *dd, int NG, int NGAL, double BOX) {

  rfftwnd_plan dp_r2c;
  dp_r2c = rfftw3d_create_plan(NG,NG,NG,(fftw_direction)-1,
			       FFTW_MEASURE | FFTW_IN_PLACE);

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
  rfftwnd_one_real_to_complex(dp_r2c,(fftw_real*)dr,NULL);
      
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
	int bin = 1+(int)( (double)(fktot-mink)/binwidth );
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
	        
	  Pg[bin] += ((dkr*dkr+dki*dki)-(NGAL_MAX))*grid_cor*grid_cor;
	  Nk[bin]++;
	        
	}
      }
    }
  }
  printf ("NGAL_MAX = %d\n", NGAL_MAX);
  printf("NGAL = %d\n",NGAL);
  for(int i=1;i<=NB;i++) {
    if(Nk[i]>0.0) Pg[i] /= (Nk[i]); else Pg[i]=0.0;
    Pg[i] *= (BOX*BOX*BOX)/((double)NGAL*(double)NGAL);
  }

  // output power spectrum values
  FILE *fout;
  if((fout = fopen("Po
