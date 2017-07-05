#include "header.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdint.h>
#include<inttypes.h>
#include<math.h>
#include<string.h>
#include<limits.h>
#include<time.h>
#include<stdarg.h>

const size_t MAXLEN=1000;//number of characters in a string (pathnames, filenames etc)

struct io_header
{
  int npart[6];                    //< number of particles of each type in this file */
  double mass[6];                      //!< mass of particles of each type. If 0, then the masses are explicitly
  double time;                         //< time of snapshot file */
  double redshift;                     //< redshift of snapshot file */
  int flag_sfr;                    //!< flags whether the simulation was including star formation */
  int flag_feedback;               //!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];              //< total number of particles of each type in this snapshot. This can be
  int flag_cooling;                /*!< flags whether cooling was included  */
  int num_files;                   /*!< number of files in multi-file snapshot */
  double BoxSize;                      /*!< box-size of simulation in case periodic boundaries were used */
  double Omega0;                       /*!< matter density in units of critical density */
  double OmegaLambda;                  /*!< cosmological constant parameter */
  double HubbleParam;                  /*!< Hubble parameter in units of 100 km/sec/Mpc */
  int flag_stellarage;             /*!< flags whether the file contains formation times of star particles */
  int flag_metals;                 /*!< flags whether the file contains metallicity values for gas and star particles */
  unsigned int npartTotalHighWord[6];      /*!< High word of the total number of particles of each type */
  int  flag_entropy_instead_u;     /*!< flags that IC-file contains entropy instead of u */
  char fill[60];                       /*!< fills to 256 Bytes */
};

struct io_header get_gadget_header(const char *fname)
{
  FILE *fp=NULL;
  char buf[MAXLEN], buf1[MAXLEN];
  int dummy;
  size_t one=1;
  struct io_header header;
  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);
  fp = fopen(buf,"r");
  if(fp == NULL)
    {
      fp = fopen(buf1,"r");
      if(fp == NULL)
	{
	  fprintf(stderr,"ERROR: Could not find snapshot file.\n neither as `%s' nor as `%s'\n",buf,buf1);
	  fprintf(stderr,"exiting..\n");
	  exit(EXIT_FAILURE);
	}
    }

  fread(&dummy, sizeof(dummy), one, fp);
  fread(&header, sizeof(header), one, fp);
  fread(&dummy, sizeof(dummy), one, fp); 
  fclose(fp);
  return header;

}

int get_Numpart(struct io_header *header)
{
  int NumPart=0;

  if(header->num_files <= 1)
    for(int i = 0; i < 6; i++)
      header->npartTotal[i] = header->npart[i];
  
  for(int i = 0; i < 6; i++) {
      NumPart += header->npartTotal[i];
      NumPart += (((int) header->npartTotalHighWord[i]) << 32);
  }
  return NumPart;
}
//function to time the code
void print_time(time_t t0,time_t t1,const char *s){                                                                      
  double timediff = difftime(t1,t0);
  double ratios[] = {24*3600.0,  3600.0,  60.0,  1};                                                                    
  char units[4][10]  = {"days", "hrs" , "mins", "secs"};                                                          
  int which = 0;
  double timeleft = timediff;                                                                                            
  double time_to_print;
  fprintf(stderr,"Time taken to execute '%s'  = ",s);
  if(timediff < ratios[2])
    fprintf(stderr,"%5d secs",(int)timediff);
  else
    while (which < 4)
      {
        time_to_print = floor(timeleft/ratios[which]);
        if (time_to_print > 1)
          {
            timeleft -= (time_to_print*ratios[which]); 
            fprintf(stderr,"%5d %s",(int)time_to_print,units[which]);
          }
        which++;
      }
  fprintf(stderr,"\n");
}
                       
void read_in_binary_gadget(double *delta, int NG, double L, char *inpath, char *snapshot_base, int snapshot, char *outfname, long int *NGAL){

  char snapshot_name[MAXLEN];
  int nfiles;
  struct io_header header,header1;
  int NumPart;
  FILE *fd=NULL;
  double BoxSize;
  double OM, OL, H0;
  int dummy;
  float pos[3];
  size_t one=1;
  nfiles = snapshot;
  int  ifile =0;
  int  index=0;
  // my_snprintf(snapshot_name,MAXLEN,"%s/%s",inpath, snapshot_base);
  sprintf(snapshot_name,"%s/%s_%03d",inpath, snapshot_base, snapshot);
  header  =  get_gadget_header(snapshot_name);
  nfiles  =  header.num_files;
  NumPart = get_Numpart(&header);
  BoxSize = header.BoxSize;
  OM =  header.Omega0;
  OL = header.OmegaLambda;
  H0 = header.HubbleParam;
  printf("cosmo params, OM =%lf, OL= %lf, H0 = %lf\n", OM, OL, H0); 
  float xmax=0.0;
  float ymax =0.0;
  float zmax =0.0;
  // int ii=0;
  while(ifile < nfiles){
    if (nfiles == 1) {
      sprintf(snapshot_name,"%s/%s_%03d",inpath,snapshot_base,snapshot);
    }
    else if (nfiles > 1) {
      sprintf(snapshot_name,"%s/%s_%03d.%d",inpath,snapshot_base,snapshot,ifile);
    }

    fd = fopen(snapshot_name,"r");
    fprintf(stderr,"Reading file `%s' ...",snapshot_name);
    const int NGMAX = NG*NG*2*(NG/2+1);
    fread(&dummy, sizeof(dummy), one, fd);
    fread(&header1, sizeof(header1), one, fd);
    fread(&dummy, sizeof(dummy), one, fd);
    fread(&dummy, sizeof(dummy), one, fd);//set it ready to read in the first of the particle positions        

    for(int k=0; k<6; k++){
      for(int npart=0; npart<header1.npart[k]; npart++){
        //read in the position from the snapshot
        fread(&pos[0],sizeof(float),3*one,fd);
        //allocate to structure
	/* gal[ii].cp[0] = pos[0];
        gal[ii].cp[1] = pos[1];
        gal[ii].cp[2] = pos[2];*/
        if (pos[0]>xmax)xmax = pos[0];
        if (pos[1]>ymax)ymax = pos[1];
        if (pos[2]>zmax)zmax = pos[2];
	// ii++;
        //index++;*/
	int xp = (int)(((pos[0])/(double)L) *(double) NG);
	int yp = (int)(((pos[1])/(double)L) *(double) NG);
	int zp = (int)(((pos[2])/(double)L) *(double) NG);
        double dx = (double)NG*(pos[0]/L) - (double)xp;
	double dy = (double)NG*(pos[1]/L) - (double)yp;
	double dz = (double)NG*(pos[2]/L) - (double)zp;
	while(xp>=NG) xp-=NG; // fold back onto NGRID
	while(yp>=NG) yp-=NG; // fold back onto NGRID
	while(zp>=NG) zp-=NG; // fold back onto NGRID

	while(xp<0) xp+=NG; // fold back onto NGRID
	while(yp<0) yp+=NG; // fold back onto NGRID
	while(zp<0) zp+=NG; // fold back onto NGRID
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

	      delta[ind] += tx*ty*tz;
	 }
    	/*int ind = zp + 2*(NG/2 +1)*(yp + NG*xp);
	if (xp>NG || xp<0) printf("gal xp = %d, galx = %lf\n", xp, pos[0]);
	if (yp>NG || yp<0)printf("gal yp = %d, galy = %lf\n", yp, pos[1]);
	if (zp>NG || zp<0) printf("gal zp = %d\n", zp, pos[2]);
	delta[ind] = delta[ind]+1.0;
        if (pos[0]>xmax)xmax = pos[0];                                                                     
        if (pos[1]>ymax)ymax = pos[1];                                                                                                     
        if (pos[2]>zmax)zmax = pos[2];*/                                                                                                 
        index++;
       }  
    }
    fprintf(stderr,"\nReading file `%s' ...done\n",snapshot_name);
    fclose(fd);
    ifile++;
  }
  fprintf(stderr,"xmax = %lf, ymax = %lf, zmax = %lf\n", xmax, ymax, zmax);
  *NGAL = NumPart;
}
   




