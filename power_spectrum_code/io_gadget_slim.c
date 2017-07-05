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


//functions specific to gadget snapshots
//struct io_header get_gadget_header(const char *fname);
//int get_Numpart(struct io_header *header);

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

//struct basic_gal *gal;  //declare a structure called gal

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//start main program
void read_in_binary_gadget(struct basic_gal *gal, char *inpath, char *snapshot_base, int snapshot, char *outfname, int *NGAL){
//int main(int argc, char **argv){

//char inpath[MAXLEN],snapshot_base[MAXLEN],snapshot_name[MAXLEN],outfname[MAXLEN];
  char snapshot_name[MAXLEN];
  int nfiles;
  struct io_header header,header1;
  int NumPart;
  FILE *fd=NULL;
  double BoxSize;
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
  float xmax=0.0;
  float ymax =0.0;
  float zmax =0.0;
  int ii=0;
  while(ifile < nfiles){
    if (nfiles == 1) {
      sprintf(snapshot_name,"%s/%s_%03d",inpath,snapshot_base,snapshot);
    }
    else if (nfiles > 1) {
      sprintf(snapshot_name,"%s/%s_%03d.%d",inpath,snapshot_base,snapshot,ifile);
    }

    fd = fopen(snapshot_name,"r");
    fprintf(stderr,"Reading file `%s' ...",snapshot_name);

    fread(&dummy, sizeof(dummy), one, fd);
    fread(&header1, sizeof(header1), one, fd);
    fread(&dummy, sizeof(dummy), one, fd);
    fread(&dummy, sizeof(dummy), one, fd);//set it ready to read in the first of the particle positions        

    for(int k=0; k<6; k++){
      for(int npart=0; npart<header1.npart[k]; npart++){
        //read in the position from the snapshot
        fread(&pos[0],sizeof(float),3*one,fd);
        //allocate to structure
        gal[ii].cp[0] = pos[0];
        gal[ii].cp[1] = pos[1];
        gal[ii].cp[2] = pos[2];
        if (pos[0]>xmax)xmax = pos[0];
        if (pos[1]>ymax)ymax = pos[1];
        if (pos[2]>zmax)zmax = pos[2];
        ii++;
        index++;
       }  
    }
    fprintf(stderr,"\nReading file `%s' ...done\n",snapshot_name);
    fclose(fd);
    ifile++;
  }
  printf("xmax = %lf, ymax = %lf, zmax = %lf\n", xmax, ymax, zmax);
  *NGAL = NumPart;
}
   
  // sprintf(outfname,"particle_positions_gaget.txt");
  //fd = fopen(outfname,"w");
  /*  for (int i=0; i<NumPart; i++)
//    if (gal[i].cp[2]<152.0 && gal[i].cp[2]>=150.0) fprintf(fd, "%lf %lf %lf\n",gal[i].cp[0], gal[i].cp[1], gal[i].cp[2]);

}*/




