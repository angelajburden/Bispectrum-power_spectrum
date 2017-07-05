//check 64 bit
#if __GNUC__
#if __x86_64__ || __ppc64__
#define ENVIRONMENT64
#else
#define ENVIRONMENT32  
#endif
#endif

//if 32 bit, enable large file macros
#ifdef ENVIRONMENT32
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE
#define _FILE_OFFSET_BITS 64
#endif

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
  int npart[6];                    /*!< number of particles of each type in this file */
  double mass[6];                      /*!< mass of particles of each type. If 0, then the masses are explicitly
					 stored in the mass-block of the snapshot file, otherwise they are omitted */
  double time;                         /*!< time of snapshot file */
  double redshift;                     /*!< redshift of snapshot file */
  int flag_sfr;                    /*!< flags whether the simulation was including star formation */
  int flag_feedback;               /*!< flags whether feedback was included (obsolete) */
  unsigned int npartTotal[6];              /*!< total number of particles of each type in this snapshot. This can be
					 different from npart if one is dealing with a multi-file snapshot. */
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
//int my_snprintf(char *buffer,size_t len,const char *format, ...);
void my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
void * my_malloc(size_t size,int N);
FILE * my_fopen(const char *fname,const char *mode);
void print_time(time_t t0,time_t t1,const char *s);
struct basic_gal {
  double cp[3];
  double vel[3];
  double id;
};
//functions specific to gadget snapshots
struct io_header get_gadget_header(const char *fname);
int get_Numpart(struct io_header *header);


//Actual utilities being called. 
/*void my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nwritten;
  nwritten = fwrite(ptr, size, nmemb, stream);
  if(nwritten != nmemb)
    {
      fprintf(stderr,"I/O error (fwrite) has occured.\n");
      fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n",nmemb,nwritten);
      exit(EXIT_FAILURE);
    }
}


void my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t nread;
  nread = fread(ptr, size, nmemb, stream);
  if(nread != nmemb)
    {
      fprintf(stderr,"I/O error (fread) has occured.\n");
      fprintf(stderr,"Instead of reading nmemb=%zu, I got nread = %zu ..exiting\n",nmemb,nread);
      exit(EXIT_FAILURE);
    }
}


void * my_malloc(size_t size,int N)
{
  size_t bytes = size*(size_t)N;
  void *x = malloc(bytes);
  if(x==NULL)
    {
      fprintf(stderr,"ERROR: Failed to allocate memory for %d elements of size %zu bytes\n exiting\n",N,size);
      exit(EXIT_FAILURE);
    }

  return x;
}
*/
FILE * my_fopen(const char *fname,const char *mode)
{
  FILE *fp=NULL;
  fp = fopen(fname,mode);
  if(fp==NULL)
    {
      fprintf(stderr,"Could not open file `%s'\n",fname);
      exit(EXIT_FAILURE);
    }
  
  return fp;
}

struct io_header get_gadget_header(const char *fname)
{
  FILE *fp=NULL;
  char buf[MAXLEN], buf1[MAXLEN];
  int dummy;
  size_t one=1;
  struct io_header header;
  //my_snprintf(buf,MAXLEN, "%s.%d", fname, 0);
  //my_snprintf(buf1,MAXLEN, "%s", fname);
  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);
  if(sizeof(struct io_header) != 256)
    {
      fprintf(stderr,"ERROR: Gadget header is not %zu bytes and not *exactly* 256 bytes..exiting\n",sizeof(struct io_header));
      exit(EXIT_FAILURE);
    }

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

  //// Don't really care which file actually succeeded (as long as one, buf or buf1, is present)
  //  my_fread(&dummy, sizeof(dummy), one, fp);
  //my_fread(&header, sizeof(header), one, fp);
  //my_fread(&dummy, sizeof(dummy), one, fp);
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


/*void print_time(time_t t0,time_t t1,const char *s)
{
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


// A real wrapper to snprintf that will exit() if the allocated buffer length 
// was not sufficient. Usage is the same as snprintf 
int my_snprintf(char *buffer,size_t len,const char *format, ...)
{
  va_list args;
  int nwritten=0;
 
  va_start(args,format);
  nwritten=vsnprintf(buffer, len, format, args );
  va_end(args);  
  if ((size_t) nwritten > len || nwritten < 0)
    {
      fprintf(stderr,"ERROR: printing to string failed (wrote %d characters while only %zu characters were allocated)\n",nwritten,len);
      fprintf(stderr,"Increase MAXLEN ..exiting\n");
      exit(EXIT_FAILURE);
    }
  return nwritten;
  }*/


struct basic_gal *gal;
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int main(int argc, char **argv){

  char inpath[MAXLEN],outpath[MAXLEN],snapshot_base[MAXLEN],snapshot_name[MAXLEN],outfname[MAXLEN];

  int nfiles,ifile;
  struct io_header header,header1;
  int NumPart;
  FILE *fd=NULL;//,*fdmass=NULL;
  double BoxSize;
  int dummy;
  float pos[3];
  size_t one=1;
  //  time_t t_codestart,t_codeend;

  //  t_codestart=time(NULL);
  if (argc !=3) {
   
      fprintf(stderr,"\n\n wrong argument(s).  Specify: `%s' \n\n",argv[0]);
      fprintf(stderr,"<inpath>      (input path without trailing slash)\n");
      fprintf(stderr,"<basename>    (snapshot files)\n");
      fprintf(stderr,"<outpath>     (output path)\n");
      fprintf(stderr,"\n exiting ..\n");
      exit(EXIT_FAILURE);
  }
  else {
   
      //copy the input directory (containing the snapshots)
    //  my_snprintf(inpath,MAXLEN,"%s",argv[1]);
      sprintf(inpath,"%s",argv[1]);
      //copy the basename of the snapshots
      //my_snprintf(snapshot_base,MAXLEN,"%s",argv[2]);
      sprintf(snapshot_base,"%s",argv[2]);
      fprintf(stderr,"Running `%s' with the parameters:\n\n",argv[0]);
      fprintf(stderr,"\t\t <inpath>   = %s \n",inpath);
      fprintf(stderr,"\t\t <basename> = %s \n",snapshot_base);
      fprintf(stderr,"\n");
  }

  // my_snprintf(snapshot_name,MAXLEN,"%s/%s",inpath, snapshot_base);
  sprintf(snapshot_name,"%s/%s",inpath, snapshot_base);
  header  =  get_gadget_header(snapshot_name);
  nfiles  =  header.num_files;
  NumPart = get_Numpart(&header);
  BoxSize = header.BoxSize;

  //allocate memory for structure gal which has the basic_gal format
  if(!(gal = (struct basic_gal*)malloc(NumPart*sizeof(struct basic_gal))-1))
    printf("memory allocation problem for galaxies\n");
 
  int ii=0;
  ifile =0;
  int index=0;
  float xmax=0.0;
  float ymax =0.0;
  float zmax =0.0;
  //  my_snprintf(snapshot_name,MAXLEN,"%s/%s",inpath, snapshot_base);
  sprintf(snapshot_name,"%s/%s",inpath, snapshot_base);
  fd = fopen(snapshot_name,"r");
  //  fd = my_fopen(snapshot_name,"r");
  fprintf(stderr,"Reading file `%s' ...",snapshot_name);
  //  my_fread(&dummy, sizeof(dummy), one, fd);
  //my_fread(&header1, sizeof(header1), one, fd);
  //my_fread(&dummy, sizeof(dummy), one, fd);
  //my_fread(&dummy, sizeof(dummy), one, fd);//set it ready to read in the first of the particle positions
  fread(&dummy, sizeof(dummy), one, fd);
  fread(&header1, sizeof(header1), one, fd);
  fread(&dummy, sizeof(dummy), one, fd);
  fread(&dummy, sizeof(dummy), one, fd);//set it ready to read in the first of the particle positions        
  for(int k=0; k<6; k++){
    for(int npart=0; npart<header1.npart[k]; npart++){
      //read in the position from the snapshot
      //my_fread(&pos[0],sizeof(float),3*one,fd);
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
  printf("xmax = %lf, ymax = %lf, zmax = %lf\n", xmax, ymax, zmax);
     
  //my_snprintf(outfname,MAXLEN,"particle_positions_gaget.txt");
  //  sprintf(outfname,"particle_positions_gaget.txt");
  //fd = my_fopen(outfname,"w");
  //fd = fopen(outfname,"w");
  // for (int i=0; i<NumPart; i++)
  // if (gal[i].cp[2]<152.0 && gal[i].cp[2]>=150.0) fprintf(fd, "%lf %lf %lf\n",gal[i].cp[0], gal[i].cp[1], gal[i].cp[2]);
    
  //t_codeend = time(NULL);
  // print_time(t_codestart,t_codeend,"read in time");
  return EXIT_SUCCESS;
}




