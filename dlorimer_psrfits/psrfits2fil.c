#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include "psrfits.h"
#include "header.h"

/*
 * This version of psrfits2fil is specifically for 20 converting
 * full stokes 20 m data on Cyborg
 */
void send_coords(double raj, double dej, double az, double za) ;
void send_double (char *name, double double_precision) ;
void send_float(char *name,float floating_point) ;
void send_int(char *name, int integer) ;
void send_long(char *name, long integer) ;
void send_string(char *string) ;
FILE *input, *output;
int swapout;

main (int argc, char **argv){
  unsigned char *data8;
  unsigned short *data16;
  int flip=0,i,j,k,l,x,status,first=1,startchan,endchan,idump=0,ndumps=1;
  int counter,bandpass=0,npersub;
  struct psrfits pf;
  double toff,fcent,tsamp,chbw;
  double rah,ram,ras,ded,dem,des,sgn;
  char filfile[1024],stem[1024], *pos;

  startchan=endchan=0;
//  startchan=695;
//  endchan=950;

  if (argc < 2) {
    printf("usage: psrfits2fil fitsfile (startchan) (endchan) (flip) (fcentMHz) (tsampus)\n");
    exit(0);
  }

  if (argc>2) {
    if (strstr(argv[2],"bandpass")!=NULL)  bandpass=1;
    startchan=atoi(argv[2]);
  }

  if (argc>3) {
    endchan=atoi(argv[3]);
  }

  if (argc>4) {
    if (strstr(argv[4],"flip")!=NULL)  flip=1;
  }

  fcent=0.0;
  if (argc>5) {
    fcent=atof(argv[5]);
  }

  tsamp=0.0;
  if (argc>6) {
    tsamp=1.0e-6*atof(argv[6]);
  }

  if ( (pos = strstr(argv[1],".fits")) ) {
    strncpy(pf.basefilename,argv[1],strlen(argv[1])-10);
    strcpy(pf.filename,argv[1]);
    strncpy(stem,argv[1],strlen(argv[1])-5);
  } else {
    puts("Your input file doesn't end in .fits");
    exit(0);
  }
  sscanf(&stem[strlen(stem) - 4], "%d", &pf.filenum);

  char basefilename[1024];
  strcpy(basefilename, pf.basefilename);
  char *filbasename = basename(basefilename);

  sprintf(filfile,"%s_%04d.fil", filbasename, pf.filenum);
  //sprintf(filfile,"%04d.fil",pf.filenum);
  printf("Output %s\n", filfile);
  output=fopen(filfile,"wb");
  if (psrfits_open(&pf) != 0) {
    fprintf(stderr,"error opening file %s\n",pf.filename);
    exit(0);
  }

  counter=0;
  while ( psrfits_read_subint(&pf)== 0) {
    if (first) {
      if (fcent == 0.0) fcent=pf.hdr.fctr;
	fprintf(stderr,"Center frequency %f MHz\n",fcent);
      chbw=pf.hdr.BW/pf.hdr.nchan;
	if (chbw>0) flip=1;
	//flip=0;
  if (flip==0) printf("No flipping of the channels!\n");
  if (flip==1) printf("Flipping channels!\n");
      /* broadcast header file */
      send_string("HEADER_START");
      send_string("rawdatafile");
      send_string(filbasename);
      send_string("source_name");
      send_string(pf.hdr.source);
      send_int("data_type",1);
      if (startchan==0) startchan=1;
      if (endchan==0) endchan=pf.hdr.nchan;
      send_int("nchans",endchan-startchan+1);
	fprintf(stderr,"Output number of channels %d\n",endchan-startchan+1);
      send_double("fch1",fcent+fabs(pf.hdr.BW)/2.+chbw/2.+(startchan-1)*chbw);
      send_double("foff",-1.0*fabs(chbw));
      send_int("nbits",pf.hdr.nbits);
      send_int("nbeams",1);
      send_int("ibeam",1);
      send_int("nifs",1);
      if (tsamp == 0.0) tsamp=pf.hdr.dt; 
	fprintf(stderr,"Sampling time %f us\n",tsamp*1.0e6);
      send_double("tsamp",tsamp);
      send_double("tstart",pf.hdr.MJD_epoch);
      send_int("telescope_id",32); /* this is going to be the 20 m code */
      send_int("machine_id",32); /* this is going to be Cyborg on 20 m */

      rah=atof(strtok(pf.hdr.ra_str,":"));
      ram=atof(strtok(NULL,":"));
      ras=atof(strtok(NULL,":"));
      src_raj=rah*10000.0+ram*100.0+ras;
      ded=atof(strtok(pf.hdr.dec_str,":"));
      dem=atof(strtok(NULL,":"));
      des=atof(strtok(NULL,":"));
      if (ded<0.0) 
	sgn=-1.0;
      else
	sgn=1.0;
      src_dej=fabs(ded)*10000.0+sgn*dem*100.0+sgn*des;
      az_start=pf.sub.tel_az;
      za_start=pf.sub.tel_zen;
      send_coords(src_raj,src_dej,az_start,za_start);
      send_string("HEADER_END");
      first=0;
      if (pf.hdr.nbits==8) {
	data8=(unsigned char *)malloc(ndumps*sizeof(char)*pf.hdr.nchan);
	npersub=pf.sub.bytes_per_subint;
      } else if (pf.hdr.nbits==16) {
	data16=(unsigned short *)malloc(ndumps*sizeof(short)*pf.hdr.nchan);
	npersub=pf.sub.bytes_per_subint/2;
      } else {
	puts("psrfits2fil currently only works with 8 or 16 bit data");
	exit(0);
      }
    }
    i=j=k=l=0;
    while (i<npersub) {
//	printf("%d %d\n",i,pf.sub.bytes_per_subint);
      if (j==0) {
	k++;
	if ((k>=startchan) && (k<=endchan)) {
	  if (pf.hdr.nbits==8) data8[l++]=pf.sub.data8[i];
	  if (pf.hdr.nbits==16) data16[l++]=pf.sub.data16[i];
	}
	if (k==pf.hdr.nchan) {
	  idump++;
	  if (idump==ndumps) {
	    if (flip) {
	      for (x=l-1;x>=0;x--)
		if (pf.hdr.nbits==8) {
		  fwrite(&data8[x],sizeof(char),1,output);
		} else {
		  fwrite(&data16[x],sizeof(short),1,output);
		}
	    } else {
	      if (pf.hdr.nbits==8) {
		fwrite(data8,sizeof(char),l,output);
	      } else {
		fwrite(data16,sizeof(short),l,output);
	      }
	    }
	    idump=l=0;
	  }
	  k=0;
	}
      }
      i++;
      if (i%pf.hdr.nchan==0) j++;
      if (j==4) j=0;
    }
    if (pf.hdr.nbits==8) free(pf.sub.data8);
    if (pf.hdr.nbits==16) free(pf.sub.data16);
    if ((idump>0)&&(pf.hdr.nbits==8))  fwrite(data8,sizeof(char),l,output);
    if ((idump>0)&&(pf.hdr.nbits==16)) fwrite(data16,sizeof(short),l,output);
    idump=0;
    if (bandpass)  exit(0);
  }
}
