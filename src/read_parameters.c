#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "macros_and_structures.h"
#include "functions.h"



//////////////////////
// GLOBAL VARIABLES //
//////////////////////
/* Struct Global Variables */
extern struct globalVariables GV;



/*
 * Function:  read_parameters
 * --------------------
 * Reads the parameter file in which are the main parameters 
 * necessary to run the code.
 *
 * The information loaded in the same order are: 
 * FILE_NAME:      File name path of the GADGET binary file.
 * OUTPUT:         Path of the output file.
 * DELTA_K:        Width of space sampled for the calculation of 
 *                 the power spectrum. The value is given in terms 
 *                 of the fundamental frequency kF. The value
 *                 should be bigger or equal than 1.
 * 
 *  param_file_name: String with the name of the parameter file.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the parameter file.
 *           -2 --> There is an error whith the settings of the 
 *                  parameter file.
 */
int read_parameters(char param_file_name[])
{
  
  FILE *cfg =NULL;    // Stream to the parameter (config) file
  int   len =LENCHAR; // Len of the read parameter
  char *buf =NULL;    // buf variables to be used to read strings variables
  char *buf1=NULL;
  char *buf2=NULL;
  char *dumb;
  
  
  if( (cfg=fopen(param_file_name,"r"))==NULL )
    {
      printf("%s not found.\n", param_file_name);
      // Value -1 means there is an error loading the param file
      return -1;
    }

  buf  = (char *) malloc( len*sizeof(char) );
  buf1 = (char *) malloc( len*sizeof(char) );
  buf2 = (char *) malloc( len*sizeof(char) );

  /* Reading FILE_NAME parameter */
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'FILE_NAME' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.FILE_NAME = strdup(buf2);
      printf("Reading from File: %s\n", GV.FILE_NAME);
    }

  /* Reading OUTPUT parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'OUTPUT' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.OUTPUT = strdup(buf2);
      printf("Output File: %s\n", GV.OUTPUT);
    }
  
  /* Reading S_KF parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 )
    {
      printf("No 'S_KF' setting in configuration file.\n");
      return -2;
    }
  else
    {
      GV.S_KF=atof(buf2);
      if(GV.S_KF >= 1.0)
	{
	  printf("Binning width in terms of the fundamental frequency kF: %g\n", GV.S_KF);
	}
      else
	{
	  printf("Invalid 'S_KF' setting in configuration file.\n");
	  return -2;
	}
    }
  
  if(dumb==NULL){}
  
  fclose(cfg);
  free(buf);
  free(buf1);
  free(buf2);
  
  return 0;
}
