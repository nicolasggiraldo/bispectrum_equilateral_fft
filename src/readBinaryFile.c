#include <math.h>
#include <stdio.h>
#include <fftw3.h>
#include "macros_and_structures.h"
#include "functions.h"

//////////////////////////////////
// GLOBAL VARIABLES FROM MAIN.C //
//////////////////////////////////
/* Global variables for the FFTW routines */
extern fftw_complex *denConX; // Density contrast in X-space
extern fftw_complex *denConK; // Density contrast in K-space
extern fftw_plan forwardPlan;
/* Struct Global Variables */
extern struct globalVariables GV;

/*
 * Function:  readBinaryFile
 * --------------------                  
 * Reads a binary file with the information of the cell and store 
 * the information in the data structure variable *part* it also 
 * returns the total number of particles.
 * 
 *  There are no arguments in the routiene.              
 *
 *  returns: Integer value.               
 *            0 --> There is no error. 
 *           -1 --> There is an error loading file
 *           -2 --> Structure cell could not be allocated. 
 */
int readBinaryFile()
{
  
  FILE *fdata;
  long id_cell;
  double density_Contrast;
  int n[3]; // Number of grids in each axis for FFT estimation
  size_t err;
  
  printf("\n-----------------------------------------------\n");
  printf("Reading file: %s\n", GV.FILE_NAME);
  
  fdata = fopen(GV.FILE_NAME,"rb");
  if(fdata == NULL)
    {
      printf("File %s cannot be open\n", GV.FILE_NAME);
      return -1;
    }

  /* Getting cosmological parameters of the simulation */
  err = fread(&GV.OMEGA_M0,    sizeof(double), 1, fdata);
  err = fread(&GV.OMEGA_L0,    sizeof(double), 1, fdata);
  err = fread(&GV.ZRS,         sizeof(double), 1, fdata);
  err = fread(&GV.HUBBLEPARAM, sizeof(double), 1, fdata);

  /* Getting simulation parameters */
  err = fread(&GV.NGRID,          sizeof(int),      1, fdata);
  err = fread(&GV.GADGET_VERSION, sizeof(int),      1, fdata);
  err = fread(&GV.L,              sizeof(double),   1, fdata);
  err = fread(&GV.NP_TOT,         sizeof(long int), 1, fdata);
  err = fread(&GV.TOTAL_MASS,     sizeof(double),   1, fdata);
  err = fread(&GV.RHO_MEAN,       sizeof(double),   1, fdata);
  err = fread(&GV.VOL_CELL,       sizeof(double),   1, fdata);
  err = fread(&GV.H,              sizeof(double),   1, fdata);
  err = fread(&(GV.SCHEME[0]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[1]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[2]),    sizeof(char),     1, fdata);
  GV.SCHEME[3] = '\0';

  GV.NGRID2     = (1L*GV.NGRID)*(1L*GV.NGRID);
  GV.NGRID3     = (1L*GV.NGRID)*(1L*GV.NGRID)*(1L*GV.NGRID);
  printf("NGRID2=%d\n",GV.NGRID);
  printf("NGRID2=%ld\n",GV.NGRID2);
  printf("NGRID3=%ld\n",GV.NGRID3);
  GV.SIM_VOL    = GV.L*GV.L*GV.L;
  GV.KF         = (2.0*M_PI) / GV.L;
  GV.DELTA_K    = GV.S_KF * GV.KF;
  GV.SHOT_NOISE = GV.VOL_CELL / GV.NP_TOT;
  GV.KN         = M_PI / GV.H;

  printf("\n-----------------------------------------------\n");
  printf("The original snapshot has a total of %ld particles\n", GV.NP_TOT);
  printf("----------------------------------------\n");
  printf(" * Redshift...     %20.8lf\n", GV.ZRS);
  printf(" * Omega0...       %20.8lf\n", GV.OMEGA_M0);
  printf(" * OmageLa...      %20.8lf\n", GV.OMEGA_L0);
  printf(" * Hubbleparam...  %20.8lf\n", GV.HUBBLEPARAM);
  printf("----------------------------------------\n");
  printf(" * Boxsize...      %20.8lf\n", GV.L);
  printf(" * Ngrid...        %20d\n",    GV.NGRID);
  printf(" * SimMass...      %20.8e\n",  GV.TOTAL_MASS);
  printf(" * Scheme...       %20s\n",    GV.SCHEME);
  printf("----------------------------------------\n");
  printf(" * kF...           %20.8lf\n", GV.KF);
  printf(" * kN...           %20.8lf\n", GV.KN);
  printf(" * DELTA_k...      %20.8lf\n", GV.DELTA_K);
  printf(" * PSshotNoise...  %20.8e\n",  GV.SHOT_NOISE);

  

  /////////////////////////////////
  //* FFTW VARIABLES ALLOCATION *//
  /////////////////////////////////
  denConK = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConK == NULL)
    {
      printf("FFTW arrays could not be allocated\n");
      return -2;
    }//if
  
  denConX = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConX == NULL)
    {
      printf("FFTW arrays could not be allocated\n");
      return -2;
    }//if
  
  
  
  /////////////////////////
  //* ESTABLISHING PLAN *//
  /////////////////////////
  // Number of grids in each axis
  n[X] = n[Y] = n[Z] = GV.NGRID;
  
  // FFTW plan for a 3D Fourier transform
  /* FFTW_MEASURE: find the optimal plan by actually computing several FFTs 
     FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
     FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
     FFTW_IN_PLACE: a plan assumes that the in and out are same */
  forwardPlan = fftw_plan_dft(3, n, denConX, denConK, FFTW_FORWARD, FFTW_ESTIMATE);
  
  
  
  ////////////////////////////////
  //* GETTING DENSITY CONTRAST *//
  ////////////////////////////////
  for(id_cell=0; id_cell<GV.NGRID3; id_cell++)
    {
      err = fread(&density_Contrast, sizeof(double), 1, fdata);
      denConX[id_cell][0] = density_Contrast;
      denConX[id_cell][1] = 0.0;
    }//for id_cell
  fclose(fdata);
  
  
  ////////////////////////////////////
  //* TAKING THE FOURIER TRANSFORM *//
  ////////////////////////////////////
  fftw_execute(forwardPlan);
  printf("\n-----------------------------------------------\n");
  printf("Fourier Transform succes.\n");
  fftw_destroy_plan(forwardPlan);
  fftw_free(denConX);
  
  
  if(err){};
  
  return 0;
}
